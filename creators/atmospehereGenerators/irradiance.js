import {mulS, dot, Transmittance,
  tableLookup,
  texture4DGetter,
  texture2DGetter,
  getIrradianceResolution,
  getTransmittanceResolution,
  GetUnitRangeFromTextureCoord,  
  clampCosine,
  phaseFunctionMie, 
  phaseFunctionRay,
  vadd,
  vmul1,
  GetTransmittanceToTopAtmosphereBoundary,
  nanCheck,
  rangeCheck,
  ComputeIndirectIrradiance,
  assert,
  clamp
} from './utils.js';

const sqrt = Math.sqrt;
const max = Math.max;
const cos = Math.cos;
const sin = Math.sin;
function vec3(x,y,z){
  return new Float32Array([x,y,z]);
}

export function getDeltaEIterativeColor(planetProps, deltaSRTexture, deltaSMTexture){
  let {radius, atmosphereHeight} = planetProps.phisical;
  let {resMu, resR, resMus, resNu}  = planetProps;
  const dphi = Math.PI / planetProps.IrradianceIntegralSamples;
  const dtheta = Math.PI / planetProps.IrradianceIntegralSamples;

  const deltaSR = tableLookup(texture4DGetter(deltaSRTexture, [resMus, resNu, resMu, resR],4), planetProps);
  const deltaSM = tableLookup(texture4DGetter(deltaSMTexture, [resMus, resNu, resMu, resR],4), planetProps);
  return ({i,j,W,H}, iteration)=>{
    let [r, muS] = getIrradianceRMuS(i, j, W, H);
    let s = vec3(max(sqrt(1.0 - muS * muS), 0.0), 0.0, muS);

    let result = vec3(0, 0, 0);
    // integral over 2.PI around x with two nested loops over w directions (theta,phi) -- Eq (15)
    for (let iphi = 0; iphi < 2 * planetProps.IrradianceIntegralSamples; ++iphi) {
        let phi = (iphi + 0.5) * dphi;
        for (let itheta = 0; itheta < planetProps.IrradianceIntegralSamples/ 2; ++itheta) {
            let theta = (itheta + 0.5) * dtheta;
            let dw = dtheta * dphi * sin(theta);
            let w = vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
            let nu = dot(s, w);
            let dRayMie;
            if (iteration == 0) {
                // first iteration is special because Rayleigh and Mie were stored separately,
                // without the phase functions factors; they must be reintroduced here
                let pr1 = phaseFunctionRay(nu);
                let pm1 = phaseFunctionMie(nu);
                let ray = deltaSR([r, w[2], muS, nu]);
                let mie = deltaSM([r, w[2], muS, nu]);

                ray = mulS(ray, pr1);
                mie = mulS(mie, pm1);
                dRayMie = vadd(ray, mie);

            } else {
              dRayMie = deltaSR([r, w[2], muS, nu]);
            }
            dRayMie = mulS(dRayMie, w[2] * dw);
            result[0] += dRayMie[0];
            result[1] += dRayMie[1];
            result[2] += dRayMie[2];
        }
    }
    if(isNaN(result[0]) || isNaN(result[1]) || isNaN(result[2]))
       debugger;
    return [...result,0];
  }


  function getIrradianceRMuS(i,j,W,H){
    let r = radius + (j - 0.5) / (H-1) * atmosphereHeight;
    let muS = -0.2 + (i - 0.5) / (W-1) * (1.0 + 0.2);
    return [r,muS]
  }
}


function ComputeDirectIrradiance(planetProps, transmittanceGetter, r, mu_s) {
  let {bottomRadius, topRadius} = planetProps;
  r = clamp(r, bottomRadius, topRadius);
  rangeCheck(r, bottomRadius, topRadius, 'R check');
  rangeCheck(mu_s,-1.0,1.0, 'MuS check');
  let {solarIrradiance} = planetProps;
  let tr = GetTransmittanceToTopAtmosphereBoundary(planetProps, transmittanceGetter, r, mu_s);
  let result = vmul1(solarIrradiance, tr);
  result = mulS(result, max(mu_s, 0.0));
  nanCheck(result,()=>{
    console.log(r, mu_s, tr, solarIrradiance);
    throw new Error("nan or null values");
  });

  return result;
}

function GetRMuSFromIrradianceTextureUv(planetProps, u, v){
  let {topRadius, bottomRadius} = planetProps;
  assert(u >= 0.0 && u <= 1.0);
  assert(v >= 0.0 && v <= 1.0);
  let [Width, Height] = getIrradianceResolution(planetProps);
  let x_mu_s = GetUnitRangeFromTextureCoord(u, Width);
  let x_r = GetUnitRangeFromTextureCoord(v, Height);
  let r = bottomRadius + x_r * (topRadius - bottomRadius);
  let mu_s = clampCosine(2.0 * x_mu_s - 1.0);
  return {mu_s, r};
}

export function getDeltaIrradiance({s,t}, planetProps, transmittanceTexture){

  let res = getTransmittanceResolution(planetProps);
  let transmittanceGetter = texture2DGetter(transmittanceTexture,res, 4);
  let {mu_s, r} = GetRMuSFromIrradianceTextureUv( planetProps, s,t);
  return ComputeDirectIrradiance(planetProps, transmittanceGetter, r, mu_s);
}

function isNotZero(v){
  for(let c = 0; c < v.length; ++c)
    if(v[c] ==0) return false;
  return true;
}
const components = 4;

export function ComputeIndirectIrradianceTexture(planetProps, irradianceTexture, deltaIrradianceTexture, singleRayleighGetter, singleMieGetter, multipleScatteringGetter, scatteringOrder){
  let [Width, Height] = getIrradianceResolution(planetProps);
  for(let i = 0; i< Width; ++i){
    for(let j = 0; j < Height; ++j){
      let s = i / (Width-1);
      let t = j / (Height-1);
      let {mu_s, r} = GetRMuSFromIrradianceTextureUv(planetProps, s,t);
      let irradiance = ComputeIndirectIrradiance(planetProps, singleRayleighGetter, singleMieGetter, multipleScatteringGetter, r, mu_s, scatteringOrder);
      let ix = components*(j * Width + i);

      nanCheck(irradiance);
      for(let c =0; c < 3; ++c){
        irradianceTexture[ix +c] += irradiance[c];
        deltaIrradianceTexture[ix+c] = irradiance[c];
      }
    }
  }
}
