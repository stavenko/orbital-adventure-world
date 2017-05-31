import {
  texture2DGetter,
  clampCosine,
  clampRadius,
  getTransmittanceResolution,
  nanCheck,
  rangeCheck,
  vmin,
  clamp,
  GetUnitRangeFromTextureCoord,
  GetTransmittanceToTopAtmosphereBoundary,
  DistanceToBottomAtmosphereBoundary,
  DistanceToTopAtmosphereBoundary,
  GetTransmittance,
  getProfileDensity,
  RayIntersectsGround,
  GetRMuMuSNuFromScatteringTextureFragCoord, vmul1, vdiv,
  assert,
  DistanceToNearestAtmosphereBoundary
} from './utils.js'

const components = 4;
const sqrt = Math.sqrt;
const min = Math.min;
let COUNTERS = {};
export function computeSingleScattering(planetProps, transmittanceTexture){
  let transmittanceRes = getTransmittanceResolution(planetProps);
  let transmittanceGetter = texture2DGetter(transmittanceTexture, transmittanceRes, components);
  let {resMu, resNu, resR, resMus} = planetProps

  let arraySize = resNu * resMu * resMus * resR * components
  let deltaRayleigh = new Float32Array(arraySize);
  let deltaMie = new Float32Array(arraySize);
  let scatteringTexture = new Float32Array(arraySize);
  let counters = new Int32Array(4);
  let sizes = new Int32Array([resMus, resNu, resMu, resR]);
  for(;counters[0] < sizes[0]; ++counters[0]){
    console.log(counters);
    for(counters[1] = 0;counters[1] < sizes[1]; ++counters[1]){
      for(counters[2] = 0;counters[2] < sizes[2]; ++counters[2]){
        for(counters[3] = 0;counters[3] < sizes[3]; ++counters[3]){
          let ix = components * index4(counters, sizes);
          let {delta_rayleigh, delta_mie} = ComputeSingleScatteringTexture(
          planetProps, transmittanceGetter, counters, sizes);
          let scattering = [...delta_rayleigh.slice(0, 3), delta_mie[0]];
          nanCheck(delta_rayleigh);
          nanCheck(delta_mie);
          nanCheck(scattering);
          for(let i=0; i<components;++i){
            deltaRayleigh[ix + i] = delta_rayleigh[i];
            deltaMie[ix + i] = delta_mie[i];
            scatteringTexture[ix + i] = scattering[i];
          }
        }
      }
    }
  }
  let S = 0;
  for(let i =0; i< scatteringTexture.length; ++i) S+= scatteringTexture[i];
  if(S < 1) throw new Error("Total sum of scattering texture is less than 1", S);
  return {deltaRayleigh, deltaMie, scatteringTexture};
}

function index4(counters, sizes){
  return ((counters[0] * sizes[1] + counters[1])*sizes[2] + counters[2])*sizes[3] + counters[3];
}

function ComputeSingleScatteringTexture(planetProps, transmittanceGetter, 
                                        counters, sizes){

  let coords = GetRMuMuSNuFromScatteringTextureFragCoord(planetProps, counters, sizes);

  return ComputeSingleScattering(planetProps, transmittanceGetter, coords);
}



function count_(val){
  if(!COUNTERS[val] ) COUNTERS[val] = 1;
  else COUNTERS[val]++;
}

function ComputeSingleScattering( planetProps, transmittanceGetter, coords){
  let {
    topRadius, bottomRadius, 
    solarIrradiance, rayleighScattering, mieScattering
  } = planetProps;

  //rayleighScattering = [...rayleighScattering, 0];
  //mieScattering = [...mieScattering, 0];
  //olarIrradiance = [..., 0];
  let {nu, mu_s, mu, r, ray_r_mu_intersects_ground} = coords;

  rangeCheck(r, bottomRadius, topRadius);
  rangeCheck(mu, -1, 1);
  rangeCheck(mu_s, -1, 1);
  rangeCheck(nu, -1, 1);

  // Number of intervals for the numerical integration.
  const SAMPLE_COUNT = 50;
  // The integration step, i.e. the length of each integration interval.
  let dx = DistanceToNearestAtmosphereBoundary(planetProps, r, mu,
          ray_r_mu_intersects_ground) / SAMPLE_COUNT;
  // Integration loop.
  let rayleigh_sum = [0,0,0,0];
  let mie_sum = [0,0,0,0];

  for (let i = 0; i <= SAMPLE_COUNT; ++i) {
    let d_i = i * dx;
    // The Rayleigh and Mie single scattering at the current sample point.
    //let rayleigh_i;
    //let mie_i;
    let {rayleigh_i, mie_i} = ComputeSingleScatteringIntegrand(planetProps, transmittanceGetter, r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground);
    // Sample weight (from the trapezoidal rule).
    let weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
    for(let j = 0; j < 4; ++j){
      rayleigh_sum[j] += rayleigh_i[j] * weight_i;
      mie_sum[j] += mie_i[j] * weight_i;
    }
  }
  let rayleigh = [];
  let mie = [];

  for(let i = 0; i<3; ++i){
    rayleigh[i] = rayleigh_sum[i] * dx * solarIrradiance[i] * rayleighScattering[i];
    mie[i] = mie_sum[i] * dx * solarIrradiance[i] * mieScattering[i];
  }
  rayleigh[3] = 0;
  mie[3] = 0;
  let S = rayleigh[0] + rayleigh[1] + rayleigh[2];
  if(S == 0) count_("zero rayleigh");
  else count_("non-zero rayleigh");

  nanCheck(rayleigh,()=>{
    console.log(rayleigh, rayleigh_sum, dx, solarIrradiance, rayleighScattering)
    throw new Error("Err");
  });
  nanCheck(mie);
  return {delta_rayleigh: rayleigh, delta_mie: mie};
}


function ComputeSingleScatteringIntegrand( planetProps, transmittanceGetter, r, mu, mu_s, nu, d, ray_r_mu_intersects_ground){

  let {topRadius, bottomRadius, muSMin, rayleighDensity, mieDensity} = planetProps;
  let r_d = clampRadius(planetProps, sqrt(d * d + 2.0 * r * mu * d + r * r));
  let mu_s_d = clampCosine((r * mu_s + d * nu) / r_d);
  let rayleigh = [0,0,0,0], mie=[0,0,0,0];

  if (RayIntersectsGround(planetProps, r_d, mu_s_d)) {
    rayleigh = [0,0,0,0];
    mie = [0,0,0,0];
  } else {
    let tr1 = GetTransmittance(planetProps, transmittanceGetter, r, mu, d,
            ray_r_mu_intersects_ground);
    let tr2 = GetTransmittanceToTopAtmosphereBoundary(
            planetProps, transmittanceGetter, r_d, mu_s_d)
    let transmittance = vmul1(tr1, tr2);
    nanCheck(transmittance, ()=>{
      console.log(transmittance, tr1, tr2);
      throw new Error('err tr');
    });

    let profileRay = getProfileDensity(
        rayleighDensity, r_d - bottomRadius)
    let profileMie =getProfileDensity(
        mieDensity, r_d - bottomRadius)
    for(let i =0 ; i < 3; ++i){
      rayleigh[i] = transmittance[i] * profileRay;
      mie[i] = transmittance[i] * profileMie;

    }
    nanCheck(rayleigh, ()=>{
      console.log(rayleigh, transmittance, profileRay);
      throw new Error('err rayleigh');
    });
  }
  nanCheck(mie);
  return {rayleigh_i: rayleigh, mie_i: mie};
}



