// import FMath from 'fmath';

const Tan126 = Math.tan(1.26*1.1);
const sqrt = Math.sqrt;
const atan = Math.atan;
const max = Math.max;
const cos = Math.cos;
const sin = Math.sin;
const exp = Math.exp;
const pow = Math.pow;
const mieG = 0.8;


export function GetTransmittanceToTopAtmosphereBoundary( planetProps, transmittanceGetter, r, mu) {
  let {bottomRadius, topRadius} = planetProps;
  r = clamp(r, bottomRadius, topRadius);
  rangeCheck(r, bottomRadius, topRadius, 'R check to top atm b');
  rangeCheck(mu, -1.0,1.0, 'MuS check');
  let [Width, Height] = getTransmittanceResolution(planetProps);
  let uv = GetTransmittanceTextureUvFromRMu(planetProps, r, mu);
  let pixel = transmittanceGetter(uv);
  if(pixel.length == 0){
    console.log(uv, Width, Height, r, mu);
    throw new Error("empty transmittance")
  }
  nanCheck(pixel);
  return pixel;
}
export function clampRadius(planetProps, r) {
  return clamp(r, planetProps.bottomRadius, planetProps.topRadius);
}

export function GetTransmittanceTextureUvFromRMu(planetProps, r, mu){
  let {topRadius, bottomRadius} = planetProps;
  r = clamp(r, bottomRadius, topRadius);
  rangeCheck(r, bottomRadius, topRadius);
  rangeCheck(mu, -1, 1);
  let [Width, Height] = getTransmittanceResolution(planetProps);

  let H = sqrt(topRadius * topRadius - bottomRadius * bottomRadius);
  // Distance to the horizon.
  let rho = safeSqrt(r * r - bottomRadius * bottomRadius);
  // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
  // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
  let d = DistanceToTopAtmosphereBoundary(planetProps, r, mu);
  let d_min = topRadius - r;
  let d_max = rho + H;
  let x_mu = (d - d_min) / (d_max - d_min);
  let x_r = rho / H;
  let u = GetTextureCoordFromUnitRange(x_mu, Width);
  let v = GetTextureCoordFromUnitRange(x_r, Height)
  //if(u <0 || u > 1 || v <0 || v > 1){
    //console.log('ALARM', x_mu, Width, x_r, Height, d, H, rho, r, mu);
    //throw new Error("Incorrect uvs");
  //}
  return vec2(clamp(u,0,1) , clamp(v, 0, 1));
}

function vec2(x,y){
  return [x,y];
}

export function DistanceToBottomAtmosphereBoundary(planetProps, r, mu) {
  let {topRadius, bottomRadius} = planetProps;
  if(r < bottomRadius) throw new Error("incorrect r")
  let discriminant = r * r * (mu * mu - 1.0) +
      bottomRadius * bottomRadius;
  let distance = clampDistance(-r * mu - safeSqrt(discriminant));
  if(isNaN(distance))
     throw `Nan results distance:${distance}`;
  return distance;
}
export function DistanceToTopAtmosphereBoundary(planetProps, r, mu) {
  let {topRadius, bottomRadius} = planetProps;
  if(r > topRadius) throw new Error("incorrect r")
  let discriminant = r * r * (mu * mu - 1.0) + topRadius * topRadius;
  let distance = clampDistance(-r * mu + safeSqrt(discriminant));
  if(isNaN(distance))
     throw `Nan results distance:${distance}`;
  return distance
}

export function GetRMuFromTransmittanceTextureUv(planetProps, u, v){
  assert(u >= 0 && u <= 1.0, 'Incorrect uv range')
  assert(v >= 0 && v <= 1.0, 'Incorrect uv range')
  let {topRadius, bottomRadius} = planetProps;

  let [Width,Height] = getTransmittanceResolution(planetProps);
  let xMu = GetUnitRangeFromTextureCoord(u,Width);
  let xR  = GetUnitRangeFromTextureCoord(v,Height);

  let H = sqrt(topRadius * topRadius - bottomRadius * bottomRadius);
  // Distance to the horizon, from which we can compute r:
  let rho = H * xR;
  let r = sqrt(rho * rho + bottomRadius * bottomRadius);
  //if(r > topRadius){
    //console.log(r, topRadius, rho, xR, H);
  //}
  // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
  // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon) -
  // from which we can recover mu:
  let d_min = topRadius - r;
  let d_max = rho + H;
  let d = d_min + xMu * (d_max - d_min);
  let mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
  mu = clampCosine(mu);
  if(isNaN(r) || isNaN(mu))
     throw `Nan results r:${r}, mu:${mu}`;
  return {mu,r}
}

export function nanCheck(arr, fn){
  for(let i=0; i< arr.length; ++i){
    if(isNaN(arr[i])) {
      if(!fn) throw new Error(`nan value in array ${serializeArray(arr)}`);
      fn();
    }
  }
}

function serializeArray(arr){
  return JSON.stringify(arr);
}
export function rangeCheck(x, m, M, message="Range check warning"){
  if(x < m || x > M)
    console.warn(`${message}: ${m} < ${x} > ${M}`);
}

export function getLayerDensity(layer, altitude) {
  let density = layer.exp_term * exp(layer.exp_scale * altitude) +
      layer.linear_term * altitude + layer.constant_term;
  return clamp(density, Number(0.0), Number(1.0));
}

export function getProfileDensity(profile, altitude) {
  return altitude < profile[0].width ?
      getLayerDensity(profile[0], altitude) :
      getLayerDensity(profile[1], altitude);
}
export function clampDistance(d) {
  return max(d, 0.0);
}
export function safeSqrt(a) {
  return sqrt(max(a, 0.0));
}

export function clampCosine(mu){
  return clamp(mu, -1.0, 1.0);
}

export function assert(cond, message = "Assertion failed") {
  if(!cond) throw new Error(message);
}

export function getIrradianceResolution(planetProps){
  let {resMus, resR} = planetProps;
  return [resMus*2, resR/2];
}

export function getTransmittanceResolution(planetProps){
  let {resMu, resR} = planetProps;
  return [resMu*2, resR*2];

}

export function GetTextureCoordFromUnitRange( x,  texture_size) {
  return 0.5 / texture_size + x * (1.0 - 1.0 / texture_size);
}

export function GetUnitRangeFromTextureCoord(u, textureSize){
  return (u - 0.5 / textureSize) / (1.0 - 1.0 / textureSize);
}

export function phaseFunctionRay(mu) {
  return (3.0 / (16.0 * Math.PI)) * (1.0 + mu * mu);
}

// Mie phase function
export function phaseFunctionMie(mu) {
  return 1.5 * 1.0 / (4.0 * Math.PI) * (1.0 - mieG*mieG) * pow(1.0 + (mieG*mieG)
  - 2.0*mieG*mu, -3.0/2.0) * (1.0 + mu * mu) / (2.0 + mieG*mieG);
}

export function vNan(v){
  for(let c =0;c <v.length; ++c){
    if(isNaN(v[c])) return true;
  }
  return false;
}

export function mulS(v, s){
  return [v[0]*s, v[1]*s, v[2] * s];
}

export function clamp(v, m, M){
  return Math.min(Math.max(v,m), M);
}

export function dot(v1, v2){
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];


}
export function vadd(v1,v2){
  let r = new Float32Array(3);
  r[0] = v1[0] + v2[0],
  r[1] = v1[1] + v2[1],
  r[2] = v1[2] + v2[2]
  return r;
}

export function vmul4(to, v1,v2){
  to[0] = v1[0] * v2[0];
  to[1] = v1[1] * v2[1];
  to[2] = v1[2] * v2[2];
  to[3] = v1[3] * v2[3];

}
export function vmul(to, v1,v2){
  to[0] = v1[0] * v2[0];
  to[1] = v1[1] * v2[1];
  to[2] = v1[2] * v2[2];
}
export function vmul1( v1,v2){
  let to = new Float32Array(3);
  to[0] = v1[0] * v2[0];
  to[1] = v1[1] * v2[1];
  to[2] = v1[2] * v2[2];
  return to;
}
export function vdiv(v1,v2){
  let r = new Float32Array(3);
  r[0] = v1[0] == 0.0? 0: v1[0] / v2[0],
  r[1] = v1[1] == 0.0? 0: v1[1] / v2[1],
  r[2] = v1[2] == 0.0? 0: v1[2] / v2[2]
  return r;
}

export function vmin(v1,v2){
  return [
    Math.min(v1[0], v2[0]),
    Math.min(v1[1], v2[1]),
    Math.min(v1[2], v2[2])
  ]
}

const TAN15= Math.tan(1.5);
const K = TAN15 / 1.15;

let STORE = new Map;

/*
vec4 texture4D(sampler3D table, float r, float mu, float muS, float nu)
{
    float H = sqrt(Rt * Rt - Rg * Rg);
    float rho = sqrt(r * r - Rg * Rg);
#ifdef INSCATTER_NON_LINEAR
    float rmu = r * mu;
    float delta = rmu * rmu - r * r + Rg * Rg;
    vec4 cst = rmu < 0.0 && delta > 0.0 ? vec4(1.0, 0.0, 0.0, 0.5 - 0.5 / float(RES_MU)) : vec4(-1.0, H * H, H, 0.5 + 0.5 / float(RES_MU));
	float uR = 0.5 / float(RES_R) + rho / H * (1.0 - 1.0 / float(RES_R));
    float uMu = cst.w + (rmu * cst.x + sqrt(delta + cst.y)) / (rho + cst.z) * (0.5 - 1.0 / float(RES_MU));
    // paper formula
    //float uMuS = 0.5 / float(RES_MU_S) + max((1.0 - exp(-3.0 * muS - 0.6)) / (1.0 - exp(-3.6)), 0.0) * (1.0 - 1.0 / float(RES_MU_S));
    // better formula
    float uMuS = 0.5 / float(RES_MU_S) + (atan(max(muS, -0.1975) * tan(1.26 * 1.1)) / 1.1 + (1.0 - 0.26)) * 0.5 * (1.0 - 1.0 / float(RES_MU_S));
#else
	float uR = 0.5 / float(RES_R) + rho / H * (1.0 - 1.0 / float(RES_R));
    float uMu = 0.5 / float(RES_MU) + (mu + 1.0) / 2.0 * (1.0 - 1.0 / float(RES_MU));
    float uMuS = 0.5 / float(RES_MU_S) + max(muS + 0.2, 0.0) / 1.2 * (1.0 - 1.0 / float(RES_MU_S));
#endif
    float lerp = (nu + 1.0) / 2.0 * (float(RES_NU) - 1.0);
    float uNu = floor(lerp);
    lerp = lerp - uNu;
    return texture3D(table, vec3((uNu + uMuS) / float(RES_NU), uMu, uR)) * (1.0 - lerp) +
           texture3D(table, vec3((uNu + uMuS + 1.0) / float(RES_NU), uMu, uR)) * lerp;
}

*/
export function tableLookup(getter, planetProps){
  let {radius, atmosphereHeight} = planetProps.phisical;
  let {resR, resMu, resNu, resMus} = planetProps;
  let Rg = radius;
  let Rt = radius + atmosphereHeight;
  return ([r, mu, muS, nu])=>{

    let H = sqrt(Rt* Rt - Rg*Rg);
    let rho = 0;
    if(r >= Rg) rho = sqrt(r*r - Rg*Rg);

    let rmu = r * mu;
    let delta = rmu * rmu - r * r + Rg * Rg;
    let cst = rmu < 0.0 && delta > 0.0 
      ? [1.0, 0.0, 0.0, 0.5 - 0.5 / resMu] 
      : [-1.0, H * H, H, 0.5 + 0.5 / resMu];

    let uR = 0.5 / resR + rho / H * (1 - 1 / resR);

    let uMu_ = cst[3] + (rmu * cst[0] + sqrt(delta + cst[1])) / (rho + cst[2]) * (0.5 - 1.0 / resMu);
    let uMu = clamp(uMu_, 0.0, 1.0);
    if(isNaN(uMu)) uMu = 0.0;
    // paper formula
    // float uMuS = 0.5 / float(RES_MU_S) + max((1.0 - exp(-3.0 * muS - 0.6)) / (1.0 - exp(-3.6)), 0.0) * (1.0 - 1.0 / float(RES_MU_S));
    // better formula
    let uMus = 0.5 / resMus 
      + (atan(max(muS, -0.1975) * Tan126) / 1.1 + (1.0 - 0.26)) * 0.5 * (1.0 - 1.0 / resMus);

     
    let lerp = (nu + 1.0) / 2.0 * (resNu - 1.0);
    let uNu = Math.max((lerp ) / resNu, 0.0);
    if(uMus < 0 || uNu < 0 || uMu < 0 || uR < 0 || uMus > 1 || uNu > 1 || uMu > 1 || uR > 1 ||isNaN(uMus) || isNaN(uNu) || isNaN(uMu) || isNaN(uR)){
      console.log([uMus, uNu, uMu, uR]);
      if(isNaN(uMu)){
        console.log('umu', uMu_, cst, rmu, delta, rho, resMu);
      }
      debugger;
    }

    return getter([uMus, uNu, uMu, uR]);

  }
}

export function texture2DGetter(texture, dimensions, components){
  let b = new Float32Array(2);
  return uv=>{
    b[0] = Math.floor((dimensions[0]-1) * uv[0]);
    b[1] = Math.floor((dimensions[1]-1) * uv[1]);
    let ix = components*(b[1] * dimensions[0] + b[0]);
    if(ix < 0 || ix > texture.length || ix+components > texture.length)
      console.log(ix, b, dimensions, components, texture.length);
    if(!texture.subarray)
      console.log(texture);
    return texture.subarray(ix, ix+components);
  }
}

export function texture4DGetter(texture, dimensions, components){
  let b = new Float32Array(4);
  return coords=>{
    vmul4(b, coords, dimensions);
    let i = Math.floor(b[0]);
    let j = Math.floor(b[1]);
    let k = Math.floor(b[2]);
    let r = Math.floor(b[3]);
    let X = i*dimensions[1] + j;
    let XX = X * dimensions[2]  + k
    let ix =  components * (XX * dimensions[3] + r);
    let color = texture.subarray(ix, ix+components);
    if(color.length < components)
      debugger;

    return color;
  }

}
export class Irradiance{
  constructor(texture, radius, atmosphereHeight, resMus, resR){
    this.texture = texture;
    this.resMus = resMus;
    this.resR = resR;
    this.radius = radius;
    this.atmosphereHeight = atmosphereHeight;
  }

  getIrradianceIx(uv){
    let W = this.resMus / 2;
    let H = this.resR * 2;
    let c0 = Math.min(Math.abs(Math.floor(uv[0]) * H), H-1);
    let c1 = Math.min(Math.abs(Math.floor(uv[1]) * W), W-1);
    return 4 * (c0 * W + c1);
  }

  getIrradianceUV(r, muS){
    let ur = (r -this.radius) / this.atmosphereHeight;
    let umus = (muS + 0.2) / 1.2
    return [umus, ur];
  }

  getIrradiance(r, mu) {
    let uv = this.getIrradianceUV(r,mu);
    let ix = this.getIrradianceIx(uv);
    let color = this.texture.subarray(ix, ix+3);
    return color;

  }
}

export class Transmittance{
  constructor(texture, radius, atmosphereHeight, resMu, resR, linear=false){
    this.texture = texture;
    this.atmosphereHeight = atmosphereHeight;
    this.radius = radius;
    this.linear = linear;
    this.resMu = resMu;
    this.resR = resR;
  }

  getTransmittance(arr){
    let v = null
    if(arr.length == 2) v = this.getTransmittanceRMU(arr);
    if(arr.length == 3) v = this.getTransmittanceRMUD(arr);
    return v;
  }

  getTransmittanceIx(uv){
    let W = this.resMu * 2;
    let H = this.resR * 2;
    let i = Math.min(Math.abs(Math.floor(uv[0] * W)), W-1);
    let j = Math.min(Math.abs(Math.floor(uv[1] * H)), H-1);
    return 4 * (j * W + i);

  }

  getTransmittanceRMU([r,mu]){
    let __t = Date.now();
    let uv = this.getTransmittenceUV(r,mu, this.radius, this.atmosphereHeight);
    let ix = this.getTransmittanceIx(uv);
    let color = this.texture.subarray(ix, ix+3);
    return color;
  }

  getTransmittanceRMUD([r, mu, d]){
    let _rSQR = r*r;
    let r1 = Math.sqrt(_rSQR + d*d  + 2*r*mu*d);
    let mu1 = (r*mu + d) / r1;
    if(mu > 0){
      let t1 = this.getTransmittanceRMU([r, mu]);
      let t2 = this.getTransmittanceRMU([r1,mu1]);
      return vmin( vdiv(t1, t2), [1,1,1]);
    }else{
      let t1 = this.getTransmittanceRMU([r1, -mu1]);
      let t2 = this.getTransmittanceRMU([r,-mu]);
      return vmin( vdiv(t1, t2), [1,1,1]);
    }
  }
  getTransmittenceUV(r, mu){
    let uR = 0, uMu;
    let Rg = this.radius;
    let Rt = Rg + this.atmosphereHeight;
    if(!this.linear){
      
      if(r >= Rg) 
        uR = Math.sqrt((r - Rg) / this.atmosphereHeight);
      uMu = Math.atan((mu + 0.15) * K) / 1.5;
    }else{
      if(r >= Rg) 
        uR = (r - Rg) / (Rt - Rg);
      uMu = (mu + 0.15) / (1.15);
    }
    return [uMu, uR];
  }
}


export function limit(r, mu, planetProperties) {
  let {atmosphereHeight, radius} = planetProperties.phisical;
  let RL = radius + atmosphereHeight + 1;
  let Rg = radius;

  let dout = -r * mu + Math.sqrt(r * r * (mu * mu - 1.0) + RL * RL);

  let delta2 = r * r * (mu * mu - 1.0) + Rg * Rg;
  if (delta2 >= 0.0) {
    let din = -r * mu - Math.sqrt(delta2);
    if (din >= 0.0) {
      dout = Math.min(dout, din);
    }
  }
  return dout;
}

export function precalculateR(planetProps) {
  let {radius, atmosphereHeight, HR, betaR} = planetProps.phisical;
  let {resMu, resNu, resMus, resR} = planetProps;
  return (k, count) => {
    let Rt = radius + atmosphereHeight;
    let _radiusSQR = radius * radius;
    let dr = (k == 0)
      ?0.01
      :(k==(count-1))
        ?-0.001
        :0.0;
    
    let w = k / (count-1);
    let r = Math.sqrt(_radiusSQR + w*w*(Rt*Rt - _radiusSQR)) + dr;
    let _rSQR = r*r;
    let dmin = Rt - r;
    let dmax = Math.sqrt(_rSQR - _radiusSQR) + Math.sqrt(Rt*Rt - _radiusSQR);
    let dminp = r-radius;
    let dmaxp = Math.sqrt(_rSQR - _radiusSQR);
    return {w, r, dmin, dmax, dminp, dmaxp, _rSQR, _radiusSQR}
  }
}

export function calculateMu(planetProps){
  let {radius, atmosphereHeight, HR, betaR} = planetProps.phisical;
  let {resMu, resNu, resMus, resR} = planetProps;
  let Rt = radius + atmosphereHeight;

  return (k, depth, precalculations) => {
    let {dmin, dmax, dminp, dmaxp, r} = precalculations;
    let y = k - 0.5;
    let mu;
    let _q = resMu / 2 - 1.0
    if(y < resMu / 2){
      let d = 1 - y / _q;
      d = Math.min(Math.max(dminp, d*dmaxp), dmaxp * 0.999);
      mu = (radius*radius - r*r - d*d) / (2*r*d);
      mu = Math.min(mu, -Math.sqrt(1.0 - Math.pow(radius/r, 2)) - 0.001);
    }else {
      let d = (y - resMu/2) / _q;
      d = Math.min(Math.max(dmin, d*dmax), dmax * 0.999);
      mu = (Rt*Rt - r*r - d*d) / (2*r*d);
    }
    return {mu};
  }
}

export function calculateMuS(planetProps){
  let {resMus} = planetProps;
  return (i, width, pre)=>{
    let x = i - 0.5;
    let muS = x / (resMus-1);
    muS = Math.tan((2 * muS -1 + 0.26)*1.1)/Tan126;
    return {muS}
  }
}
export function calculateNu(planetProps){
  let {resNu} = planetProps;
  return (k, depth, pre)=>{
    let z = k -0.5;
    let nu = -1.0 + z / (resNu - 1) * 2;
    return {nu};
  }
}

export const precalcs = {
  count: precalculateR,
  width: calculateMuS,
  height: calculateNu,
  depth: calculateMu
}
