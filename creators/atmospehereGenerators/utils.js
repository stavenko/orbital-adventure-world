// import FMath from 'fmath';

const Tan126 = Math.tan(1.26*1.1);
const sqrt = Math.sqrt;
const atan = Math.atan;
const max = Math.max;
const min = Math.min;
const cos = Math.cos;
const sin = Math.sin;
const exp = Math.exp;
const pow = Math.pow;
const PI = Math.PI;
const mieG = 0.8;

const LerpIncrements2 = [
  new Int32Array([0,0]),
  new Int32Array([0,1]),
  new Int32Array([1,0]),
  new Int32Array([1,1]),
];
const LerpIncrements3 = [
  new Int32Array([0,0,0]),
  new Int32Array([0,0,1]),
  new Int32Array([0,1,0]),
  new Int32Array([0,1,1]),
  new Int32Array([1,0,0]),
  new Int32Array([1,0,1]),
  new Int32Array([1,1,0]),
  new Int32Array([1,1,1]),
];
const LerpIncrements4 = [
  new Int32Array([0,0,0,0]),
  new Int32Array([0,0,0,1]),
  new Int32Array([0,0,1,0]),
  new Int32Array([0,0,1,1]),
  new Int32Array([0,1,0,0]),
  new Int32Array([0,1,0,1]),
  new Int32Array([0,1,1,0]),
  new Int32Array([0,1,1,1]),
  new Int32Array([1,0,0,0]),
  new Int32Array([1,0,0,1]),
  new Int32Array([1,0,1,0]),
  new Int32Array([1,0,1,1]),
  new Int32Array([1,1,0,0]),
  new Int32Array([1,1,0,1]),
  new Int32Array([1,1,1,0]),
  new Int32Array([1,1,1,1]),
];
let LerpIncrements = [[],
  [[0], [1]],
  LerpIncrements2,
  LerpIncrements3,
  LerpIncrements4
]

export function normalize(v){
  let u = new Float32Array(v.length);
  let l2 = 0;
  for(let i=0; i< v.length; ++i) l2 += v[i]*v[i];
  let l = sqrt(l2);
  for(let i =0; i<v.length; ++i) u[i] = v[i] / l
  return u;
}

export function index4(counters, sizes){
  return ((counters[0] * sizes[1] + counters[1])*sizes[2] + counters[2])*sizes[3] + counters[3];
}

function vmul(t, a, b){
  for(let i =0; i< t.length; ++i) t[i] = a[i] + b[i];
}

export function ComputeMultipleScattering( planetProps, transmittanceGetter, scatteringDensityGetter, counters, sizes){
  //Length r;
  //Number mu;
  //Number mu_s;
  //bool ray_r_mu_intersects_ground;
  let {r,mu, mu_s, nu, ray_r_mu_intersects_ground} =GetRMuMuSNuFromScatteringTextureFragCoord(planetProps, counters, sizes);
  return ComputeMultipleScattering(atmosphere, transmittance_texture,
      scattering_density_texture, r, mu, mu_s, nu,
      ray_r_mu_intersects_ground);
}
export function DistanceToNearestAtmosphereBoundary(planetProps, r, mu, ray_r_mu_intersects_ground) {
  if (ray_r_mu_intersects_ground) {
    return DistanceToBottomAtmosphereBoundary(planetProps, r, mu);
  } else {
    return DistanceToTopAtmosphereBoundary(planetProps, r, mu);
  }
}

function ComputeMultipleScattering( planetProps, transmittanceGetter, scatteringDensityGetter, r, mu, mu_s, nu, ray_r_mu_intersects_ground) {

  // Number of intervals for the numerical integration.
  const SAMPLE_COUNT = 50;
  // The integration step, i.e. the length of each integration interval.
  let dx = DistanceToNearestAtmosphereBoundary( planetProps, r, mu, ray_r_mu_intersects_ground) / SAMPLE_COUNT;
  // Integration loop.
  let rayleigh_mie_sum = [0,0,0,0]
  for (let i = 0; i <= SAMPLE_COUNT; ++i) {
    let d_i = i * dx;

    // The r, mu and mu_s parameters at the current integration point (see the
    // single scattering section for a detailed explanation).
    let r_i = clampRadius(planetProps, sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r));
    let mu_i = clampCosine((r * mu + d_i) / r_i);
    let mu_s_i = clampCosine((r * mu_s + d_i * nu) / r_i);

    // The Rayleigh and Mie multiple scattering at the current sample point.
    let scat = GetScattering_(
      planetProps,
      scatteringDensityGetter,
      r_i, mu_i, mu_s_i, nu, ray_r_mu_intersects_ground)
    let tr = GetTransmittance(planetProps, transmittanceGetter, r, mu, d_i, ray_r_mu_intersects_ground);

    let rayleigh_mie_i = [0,0,0,0];
    for(let cc = 0; cc < scat.length; ++cc)
      rayleigh_mie_i[cc] = scat[cc] * tr[cc] * dx;
    // Sample weight (from the trapezoidal rule).
    let weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
    for(let cc = 0; cc < rayleigh_mie_sum.length; ++cc)
    rayleigh_mie_sum[cc] += rayleigh_mie_i[cc] * weight_i;
  }
  return rayleigh_mie_sum;
}

function GetScattering_(planetProps, textureGetter, r,  mu,  mu_s,  nu, ray_r_mu_intersects_ground){
  let uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
      planetProps, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  let pixel = textureGetter(uvwz);
  nanCheck(pixel)
  return pixel;
}
export function RayleighPhaseFunction(nu) {
  let k = 3.0 / (16.0 * PI );
  return k * (1.0 + nu * nu);
}

export function MiePhaseFunction(g, nu) {
  let k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
  return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}


export function GetScattering(planetProps, textures, r,  mu,  mu_s,  nu, ray_r_mu_intersects_ground, scatteringOrder) {
  if(scatteringOrder == 1){
    let {singleRayGetter, singleMieGetter} = textures;
    let rayleigh = GetScattering_( planetProps, singleRayGetter, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    let mie = GetScattering_( planetProps, singleMieGetter, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    let _ray = mulS(rayleigh, RayleighPhaseFunction(nu));
    let _mie = mulS(mie, MiePhaseFunction(planetProps.miePhaseFunctionG, nu));
    for(let i = 0; i< _ray.length; ++i){
      _ray[i] += _ray[i];
    }
    return _ray;
  }else{
    let {scattering} = textures;
    return GetScattering_(planetProps, scattering, r,  mu,  mu_s,  nu, ray_r_mu_intersects_ground)
  }
}

function vec3(a,b,c){return new Float32Array([a,b,c])}

export function ComputeIndirectIrradiance(planetProps, single_rayleigh_scattering_texture, single_mie_scattering_texture, multiple_scattering_texture, r, mu_s, scattering_order) {

  const pi = Math.PI;
  const SAMPLE_COUNT = 32;
  const dphi = pi / SAMPLE_COUNT;
  const dtheta = pi / SAMPLE_COUNT;

  let result = [0,0,0,0];

  let omega_s = vec3(sqrt(1.0 - mu_s * mu_s), 0.0, mu_s);
  for (let j = 0; j < SAMPLE_COUNT / 2; ++j) {
    let theta = (Number(j) + 0.5) * dtheta;
    let ray_r_theta_intersects_ground = RayIntersectsGround(planetProps, clampRadius(planetProps, r), cos(theta));
    for (let i = 0; i < 2 * SAMPLE_COUNT; ++i) {
      let phi = (i + 0.5) * dphi;
      let omega = vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
      let domega = (dtheta ) * (dphi ) * sin(theta);

      let nu = dot(omega, omega_s);
      let scattering = GetScattering(
        planetProps,
        {
          singleRayGetter: single_rayleigh_scattering_texture,
          singleMieGetter: single_mie_scattering_texture,
          scattering: multiple_scattering_texture
        },
        r, omega[2],
        mu_s,
        nu,
        ray_r_theta_intersects_ground,
        scattering_order)
      for(let cc =0; cc<3; ++cc)
        result[cc] +=  scattering[cc] * omega[2] * domega;
    }
  }
  return result;
}

function GetScatteringTextureUvwzFromRMuMuSNu( planetProps, r, mu, mu_s, nu, ray_r_mu_intersects_ground) {

  let {topRadius, bottomRadius} = planetProps;
  let {resR, resMus, resNu, resMu, muSMin} = planetProps;

  // Distance to top atmosphere boundary for a horizontal ray at ground level.
  let H = sqrt(topRadius * topRadius - bottomRadius * bottomRadius);
  // Distance to the horizon.
  let rho = safeSqrt(r * r - bottomRadius * bottomRadius);
  let u_r = GetTextureCoordFromUnitRange(rho / H, resR);

  // Discriminant of the quadratic equation for the intersections of the ray
  // (r,mu) with the ground (see RayIntersectsGround).
  let r_mu = r * mu;
  let discriminant = r_mu * r_mu - r * r + bottomRadius * bottomRadius;
  let u_mu;
  if (ray_r_mu_intersects_ground) {
    // Distance to the ground for the ray (r,mu), and its minimum and maximum
    // values over all mu - obtained for (r,-1) and (r,mu_horizon).
    let d = -r_mu - safeSqrt(discriminant);
    let d_min = r - bottomRadius;
    let d_max = rho;
    u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 :
        (d - d_min) / (d_max - d_min), resMu / 2);
    if(isNaN(u_mu)){
      console.log(d, d_min, d_max, (d - d_min) / (d_max - d_min),  GetTextureCoordFromUnitRange( (d - d_min) / (d_max - d_min), resMu / 2));
      throw new Error("==============");
    }
  } else {
    // Distance to the top atmosphere boundary for the ray (r,mu), and its
    // minimum and maximum values over all mu - obtained for (r,1) and
    // (r,mu_horizon).
    let d = -r_mu + safeSqrt(discriminant + H * H);
    let d_min = topRadius - r;
    let d_max = rho + H;
    u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
        (d - d_min) / (d_max - d_min), resMu / 2);
    if(isNaN(u_mu)){
      console.log(r, mu, r_mu, discriminant, H, bottomRadius, rho);
      console.log(d, d_min, d_max, (d - d_min) / (d_max - d_min),  GetTextureCoordFromUnitRange( (d - d_min) / (d_max - d_min), resMu / 2));
      throw new Error("--------");
    }
  }

  let d = DistanceToTopAtmosphereBoundary( planetProps, bottomRadius, mu_s);
  let d_min = topRadius - bottomRadius;
  let d_max = H;
  let a = (d - d_min) / (d_max - d_min);
  let A =
      -2.0 * muSMin * bottomRadius / (d_max - d_min);
  let u_mu_s = GetTextureCoordFromUnitRange(
      max(1.0 - a / A, 0.0) / (1.0 + a), resMus);

  let u_nu = (nu + 1.0) / 2.0;
  nanCheck(vec4(u_nu, u_mu_s, u_mu, u_r),()=>{
    console.log(u_mu, u_r, resMu);
    throw new Error("not good")
  });

  return vec4(u_nu, u_mu_s, u_mu, u_r);
}
function vec4(a,b,c,d){ return new Float32Array([a,b,c,d]);}

export function GetTransmittance( planetProps, transmittanceGetter, r, mu, d, ray_r_mu_intersects_ground) {
  let {topRadius, bottomRadius, muSMin} = planetProps;

  rangeCheck(r, bottomRadius,topRadius);
  rangeCheck(mu, -1.0, 1.0);
  if(d < 0) throw new Error('Negative d:' + d);

  let r_d = clampRadius(planetProps, sqrt(d * d + 2.0 * r * mu * d + r * r));
  let mu_d = clampCosine((r * mu + d) / r_d);

  if (ray_r_mu_intersects_ground) {
    let t1 = GetTransmittanceToTopAtmosphereBoundary(
            planetProps, transmittanceGetter, r_d, -mu_d);
    let t2 = GetTransmittanceToTopAtmosphereBoundary(
            planetProps, transmittanceGetter, r, -mu);
    let t = vdiv(t1, t2)
    let tr =  vmin( t , [1,1,1]);
    nanCheck(tr, ()=>{
      console.log(tr, t1, t2, t);
      throw new Error('err intersect');
    });
    return tr;
  } else {
    let t1 = GetTransmittanceToTopAtmosphereBoundary(
            planetProps, transmittanceGetter, r, mu);
    let t2 = GetTransmittanceToTopAtmosphereBoundary(
            planetProps, transmittanceGetter, r_d, mu_d);
    let t = vdiv(t1, t2);
    let tr = vmin(t, [1,1,1]);
    nanCheck(tr, ()=>{
      console.log(tr, t1, t2, t);
      throw new Error('err not intersect');
    });
    return tr;
  }
}
export function GetRMuMuSNuFromScatteringTextureUvwz(planetProps, uvwz) {
  // x,0 - mus
  // y,1 - nu
  // z,2 - mu
  // w,3 - r
  let {topRadius, bottomRadius, muSMin} = planetProps;
  let {resMu, resNu, resMus, resR} = planetProps;
  //assert(uvwz[0] >= 0.0 && uvwz[0] <= 1.0);
  //assert(uvwz[1] >= 0.0 && uvwz[1] <= 1.0);
  //assert(uvwz[2] >= 0.0 && uvwz[2] <= 1.0);
  //assert(uvwz[3] >= 0.0 && uvwz[3] <= 1.0);
  if(!muSMin) throw new Error("mu s min is not correct");
  let ray_r_mu_intersects_ground;

  // Distance to top atmosphere boundary for a horizontal ray at ground level.
  let H = sqrt(topRadius * topRadius - bottomRadius * bottomRadius);
  // Distance to the horizon.
  let rho = H * GetUnitRangeFromTextureCoord(uvwz[3], resR);
  let r = sqrt(rho * rho + bottomRadius * bottomRadius);
  let mu;

  if (uvwz[2] < 0.5) {
    // Distance to the ground for the ray (r,mu), and its minimum and maximum
    // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
    // we can recover mu:
    let d_min = r - bottomRadius;
    let d_max = rho;
    let d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
        1.0 - 2.0 * uvwz[2], resMu / 2);
    mu = d == 0.0 ? -1.0 :
        clampCosine(-(rho * rho + d * d) / (2.0 * r * d));
    ray_r_mu_intersects_ground = true;
  } else {
    // Distance to the top atmosphere boundary for the ray (r,mu), and its
    // minimum and maximum values over all mu - obtained for (r,1) and
    // (r,mu_horizon) - from which we can recover mu:
    let d_min = topRadius - r;
    let d_max = rho + H;
    let d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
        2.0 * uvwz[2] - 1.0, resMu / 2);
    mu = d == 0.0  ? 1.0 :
        clampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
    ray_r_mu_intersects_ground = false;
  }

  let x_mu_s = GetUnitRangeFromTextureCoord(uvwz[0], resMus);
  let d_min = topRadius - bottomRadius;
  let d_max = H;
  let A = -2.0 * muSMin * bottomRadius / (d_max - d_min);
  let a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
  let d = d_min + min(a, A) * (d_max - d_min);
  let mu_s = d == 0.0  ? 1.0 :
     clampCosine((H * H - d * d) / (2.0 * bottomRadius * d));

  let nu = clampCosine(uvwz[1] * 2.0 - 1.0);
  nanCheck([mu, mu_s, nu, r]);
  return {mu, mu_s, nu, r, ray_r_mu_intersects_ground}
}

export function GetRMuMuSNuFromScatteringTextureFragCoord(planetProps, counters, sizes) {
  let uvwz = new Float32Array(4);
  for(let i =0; i<counters.length;  ++i){
    uvwz[i] = counters[i] / sizes[i];
  }

  let C = GetRMuMuSNuFromScatteringTextureUvwz(planetProps, uvwz);
  // Clamp nu to its valid range of values, given mu and mu_s.
  let {nu, mu, mu_s} = C;
  C.nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
      mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
  return C;
}
export function RayIntersectsGround(planetProps, r, mu) {
  let {topRadius, bottomRadius, muSMin} = planetProps;
  if(r < bottomRadius) throw new Error("r is negative", r);
  rangeCheck(mu, -1, 1);
  return mu < 0.0 && r * r * (mu * mu - 1.0) +
      bottomRadius * bottomRadius >= 0.0 ;
}

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

export function GetIrradiance( planetProps, irradianceGetter, r, mu_s) {
  let uv = GetIrradianceTextureUvFromRMuS(planetProps, r, mu_s);

  nanCheck(uv);
  return irradianceGetter(uv);
}

export function GetIrradianceTextureUvFromRMuS(planetProps, r, mu_s) {
  let {bottomRadius, topRadius} = planetProps;
  let x_r = (r - bottomRadius) / (topRadius - bottomRadius);
  let x_mu_s = mu_s * 0.5 + 0.5;
  let [Wi, Hi] = getIrradianceResolution(planetProps);
  let x = GetTextureCoordFromUnitRange(x_mu_s, Wi);
  let y = GetTextureCoordFromUnitRange(x_r, Hi);
  nanCheck([x,y],()=>{
    console.log(mu_s);
    throw new Error("aaaa")
  })
  return vec2(clamp(x, 0, 1) , clamp(y, 0, 1));
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
  if(v.length == 3)
    return [v[0]*s, v[1]*s, v[2] * s];
  if(v.length == 4)
    return [v[0]*s, v[1]*s, v[2] * s, v[3] *s];
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
  for(let c = 0; c< to.length; ++c){
    to[c] = v1[c] * v2[c];
  }
  //to[1] = v1[1] * v2[1];
  //to[2] = v1[2] * v2[2];
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

export function texture2DGetter(texture, dimensions, components, interpolated = true){
  let b = new Float32Array(2);
  let dims = dimensions.map(x=>x-1)
  return uv=>{
    if(interpolated){
      let coords = getInterpolatedCoords(uv, dimensions);
      let pixel = new Float32Array(components);
      for(let i =0;i<coords.length; ++i){
        let ix = index(coords[i].coords);
        let p = texture.subarray(ix, ix+components);
        for(let c =0; c < components; ++c)
          pixel[c] += coords[i].lerper * p[c];
      }
      return pixel;
    }else{
      b[0] = Math.floor((dimensions[0]-1) * uv[0]);
      b[1] = Math.floor((dimensions[1]-1) * uv[1]);
      let ix = components*(b[1] * dimensions[0] + b[0]);
      if(!texture.subarray)
        console.log(texture);
      return texture.subarray(ix, ix+components);

    }
  }

  function index(b){
    let ix = components*(b[1] * dimensions[0] + b[0]);
    if(ix < 0 || ix > texture.length || ix+components > texture.length)
      console.log(ix, b, dimensions, components, texture.length);
    return ix;

  }
}
export function texture4DGetter(texture, dimensions, components, interpolated = true){
  let b = new Float32Array(4);
  let index = ixGetter(components);
  return coords=>{
    if(interpolated){
      let intr = getInterpolatedCoords(coords, dimensions);
      let pixel = new Float32Array(components);
      for(let i =0; i < intr.length; ++i){
        let cc =  intr[i].coords;
        let ix = index(...intr[i].coords);
        let p = texture.subarray(ix, ix+components);
        for(let c =0; c< components; ++c)
          pixel[c] += intr[i].lerper * p[c];
      }
      return pixel;
    }else{
      vmul4(b, coords, dimensions);
      let i = Math.floor(b[0]);
      let j = Math.floor(b[1]);
      let k = Math.floor(b[2]);
      let r = Math.floor(b[3]);
      ix = index(i,j,k,r);
      let color = texture.subarray(ix, ix+components);
      if(color.length < components){
        throw new Error("color too short "+  ix + ", " + texture.length + ", "+ components);
      }
      return color;
    }
  }

  function ixGetter(components){
    return  (i,j,k,r)=>{
      let X = i*dimensions[1] + j;
      let XX = X * dimensions[2]  + k
      let ix =  components * (XX * dimensions[3] + r);
      return ix;
    }
  }


}



function getInterpolatedCoords(uv, dimensions, SZ){
  if(!SZ) SZ = uv.length;
  let coords = new Float32Array(SZ);
  let dims = dimensions.map(x=>x-1);
  vmul(coords, uv, dims);
  let tCoords  = coords.map(x=>Math.floor(x));
  let lerpCoords = new Float32Array(SZ);
  for(let i = 0; i < SZ; ++i) lerpCoords[i] = coords[i] - tCoords[i];
  let textureSelectors = [];
  let LI = LerpIncrements[SZ];
  for(let i = 0; i< LI.length; ++i){
    let koefs = LI[i];
    let C = new Float32Array(SZ);
    let lerper = 1;
    for(let j = 0; j < koefs.length; ++j){
      let increaser = koefs[j];
      if(tCoords[j] == dims[j]){ // final pixel assumed backwards
        C[j] = tCoords[j] - increaser;
        if(!increaser) lerper *= lerpCoords[j];
        else lerper *= 1 - lerpCoords[j];
      }else{
        C[j] = tCoords[j] + increaser;
        if(!increaser) lerper *= 1 - lerpCoords[j];
        else lerper *= lerpCoords[j];
      }
    }
    textureSelectors.push({coords: C, lerper});
  }
  return textureSelectors;
}

function getInterpolated2DCoords(uv, dimensions){
  let coords = new Float32Array(2);
  let dims = dimensions.map(x=>x-1);
  vmul2(coords, uv, dims);
  let tCoords  = coords.map(x=>Math.floor(x));
  let lerpCoords = new Float32Array(2);
  for(let i = 0; i < 2; ++i) lerpCoords[i] = coords[i] - tCoords[i];
  let textureSelectors = [];
  for(let i = 0; i< LerpIncrements.length; ++i){
    let koefs = LerpIncrements2[i];
    let C = new Float32Array[2];
    let lerper = 1;
    for(let j = 0; j < koefs.length; ++j){
      let increaser = koefs[j];
      if(tCoords[j] == dims[j]){ // final pixel assumed backwards
        C[j] = tCoords[j] - increaser;
        if(!increaser) lerper *= lerpCoords[j];
        else lerper *= 1 - lerpCoords[j];
      }else{
        C[j] = tCoords[j] + increaser;
        if(!increaser) lerper *= 1 - lerpCoords[j];
        else lerper *= lerpCoords[j];
      }
    }
    textureSelectors.push({coords: C, lerper});
  }
  return textureSelectors;


}

function getInterpolated4DCoords(uvwt, dimensions){
  let coords = new Float32Array(4);
  let dims = dimensions.map(x=>x-1);
  vmul4(coords, uvwt, dims);
  let tCoords  = coords.map(x=>Math.floor(x));
  let lerpCoords = new Float32Array(4);
  for(let i = 0; i < 4; ++i) lerpCoords[i] = coords[i] - tCoords[i];
  let textureSelectors = [];
  for(let i = 0; i< LerpIncrements.length; ++i){
    let koefs = LerpIncrements[i];
    let C = new Float32Array[4];
    let lerper = 1;
    for(let j = 0; j < koefs.length; ++j){
      let increaser = koefs[j];
      if(tCoords[j] == dims[j]){ // final pixel assumed backwards
        C[j] = tCoords[j] - increaser;
        if(!increaser) lerper *= lerpCoords[j];
        else lerper *= 1 - lerpCoords[j];
      }else{
        C[j] = tCoords[j] + increaser;
        if(!increaser) lerper *= 1 - lerpCoords[j];
        else lerper *= lerpCoords[j];
      }
    }
    textureSelectors.push({coords: C, lerper});
  }
  return textureSelectors;

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
