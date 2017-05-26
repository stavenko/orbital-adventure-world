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

  getProfileDensity,
  vmul1, vdiv,
  assert
} from './utils.js'

const components = 4;
const sqrt = Math.sqrt;
const min = Math.min;
export function computeSingleScattering(planetProps, transmittanceTexture){
  let transmittanceRes = getTransmittanceResolution(planetProps);
  let transmittanceGetter = texture2DGetter(transmittanceTexture, transmittanceRes, components);
  let {resMu, resNu, resR, resMus} = planetProps

  let arraySize = resNu * resMu * resMus * resR * components
  let deltaRayleight = new Float32Array(arraySize);
  let deltaMie = new Float32Array(arraySize);
  let scattering = new Float32Array(arraySize);
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
            deltaRayleight[ix + i] = delta_rayleigh[i];
            deltaMie[ix + i] = delta_mie[i];
            scattering[ix + i] = scattering[i];
          }
        }
      }
    }
  }
  return {deltaRayleight, deltaMie, scatteringTexture:scattering};
}

function index4(counters, sizes){
  return ((counters[0] * sizes[1] + counters[1])*sizes[2] + counters[2])*sizes[3] + counters[3];
}

function ComputeSingleScatteringTexture(planetProps, transmittanceGetter, 
                                        counters, sizes){

  let coords = GetRMuMuSNuFromScatteringTextureFragCoord(planetProps, counters, sizes);

  return ComputeSingleScattering(planetProps, transmittanceGetter, coords);
}

function GetRMuMuSNuFromScatteringTextureFragCoord(planetProps, counters, sizes) {
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

function GetRMuMuSNuFromScatteringTextureUvwz(planetProps, uvwz) {
  // x,0 - mus
  // y,1 - nu
  // z,2 - mu
  // w,3 - r
  let {topRadius, bottomRadius, muSMin} = planetProps;
  let {resMu, resNu, resMus, resR} = planetProps;
  assert(uvwz[0] >= 0.0 && uvwz[0] <= 1.0);
  assert(uvwz[1] >= 0.0 && uvwz[1] <= 1.0);
  assert(uvwz[2] >= 0.0 && uvwz[2] <= 1.0);
  assert(uvwz[3] >= 0.0 && uvwz[3] <= 1.0);
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
  nanCheck(rayleigh,()=>{
    console.log(rayleigh, rayleigh_sum, dx, solarIrradiance, rayleighScattering)
    throw new Error("Err");
  });
  nanCheck(mie);
  return {delta_rayleigh: rayleigh, delta_mie: mie};
}

function DistanceToNearestAtmosphereBoundary(planetProps, r, mu, ray_r_mu_intersects_ground) {
  if (ray_r_mu_intersects_ground) {
    return DistanceToBottomAtmosphereBoundary(planetProps, r, mu);
  } else {
    return DistanceToTopAtmosphereBoundary(planetProps, r, mu);
  }
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

function RayIntersectsGround(planetProps, r, mu) {
  let {topRadius, bottomRadius, muSMin} = planetProps;
  if(r < bottomRadius) throw new Error("r is negative");
  rangeCheck(mu, -1, 1);
  return mu < 0.0 && r * r * (mu * mu - 1.0) +
      bottomRadius * bottomRadius >= 0.0 ;
}

function GetTransmittance( planetProps, transmittanceGetter, r, mu, d, ray_r_mu_intersects_ground) {
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

