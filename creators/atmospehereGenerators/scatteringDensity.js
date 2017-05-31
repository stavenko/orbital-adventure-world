import {
  GetRMuMuSNuFromScatteringTextureFragCoord,
  DistanceToBottomAtmosphereBoundary,
  GetTransmittance,
  dot,
  RayIntersectsGround,
  GetScattering, 
  index4,
  nanCheck,
  mulS,
  vadd,
  normalize,
  GetIrradiance,
  getProfileDensity,
  RayleighPhaseFunction,
  MiePhaseFunction

} from './utils';

const sqrt = Math.sqrt;
const atan = Math.atan;
const max = Math.max;
const min = Math.min;
const cos = Math.cos;
const sin = Math.sin;
const exp = Math.exp;
const pow = Math.pow;
const PI = Math.PI;
const components = 4;
function vec3(z,x,c){ return [z,x,c]};

export function computeScatteringDensity(planetProps, scatteringDensity, transmittanceGetter, singleRayleighScatteringGetter, singleMieScatteringGetter, multipleScatteringGetter, irradianceGetter, scattering_order){
  let {resMu, resNu, resR, resMus} = planetProps;

  let counters = new Int32Array(4);
  let sizes = new Int32Array([resMus, resNu, resMu, resR]);
  for(;counters[0] < sizes[0]; ++counters[0]){
    console.log(counters);
    for(counters[1] = 0;counters[1] < sizes[1]; ++counters[1]){
      for(counters[2] = 0;counters[2] < sizes[2]; ++counters[2]){
        for(counters[3] = 0;counters[3] < sizes[3]; ++counters[3]){
          let ix = components * index4(counters, sizes);
          let scattering = ComputeScatteringDensity(
            planetProps, 
            transmittanceGetter, 
            singleRayleighScatteringGetter, 
            singleMieScatteringGetter, 
            multipleScatteringGetter,
            irradianceGetter, 
            scattering_order, counters, sizes);
          nanCheck(scattering);
          for(let i=0; i<components;++i){
            scatteringDensity[ix + i] = scattering[i];
          }
        }
      }
    }
  }
  let S = 0;
  for(let i =0; i< scatteringDensity.length; ++i) S+= scatteringDensity[i];
  if(S < 1) throw new Error("Total sum of scattering density is less than 1", S);
  return null;
}


function ComputeScatteringDensity(planetProps, transmittanceGetter, singleRayleighGetter, singleMieGetter, multipleScatteringGetter, irradianceGetter, scattering_order, counters, sizes){
  let coords = GetRMuMuSNuFromScatteringTextureFragCoord(planetProps, counters, sizes);

  let {mu, nu, mu_s, r, ray_r_mu_intersects_ground} = coords;
  nanCheck([mu, nu, r, mu_s]);

  // Compute unit direction vectors for the zenith, the view direction omega and
  // and the sun direction omega_s, such that the cosine of the view-zenith
  // angle is mu, the cosine of the sun-zenith angle is mu_s, and the cosine of
  // the view-sun angle is nu. The goal is to simplify computations below.
  let zenith_direction = vec3(0.0, 0.0, 1.0);
  let omega = vec3(sqrt(1.0 - mu * mu), 0.0, mu);
  let sun_dir_x = omega[0] == 0.0 ? 0.0 : (nu - mu * mu_s) / omega[0];
  let sun_dir_y = sqrt(max(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
  let omega_s = vec3(sun_dir_x, sun_dir_y, mu_s);

  const SAMPLE_COUNT = 16;
  const pi = Math.PI;
  const dphi = pi / SAMPLE_COUNT;
  const dtheta = pi / SAMPLE_COUNT;
  let rayleigh_mie = [0, 0, 0];

  // Nested loops for the integral over all the incident directions omega_i.
  for (let l = 0; l < SAMPLE_COUNT; ++l) {
    let theta = (l + 0.5) * dtheta;
    let cos_theta = cos(theta);
    let sin_theta = sin(theta);
    let ray_r_theta_intersects_ground = RayIntersectsGround(planetProps, r, cos_theta);

    // The distance and transmittance to the ground only depend on theta, so we
    // can compute them in the outer loop for efficiency.
    let distance_to_ground = 0.0 ;
    let transmittance_to_ground = vec3(0.0, 0, 0);
    let ground_albedo = vec3(0, 0, 0);
    if (ray_r_theta_intersects_ground) {
      distance_to_ground = DistanceToBottomAtmosphereBoundary(planetProps, r, cos_theta);
      transmittance_to_ground =
          GetTransmittance(planetProps, transmittanceGetter, r, cos_theta,
              distance_to_ground, ray_r_theta_intersects_ground);
      ground_albedo = [
        planetProps.groundAlbedo,
        planetProps.groundAlbedo,
        planetProps.groundAlbedo
      ];
    }
    nanCheck(ground_albedo, ()=>{
      console.log(ground_albedo, planetProps);
      throw new Error("aaaa");
    })


    for (let m = 0; m < 2 * SAMPLE_COUNT; ++m) {
      let phi = (m + 0.5) * dphi;
      let omega_i = vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
      let domega_i = dtheta * dphi * sin(theta);

      // The radiance L_i arriving from direction omega_i after n-1 bounces is
      // the sum of a term given by the precomputed scattering texture for the
      // (n-1)-th order:
      let nu1 = dot(omega_s, omega_i);
      let incident_radiance = GetScattering(
        planetProps,
        {
          singleRayGetter:singleRayleighGetter, 
          singleMieGetter,
          scattering: multipleScatteringGetter
        }, 
        r, omega_i[2], mu_s, nu1, 
        ray_r_theta_intersects_ground, 
        scattering_order - 1);

      // and of the contribution from the light paths with n-1 bounces and whose
      // last bounce is on the ground. This contribution is the product of the
      // transmittance to the ground, the ground albedo, the ground BRDF, and
      // the irradiance received on the ground after n-2 bounces.
      let zd = mulS(zenith_direction, r);
      let od = mulS(omega_i, distance_to_ground);
      let zdod = vadd(zd, od);
      let ground_normal = normalize(zdod);


      let ground_irradiance = GetIrradiance(
          planetProps, irradianceGetter, planetProps.bottomRadius,
          dot(ground_normal, omega_s));

      //console.log(transmittance_to_ground, ground_albedo );
      for(let ii = 0; ii< 3; ++ii){
        incident_radiance[ii] += transmittance_to_ground[ii] * ground_albedo[ii] * (1/PI) * ground_irradiance[ii];
      }
      nanCheck(incident_radiance);

      // The radiance finally scattered from direction omega_i towards direction
      // -omega is the product of the incident radiance, the scattering
      // coefficient, and the phase function for directions omega and omega_i
      // (all this summed over all particle types, i.e. Rayleigh and Mie).
      let nu2 = dot(omega, omega_i);
      let rayleigh_density = getProfileDensity(
          planetProps.rayleighDensity, r - planetProps.bottomRadius);
      let mie_density = getProfileDensity(
          planetProps.mieDensity, r - planetProps.bottomRadius);
      for(let ii =0 ;ii < 3; ++ii){
        rayleigh_mie[ii] += incident_radiance[ii] * (
            planetProps.rayleighScattering[ii] * rayleigh_density *
                RayleighPhaseFunction(nu2) +
            planetProps.mieScattering[ii] * mie_density *
                MiePhaseFunction(planetProps.miePhaseFunctionG, nu2)) *
            domega_i;
      }
      nanCheck(rayleigh_mie, ()=>{
        console.log( incident_radiance , planetProps.rayleighScattering , rayleigh_density ,
                RayleighPhaseFunction(nu2) , planetProps.mieScattering , mie_density ,
                MiePhaseFunction(planetProps.miePhaseFunctionG, nu2) ,planetProps, nu2, domega_i);
        throw new Error("-----------");
      });
    }
  }
  return rayleigh_mie;
}

