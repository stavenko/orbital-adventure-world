import {mulS, 
  limit, 
  getTransmittanceResolution, 
  unitRangeFromTextureCoord,
  clampCosine,
  safeSqrt,
  clampDistance,
  getProfileDensity,
  assert,
  clamp,
  rangeCheck,
  GetRMuFromTransmittanceTextureUv,
  DistanceToTopAtmosphereBoundary
} from './utils.js';

const sqrt = Math.sqrt;
const exp = Math.exp;




function computeOpticalLengthToTopAtmosphereBoundary(planetProps, profile, r, mu) {
  let {topRadius, bottomRadius} = planetProps;
  r = clamp(r, bottomRadius, topRadius);
  rangeCheck(r, bottomRadius, topRadius, 'R check in optical depth');
  rangeCheck(mu,-1.0,1.0, 'Mu check');
  // Number of intervals for the numerical integration.
  const SAMPLE_COUNT = 500;
  // The integration step, i.e. the length of each integration interval.
  let dx = DistanceToTopAtmosphereBoundary(planetProps, r, mu) / SAMPLE_COUNT;
  // Integration loop.
  let result = 0.0;
  for (let i = 0; i <= SAMPLE_COUNT; ++i) {
    let d_i = i * dx;
    // Distance between the current sample point and the planet center.
    let r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
    // Number density at the current sample point (divided by the number density
    // at the bottom of the atmosphere, yielding a dimensionless number).
    let y_i = getProfileDensity(profile, r_i - bottomRadius);
    // Sample weight (from the trapezoidal rule).
    let weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
    result += y_i * weight_i * dx;
  }
  return result;
}

function computeTransmittanceToTopAtmospehereBoundary(planetProps, r, mu){
  let {bottomRadius, topRadius} = planetProps;
  // rangeCheck(r, bottomRadius, topRadius);
  assert(mu >= -1 && mu <= 1, 'Incorrect Mu');
  let {
    rayleighScattering, 
    mieExtinction, 
    absorptionExtinction,
    rayleighDensity,
    mieDensity,
    absorptionDensity
  } = planetProps;

  let rayd = computeOpticalLengthToTopAtmosphereBoundary(planetProps, rayleighDensity, r, mu);
  let mie = computeOpticalLengthToTopAtmosphereBoundary(planetProps, mieDensity, r, mu);
  let absorp = computeOpticalLengthToTopAtmosphereBoundary(planetProps, absorptionDensity, r, mu);
  let ray = mulS(rayleighScattering, rayd);
  mie = mulS(mieExtinction, mie);
  absorp = mulS(absorptionExtinction, absorp);
  return [
    exp(-(ray[0] + mie[0] + absorp[0])),
    exp(-(ray[1] + mie[1] + absorp[1])),
    exp(-(ray[2] + mie[2] + absorp[2])),
    0
  ];
}

export function getTransmittenceColor({s,t}, planetProps){

  let {r,mu} = GetRMuFromTransmittanceTextureUv(planetProps, s,t);
  let result = computeTransmittanceToTopAtmospehereBoundary(planetProps, r, mu);

  return result;
}

