import {
  ComputeMultipleScattering, 
  index4,
  nanCheck
} from './utils.js'
const components = 4;
export function ComputeMultipleScatteringTexture(planetProps, delta_multiple_scattering, scattering, transmittanceGetter, scatteringDensityGetter, scattering_order){

  let {resMu, resNu, resR, resMus} = planetProps;

  let counters = new Int32Array(4);
  let sizes = new Int32Array([resMus, resNu, resMu, resR]);
  for(;counters[0] < sizes[0]; ++counters[0]){
    console.log(counters);
    for(counters[1] = 0;counters[1] < sizes[1]; ++counters[1]){
      for(counters[2] = 0;counters[2] < sizes[2]; ++counters[2]){
        for(counters[3] = 0;counters[3] < sizes[3]; ++counters[3]){
          let ix = components * index4(counters, sizes);
          let scatteringComponent = ComputeMultipleScattering(
            planetProps,
            transmittanceGetter,
            scatteringDensityGetter, counters, sizes);
          nanCheck(scatteringComponent);
          for(let c=0; c<components;++c){
            delta_multiple_scattering[ix+c] =  scatteringComponent[c];
            scattering[ix+c] = scatteringComponent[c];
          }
        }
      }
    }
  }
}
