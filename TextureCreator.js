import {ClassicalNoise} from './PerlinNoise.js';

const TextureSize = 2048


function createHeightMapTexture(planet, face, lod, tile){
  let division = Math.pow(2, lod);
  let S = Math.floor(tile / division);
  let T = tile % division;

  let s = S / division;
  let t = T / division;

  let generator = new ClassicalNoise;

  for(let i =0; i < TextureSize; ++i){
    for(let j =0; j < TextureSize; ++j){
      let ts = i / TextureSize / division;
      let tt = j / TextureSize / division;
      let normal = stToNormal(s+ts, t+tt, face)
      let noiseLevel =  generator.noise(...normal);
      console.log(noiseLevel);
    }
  }
}

export function generateTexture(planet, params, callback){
  // assume we will want to make this in worker
  
}

function normalize(n){
  let l = Math.sqrt(n[0] * n[0] +  n[1] * n[1] + n[2] * n[2]);
  return n.map(i=>i/l);
}
function stToNormal(s,t, face){
  let ss = s * 2 - 1;
  let tt = t * 2 - 1;

  return [
    [-1, -tt, ss], //back
    [ss, -1, -tt], // left
    [1, -tt,-ss], // front
    [ss,  1, tt], // right
    [ss, -tt, 1], // top
    [-ss, -tt, -1], // bottom
  ][face];
}
