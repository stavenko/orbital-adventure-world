import fs from 'fs';
import zlib from 'zlib';
import path from 'path';
import config from '../config.json';
import {ClassicalNoise} from '../PerlinNoise.js'; 
import {getTextureFilename} from '../utils.js';
import {writeFileInDir} from '../fsUtils.js';

let TextureSize = 2048;

export function create(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = params;
  let {radius} = planet;
  let surfaceArea = 4.0 * Math.PI * Math.pow(radius, 2);
  let division = Math.pow(2, lod);
  let S = Math.floor(tile / division);
  let T = tile % division;

  let s = S / division;
  let t = T / division;

  let generator = new ClassicalNoise(planet.table);

  let ab = new Buffer(TextureSize*TextureSize * 4);

  for(let i =0; i < TextureSize; ++i){
    for(let j =0; j < TextureSize; ++j){
      let ts = i / TextureSize / division;
      let tt = j / TextureSize / division;
      let normal = stToNormal(s+ts, t+tt, face)
      let [x,y,z] = normal;
      let length = Math.sqrt(x*x + y*y + z*z);
      normal = normal.map(x=>x/length*2);

      let noiseLevel =  (generator.noise(...normal)+1)/2;

      let theta = Math.atan2(y, x) / (Math.PI);
      let r = Math.hypot(x, y) ;
      let phi = Math.atan2(z, r)/(Math.PI/2);

      let ix = (j*TextureSize + i)*4;
      ab[ix] = (noiseLevel+1)*0.5*255;// (theta + 1)/2 * 255;
      ab[ix+1] =(noiseLevel+1)*0.5*255 //; Math.abs(phi) * 255;
      ab[ix+2] = (noiseLevel+1)*0.5*255;
      ab[ix+3] = 255;
    }
  }
  let gz = zlib.createGzip();

  let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});

  zlib.deflate(ab, {level:9}, (err, buffer)=>{
    if(err) throw err;
    console.log("writing file to ", filePath);
    writeFileInDir(filePath, buffer, err=>{
      if(err) throw err;
      callback(null);
      
    });
  })

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
