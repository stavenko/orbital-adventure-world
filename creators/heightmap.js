import fs from 'fs';
import zlib from 'zlib';
import path from 'path';
import config from '../config.json';
import {ClassicalNoise} from '../PerlinNoise.js'; 
import {getTextureFilename} from '../utils.js';
import {writeFileInDir} from '../fsUtils.js';

let TextureSize = 2048;

function reducer(x){
  return 1/Math.pow(2,x);
}

export function create(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = params;
  console.log('create' ,lod, face, tile);
  let {radius} = planet;
  let surfaceArea = 4.0 * Math.PI * Math.pow(radius, 2);
  let division = Math.pow(2, lod);
  let S = Math.floor(tile / division);
  let T = tile % division;

  let s = S / division;
  let t = T / division;

  let generator = new ClassicalNoise(planet.table);
  let ab = new Buffer(TextureSize*TextureSize * 4);

  let levels = [0,1,2,3,4,5,6,7,8];

  let _arrayStart = Date.now();
  let _avgNoise = [];
  for(let i =0; i < TextureSize; ++i){
    for(let j =0; j < TextureSize; ++j){
      let ts = i / TextureSize / division;
      let tt = j / TextureSize / division;
      let normal = stToNormal(s+ts, t+tt, face)
      let [x,y,z] = normal;
      let normalLength = Math.sqrt(x*x + y*y + z*z);
      // normal = normal.map(x=>x/length);
      let heightValue = 0;

      let _noiseStart = Date.now();
      for(let cc = 0; cc < levels.length; ++cc){
        let l = levels[cc];
        let ll = Math.pow(2,l+lod)/normalLength;
        let noiseLevel =  generator.noise(normal[0]*ll, normal[1]*ll, normal[2]*ll);
        noiseLevel *= reducer(l);
        heightValue += noiseLevel;
      }

      _avgNoise.push(Date.now() - _noiseStart);

      let ix = (j*TextureSize + i)*4;
      ab[ix] = (heightValue+1)*0.5*255;// (theta + 1)/2 * 255;
      ab[ix+1] =(heightValue+1)*0.5*255 //; Math.abs(phi) * 255;
      ab[ix+2] = (heightValue+1)*0.5*255;
      ab[ix+3] = 255;
    }
  }
  let l = _avgNoise.length;
  let n = _avgNoise.reduce((a,b)=>a+b);
  console.log('avg noise calc', n/l);
  console.log("array time", Date.now() - _arrayStart);
  let gz = zlib.createGzip();

  let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});

  let _zipStart = Date.now();
  zlib.deflate(ab, {level:9}, (err, buffer)=>{
    if(err) throw err;
    // console.log("writing file to ", filePath);
    writeFileInDir(filePath, buffer, err=>{

      if(err) throw err;
      console.log("zipping done", Date.now() - _zipStart)
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
