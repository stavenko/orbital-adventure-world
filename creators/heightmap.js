import fs from 'fs';
import zlib from 'zlib';
import path from 'path';
import config from '../config.json';
import {ClassicalNoise} from '../PerlinNoise.js'; 
import {getTextureFilename} from '../utils.js';
import {writeFileInDir} from '../fsUtils.js';
import * as zipper from './zipper';

let TextureSize = 2048;

function reducer(x){
  return 1/Math.pow(2,x);
}

const SAMPLES = 8;

function createZeroLod(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = params;
  console.log('create' ,`lod:${lod}, face:${face}, tile:${tile}`);
  let {radius} = planet;
  let surfaceArea = 4.0 * Math.PI * Math.pow(radius, 2);
  let division = Math.pow(2, lod);
  let T = Math.floor(tile / division);
  let S = tile % division;

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
      for(let cc = 0; cc < SAMPLES; ++cc){
        let l = cc+lod;
        let ll = Math.pow(2,l+lod)/normalLength;
        let noiseLevel =  generator.noise(normal[0]*ll, normal[1]*ll, normal[2]*ll);
        noiseLevel *= 1/Math.pow(2, l);
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
  let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});
  zipper.deflateTo(filePath, ab, callback);
}

export function create(input, callback){
  let {lod, face, tile} = input.params;
  if(lod == 0) return createZeroLod(input, callback);
  console.log('==create' ,`lod:${lod}, face:${face}, tile:${tile}`);

  let prevLod = lod - 1;
  let division = Math.pow(2, lod);
  let prevDivision = Math.pow(2, lod-1);
  let T = Math.floor(tile / division);
  let S = tile % division;
  let s = S / division;
  let t = T / division;
  console.log('lod', lod, prevLod, T, S, tile);
  S = s * prevDivision;
  T = t * prevDivision;
  let tileS = S - Math.floor(S);
  let tileT = T - Math.floor(T);
  let tileNum = Math.floor(T) * prevDivision + Math.floor(S);

  let filePath = getTextureFilename({planetUUID:input.planet.uuid, textureType:'height', lod:lod-1, tile:tileNum, face});

   
  zipper.inflateFrom(filePath, createTileUppersLods(input,{
    tile: tileNum,
    s:Math.floor(S) / prevDivision,
    t:Math.floor(T) / prevDivision,
    tileSize: 1/prevDivision,
    division: prevDivision,
  }, callback))


}

const EPSILON = 1e-6;

function eq(a,b){
  return Math.abs(a-b) < EPSILON;
}

function createTileUppersLods(input, prevTile, callback){
  return prevousTile=>{
    let {planet, params} = input;
    let {lod, face, tile} = params;
    //if(lod == 2) {
      //console.log("skip 2");
      //return callback(null);
    //}
    console.log('create' ,`lod:${lod}, face:${face}, tile:${tile}`);
    let division = Math.pow(2, lod);
    let T = Math.floor(tile / division);
    let S = tile % division;


    let s = S / division; // this texture start
    let t = T / division; 
    console.log('create' ,`T:${T}, S:${S}, t:${t}, s:${s}`);


    let generator = new ClassicalNoise(planet.table);
    let ab = new Buffer(TextureSize*TextureSize * 4);

    let addColor = [0,0,0,0]

    for(let i =0; i < TextureSize; ++i){
      for(let j =0; j < TextureSize; ++j){
        let ts = j / TextureSize / division;
        let tt = i / TextureSize / division;
        let normal = stToNormal(s+ts, t+tt, face)
        let [x,y,z] = normal;
        let normalLength = Math.sqrt(x*x + y*y + z*z);
        let heightValue = 0;

        // calculate prev texture coords;
        let I = Math.floor((t+tt-prevTile.t) * prevTile.division * TextureSize);
        let J = Math.floor((s+ts-prevTile.s) * prevTile.division * TextureSize);
        let prevIx = (J * TextureSize + I )*4;
        
        
        let ix = (j*TextureSize + i) * 4;

        ab[ix] = prevousTile[prevIx]+addColor[0];
        ab[ix+1] = prevousTile[prevIx+1]+addColor[1];
        ab[ix+2] = prevousTile[prevIx+2]+addColor[2];
        ab[ix+3] = 255;
      }
    }

    let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});

    zipper.deflateTo(filePath, ab, (...args)=>{
      console.log("done");
      callback(...args)
    });
  }
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
