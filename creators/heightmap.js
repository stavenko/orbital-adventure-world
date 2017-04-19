import fs from 'fs';
import zlib from 'zlib';
import path from 'path';
import config from '../config.json';
import {ClassicalNoise} from '../PerlinNoise.js'; 
import {getTextureFilename} from '../utils.js';
import {writeFileInDir} from '../fsUtils.js';
import {stToNormal} from './lookup.js';
import * as zipper from './zipper';

let TextureSize = 512;

function reducer(x){
  return 1/Math.pow(2,x);
}

const SAMPLES = 8;

function createZeroLod(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = params;
  console.log('create heightMap' ,`lod:${lod}, face:${face}, tile:${tile}`);
  let {radius} = planet;
  let surfaceArea = 4.0 * Math.PI * Math.pow(radius, 2);
  let division = Math.pow(2, lod);
  let T = Math.floor(tile / division);
  let S = tile % division;

  let s = S / division;
  let t = T / division;

  let generator = new ClassicalNoise(planet.table);
  let ab = new Float32Array(TextureSize * TextureSize);

  let levels = [0,1,2,3,4,5,6,7,8];

  let _arrayStart = Date.now();
  let _avgNoise = [];
  let SampleFrom = lod * SAMPLES;
  let SampleTo = SampleFrom + SAMPLES;
  for(let i =0; i < TextureSize; ++i){
    for(let j =0; j < TextureSize; ++j){
      let ts = i / TextureSize / division;
      let tt = j / TextureSize / division;
      let normal = stToNormal(s+ts, t+tt, face)
      let [x,y,z] = normal;
      let normalLength = Math.sqrt(x*x + y*y + z*z);
      let heightValue = 0;

      let _noiseStart = Date.now();
      
      for(let cc = SampleFrom; cc < SampleTo; ++cc){
        let l = cc;
        let ll = Math.pow(2,l)/normalLength;
        let noiseLevel =  generator.noise(normal[0]*ll, normal[1]*ll, normal[2]*ll);
        noiseLevel *= 1/Math.pow(2,l);;
        heightValue += noiseLevel;
      }

      _avgNoise.push(Date.now() - _noiseStart);

      let ix = (j*TextureSize + i); 
      ab[ix] = heightValue;
    }
  }
  let l = _avgNoise.length;
  let n = _avgNoise.reduce((a,b)=>a+b);
  console.log('avg noise calc', n/l);
  console.log("array time", Date.now() - _arrayStart);
  let buffer = Buffer.from(ab.buffer);
  console.log("buffer size", buffer.length);
  let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});
  zipper.deflateTo(filePath, buffer, callback);
}

export function create(input, callback){
  let {lod, face, tile} = input.params;
  if(lod == 0) return createZeroLod(input, callback);
  console.log('==create heightMap' ,`lod:${lod}, face:${face}, tile:${tile}`);

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

   
  let tup = createTileUppersLods(input,{
    tile: tileNum,
    s:Math.floor(S) / prevDivision,
    t:Math.floor(T) / prevDivision,
    tileSize: 1/prevDivision,
    division: prevDivision,
  }, callback);

  tup(new Buffer(TextureSize*TextureSize*4));

  //zipper.inflateFrom(filePath, createTileUppersLods(input,{
    //tile: tileNum,
    //s:Math.floor(S) / prevDivision,
    //t:Math.floor(T) / prevDivision,
    //tileSize: 1/prevDivision,
    //division: prevDivision,
  //}, callback))


}

const EPSILON = 1e-6;

function eq(a,b){
  return Math.abs(a-b) < EPSILON;
}

function createTileUppersLods(input, prevTile, callback){
  return prevousTileBuffer=>{
    let {planet, params} = input;
    let {lod, face, tile} = params;
    let prevousTile = new Float32Array(prevousTileBuffer.buffer);
    //if(lod == 2) {
      //console.log("skip 2");
      //return callback(null);
    //}
    console.log('create' ,`lod:${lod}, face:${face}, tile:${tile}`);
    let division = Math.pow(2, lod);
    let S = Math.floor(tile / division);
    let T = tile % division;


    let s = S / division; // this texture start
    let t = T / division; 
    console.log('create' ,`T:${T}, S:${S}, t:${t}, s:${s}`);


    let generator = new ClassicalNoise(planet.table);
    let ab = new Float32Array(TextureSize*TextureSize);

    let addColor = [0,0,0,0]
    let SampleFrom = 0;//  lod * SAMPLES;
    let SampleTo = SAMPLES+lod; //SampleFrom + SAMPLES;

    for(let i =0; i < TextureSize; ++i){
      for(let j =0; j < TextureSize; ++j){
        let ts = i / TextureSize / division;
        let tt = j / TextureSize / division;
        let normal = stToNormal(s+ts, t+tt, face)
        let [x,y,z] = normal;
        let normalLength = Math.sqrt(x*x + y*y + z*z);
        let heightValue = 0;

        // calculate prev texture coords;
        let I = Math.floor((t+tt-prevTile.t) * prevTile.division * TextureSize);
        let J = Math.floor((s+ts-prevTile.s) * prevTile.division * TextureSize);
        let prevIx = (J * TextureSize + I ); // *4;

        for(let cc = SampleFrom; cc < SampleTo; ++cc){
          let l = cc;
          let ll = Math.pow(2,l)/normalLength;
          let noiseLevel =  generator.noise(normal[0]*ll, normal[1]*ll, normal[2]*ll);
          noiseLevel *= 1/Math.pow(2,l);
          heightValue += noiseLevel;
        }
        
        let ix = (j*TextureSize + i); // * 4;
        // console.log(heightValue/prevousTile[prevIx]);

        ab[ix] = heightValue; //+addColor[0];
        //ab[ix] = prevousTile[prevIx] + heightValue; //+addColor[0];
        //ab[ix] = prevousTile[prevIx]+addColor[0];
        //ab[ix+1] = prevousTile[prevIx+1]+addColor[1];
        //ab[ix+2] = prevousTile[prevIx+2]+addColor[2];
        //ab[ix+3] = 255;
      }
    }

    let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});

    let buffer = Buffer.from(ab.buffer);

    zipper.deflateTo(filePath, buffer, (...args)=>{
      console.log("done");
      callback(...args)
    });
  }
}



