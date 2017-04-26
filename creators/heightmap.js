import fs from 'fs';
import zlib from 'zlib';
import path from 'path';
import config from '../config.json';
import {ClassicalNoise} from '../PerlinNoise.js'; 
import {getTextureFilename} from '../utils.js';
import {writeFileInDir} from '../fsUtils.js';
import {stToNormal, calculateTileProperties} from './lookup.js';
import * as zipper from './zipper';

let TextureSize = 512;

function reducer(x){
  return 1/Math.pow(2,x);
}

function testHeightValue(i,j, face){
  let div = Math.floor(TextureSize / (face+1));
  let remainder = i % div;
  let value = 1-Math.abs((remainder / div)*2-1);
  if(i > div) return 0.0;
  return value;
}

const SAMPLES = 8;

function createZeroLod(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = params;
  console.log('create heightMap' ,`lod:${lod}, face:${face}, tile:${tile}`);
  let {radius} = planet;
  let surfaceArea = 4.0 * Math.PI * Math.pow(radius, 2);
  let division = Math.pow(2, lod);
  let thisTileProps = calculateTileProperties(face, lod, tile);

  let generator = new ClassicalNoise(planet.table);
  let ab = new Float32Array(TextureSize * TextureSize);

  let levels = [0,1,2,3,4,5,6,7,8];

  let _arrayStart = Date.now();
  let _avgNoise = [];
  let SampleFrom = lod * SAMPLES;
  let SampleTo = SampleFrom + SAMPLES;
  for(let i =0; i < TextureSize; ++i){
    for(let j =0; j < TextureSize; ++j){
      let tt = i / TextureSize / division;
      let ts = j / TextureSize / division;
      let normal = stToNormal(thisTileProps.s+ts, thisTileProps.t+tt, face)
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

      //***********
      // heightValue = testHeightValue(i,j, face);

      //**************
      // heightValue = face;
      /////////////////
      //
      // if(j  == 128){
        //heightValue = (face -1)%6
      //}
      //if(i  == 128){
        //heightValue = (face + 1)%6
      //}

      let ix = (j*TextureSize + i); 
      ab[ix] = heightValue;
    }
  }
  let buffer = Buffer.from(ab.buffer);
  let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});
  zipper.deflateTo(filePath, buffer, callback);
}

export function isExists(forPlanet, params, callback){
  let filename = getTextureFilename({
    planetUUID: forPlanet.uuid,
    textureType: 'height',
    lod: params.lod, 
    tile: params.tile, 
    face: params.face
  });
  fs.access(filename, err=>{
    if(!err) return callback(true);
    callback(false);
  })
}

export function getRequirements(){
  return [];
}

export function create(input, callback){
  let {lod, face, tile} = input.params;


  if(lod == 0) return createZeroLod(input, callback);



  let prevLod = lod - 1;

  let thisTileProps = calculateTileProperties(face, lod, tile);


  let division = Math.pow(2, lod);
  let prevDivision = Math.pow(2, lod-1);
  let S = thisTileProps.s * prevDivision;
  let T = thisTileProps.t * prevDivision;
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
    let division = Math.pow(2, lod);

    let thisTileProp = calculateTileProperties(face, lod, tile);

    let generator = new ClassicalNoise(planet.table);
    let ab = new Float32Array(TextureSize*TextureSize);

    let addColor = [0,0,0,0]
    let SampleFrom = 0;
    let SampleTo = SAMPLES+lod; 

    for(let i =0; i < TextureSize; ++i){
      for(let j =0; j < TextureSize; ++j){
        let ts = j / TextureSize / division;
        let tt = i / TextureSize / division;
        let normal = stToNormal(thisTileProp.s+ts, thisTileProp.t+tt, face)
        let [x,y,z] = normal;
        let normalLength = Math.sqrt(x*x + y*y + z*z);
        let heightValue = 0;

        // calculate prev texture coords;
        let I = Math.floor((thisTileProp.t+tt-prevTile.t) * prevTile.division * TextureSize);
        let J = Math.floor((thisTileProp.s+ts-prevTile.s) * prevTile.division * TextureSize);
        let prevIx = (J * TextureSize + I ); // *4;

        for(let cc = SampleFrom; cc < SampleTo; ++cc){
          let l = cc;
          let ll = Math.pow(2,l)/normalLength;
          let noiseLevel =  generator.noise(normal[0]*ll, normal[1]*ll, normal[2]*ll);
          noiseLevel *= 1/Math.pow(2,l);
          heightValue += noiseLevel;
        }
        
        let ix = (j*TextureSize + i); 

        ab[ix] = heightValue; 
      }
    }

    let filePath = getTextureFilename({planetUUID:planet.uuid, textureType:'height', lod, tile, face});

    let buffer = Buffer.from(ab.buffer);

    zipper.deflateTo(filePath, buffer, (...args)=>{
      callback(...args)
    });
  }
}



