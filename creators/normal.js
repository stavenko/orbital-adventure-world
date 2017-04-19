import async from 'async';
import {getTileProps, transformToCorrectFace, stToNormal} from './lookup.js'
import {getTextureFilename} from '../utils.js';

let TextureSize = 512;
export function calculateNormalMap(planet, lod, face,tile){
  let division = Math.pow(2, lod);
  let T = Math.floor(tile / division);
  let S = tile % division;

  let s = S / division;
  let t = T / division;

  let filesOpened = {};

  let thisTileProps = calculateTileProperties(face, lod, tile);

  let normalMap = new Buffer(TextureSize * TextureSize * 3);
  for(let i = 0; i < TextureSize; ++i){
    for(let j = 0; j < TextureSize; ++j){
      let vs = [[i-1, j],
       [i+1, j],

       [i, j-1],
       [i, j+1]]

       async.map(vs, lookupAt(tileProps), results=>{

       })
    }
  }

  function lookupAt(tileProps){
    return ([i,j], next)=>{
      let ts = i / TextureSize / division;
      let tt = j / TextureSize / division;
      let {t,s,face} = tileProps;
      let coords = {...transformToCorrectFace(s+ts, t+tt, face)};

      let normal = stToNormal(coords.s, coords.t, coords.face);

      getHeightAt(getTileProps(coords, lod), coords, height=>{
        console.log(normal, height);
        next(null, {normal, height});
      });
    }
  }

  function getHeightAt(tileProp,  callback){
    let filename = getTextureFilename({
      planetUUID: planet.uuid,
      textureType: 'height',
      lod: tileProp.lod, 
      tile: tileProp.tile, 
      face: tileProps.face
    });

    if(filesOpened[filename]){
      let buffer = filesOpened[filename];
      let height = readHeight(buffer, tileProps);
      callback(height);
    }


    readFileIfExists(filename,(err, content)=>{
      if (err) throw err;
      filesOpened[filename] = floatArray(content)
      let buffer = filesOpened[filename];
      let height = readHeight(buffer, tileProps);
      callback(height);
    })

  }

  function readHeight(buffer, tileProps){
    let {lod} = tileProps;
    let division = Math.pow(2, lod);
    let size = 1/division;

    let I = Math.floor(tileProps.inTileT / TextureSize)
    let J = Math.floor(tileProps.inTileS / TextureSize);
    let ix = J*division + I;
    return buffer[ix];
  }
}

function floatArray(buf){
  return new Float32Array(buf.buffer);
}

