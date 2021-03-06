import async from 'async';
import {getTileProps, 
  calculateTileProperties, 
  normalToST, transformToCorrectFace, 
  stToNormal } from './lookup.js';
import {generateSystemId, getTextureFilename} from '../utils.js';
import {inflateFromIfExists} from './zipper.js';
import * as zipper from './zipper';
import fs from 'fs';

let TextureSize = 512;

const COLORS = [
  [255, 0,0],
  [255, 255,0],

  [0, 255,0],
  [0, 255,255], // 3

  [0, 0,255],
  [255, 0, 255],
]
let filesOpened = {};

function getTileShifts(params){
  let {lod, face, tile} = params;
  let division = Math.pow(2, lod);
  let size = 1.0 / division;
  let queryQueue = {};
  let dv = 1/TextureSize/division;

  let thisTileProps = calculateTileProperties(face, lod, tile);
  let tileShifts = [
    [thisTileProps.s-dv, thisTileProps.t],
    [thisTileProps.s-dv, thisTileProps.t+size+dv],
    [thisTileProps.s-dv, thisTileProps.t-dv],

    [thisTileProps.s+dv, thisTileProps.t+dv],
    [thisTileProps.s+dv, thisTileProps.t-dv],

    [thisTileProps.s, thisTileProps.t],
    [thisTileProps.s, thisTileProps.t-dv],
    [thisTileProps.s, thisTileProps.t+size+dv],

    [thisTileProps.s+size+dv, thisTileProps.t+size+dv],
    [thisTileProps.s+size+dv, thisTileProps.t],
    [thisTileProps.s+size, thisTileProps.t+dv],
    [thisTileProps.s+size+dv, thisTileProps.t-dv],

  ]

  console.log('====================TILESHIFTS=================');
  console.log(tileShifts);
  console.log('===============================================');
  return tileShifts;

}
export function isExists(forPlanet, params, callback){
  let filename = getTextureFilename({
    planetUUID: forPlanet.uuid,
    textureType: 'normal',
    lod: params.lod, 
    tile: params.tile, 
    face: params.face
  });
  fs.access(filename, err=>{
    if(!err) return callback(true);
    callback(false);
  })
}

export function getRequirements(params){
  let {lod, face, tile} = params;
  let thisTileProps = calculateTileProperties(face, lod, tile);

  let tiles = [];
  let indexesI = [0,1,2,3, 512-3,512-2,512-1];
  let indexesJ = [0,1,2,3, 512-3,512-2,512-1];
  let L = indexesI.length;
  for(let ii = 0; ii <L; ++ii){
    for(let jj = 0; jj <L; ++jj){
      let i = indexesI[ii];
      let j = indexesJ[jj];
      let vs = [
        [i-1, j],
        [i+1, j],
        [i, j-1],
        [i, j+1]
      ]
      for(let d =0; d< vs.length; ++d){
        let [di, dj] = vs[d];
        let normal = getLookupNormal(di, dj, thisTileProps); 
        let coords = normalToST(normal);
        let tp = getTileProps(coords, thisTileProps.lod);
        let ix = tiles.findIndex(({face, tile})=> {
          return face == tp.face && tile == tp.tile 
        });
        if(ix === -1) tiles.push({...tp, textureType:'height'});
      }
    }

  }
  return tiles;
}

export function create(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = input.params;

  let thisTileProps = calculateTileProperties(face, lod, tile);
  let normalMap = new Buffer(TextureSize * TextureSize * 3);
  async.map(getRequirements(params), getTextures(thisTileProps), (err, results)=>{
    if(err) throw err;

    console.log('create normalmap' ,`lod:${lod}, face:${face}, tile:${tile}`);
    let __t = Date.now();
    for(let i = 0; i < TextureSize; ++i){
      for(let j = 0; j < TextureSize; ++j){
        let vs = [
          [i-1, j],
          [i+1, j],
          [i, j-1],
          [i, j+1],
          [i, j]
        ];

        let heights = [];
        // let HHH = [];
        for(let d =0; d< vs.length; ++d){
          let [di, dj] = vs[d];
          let {normal, height} = lookupAt(di, dj, thisTileProps, d);
          normal = normalize(normal);
          let mheight = height + 10;
          let hgt = [normal[0] * mheight, normal[1] *mheight, normal[2] *mheight];
          heights.push(hgt);
        }
        let v1 = [
          heights[0][0] - heights[1][0],
          heights[0][1] - heights[1][1],
          heights[0][2] - heights[1][2]];

        let v2 = [
          heights[2][0] - heights[3][0],
          heights[2][1] - heights[3][1],
          heights[2][2] - heights[3][2]];

        let pnormal = normalize(cross((v1), (v2)));
        let direction = dot(pnormal, heights[4]);

        if(direction < 0){
          pnormal = [pnormal[0] * -1, pnormal[1]*-1, pnormal[2]*-1];
        }
        
        let ix = 3*(j * TextureSize + i);

        pnormal = ntc(pnormal);

        normalMap[ix]   = pnormal[0];
        normalMap[ix+1] = pnormal[1];
        normalMap[ix+2] = pnormal[2];
      }
    }
    console.log('time spent', Date.now() - __t);
    let filePath = getTextureFilename({
      planetUUID: planet.uuid,
      textureType: 'normal',
      lod: thisTileProps.lod, 
      tile: thisTileProps.tile, 
      face: thisTileProps.face
    })
    zipper.deflateTo(filePath, normalMap, ()=>{
      callback(null)
    });
  })


  function ntc([x,y,z]){
    return [(x+1)/2*255, (y+1)/2*255, (z+1)/2*255];
  }
  function getTextures(mainTileProps){
    return ({face, tile, textureType}, next)=>{
      let {lod} = mainTileProps;
      let filename = getTextureFilename({
        planetUUID: planet.uuid,
        textureType: 'height',
        lod: mainTileProps.lod, 
        tile: tile, 
        face: face
      });

      inflateFromIfExists(filename, (err, content)=>{
        if (err) return next(err);
        filesOpened[filename] = floatArray(content)
        next(null);
      })
    }
  }


  function lookupAt(i, j, tileProps, d){

    let normal = getLookupNormal(i, j, tileProps);
    let coords = normalToST(normal);
    let newTileProps = getTileProps(coords, tileProps.lod)
    let height = getHeightAt(newTileProps);
    return {normal, height, i, j}
  }

  function getHeightAt(tileProp, __ish){
    let filename = getTextureFilename({
      planetUUID: planet.uuid,
      textureType: 'height',
      lod: tileProp.lod, 
      tile: tileProp.tile, 
      face: tileProp.face
    });

    let data = filesOpened[filename];
    if(!data){
      let c = __ish.coords
      throw `Filebuffer for lod:${tileProp.lod} face:${tileProp.face} tile:${tileProp.tile} is not loaded 
      new normal: ${__ish.normal[0]}, ${__ish.normal[1]}, ${__ish.normal[2]}
      in new coords: s:${__ish.coords.s}, t:${__ish.coords.t}, face:${__ish.coords.face}, 
      let's get height from normal from: s:${__ish.s+__ish.ts}, t:${__ish.t + __ish.tt}, ts:${__ish.ts}, tt:${__ish.tt}, face:${__ish.face} ixj:${__ish.i}X${__ish.j}
      `;
    }
    return readHeight(data, tileProp);
  }

  function readHeight(buffer, tileProps){
    let {lod} = tileProps;
    let division = Math.pow(2, lod);
    let size = 1/division;

    let I = Math.floor(tileProps.inTileT * TextureSize)
    let J = Math.floor(tileProps.inTileS * TextureSize);
    let ix = J*TextureSize + I;
    
    try{
      let hv = buffer[ix];
      return   hv;
    }catch(e){
      console.log(tileProps, I, J);
    }
  }
}

function getLookupNormal(i, j, tileProps){
  let division = Math.pow(2, tileProps.lod);
  let tt = i / TextureSize / division;
  let ts = j / TextureSize / division;
  let {t,s,face} = tileProps;
  let normal = stToNormal(ts+s, tt+t, face);
  return normal
}

function floatArray(buf){
  return new Float32Array(buf.buffer);
}

function normalize(v){
  let l = Math.hypot(...v);
  return [v[0]/l, v[1]/l, v[2]/l];
}

function dot(v, u){
  return v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
}

function cross(u, v){
  let x = u[0],
      y = u[1],
      z = u[2];

  return [
    y*v[2] - z * v[1],
    z*v[0] - x * v[2],
    x*v[1] - y * v[0],
  ]
}
