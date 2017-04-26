import async from 'async';
import {getTileProps, calculateTileProperties, normalToST, transformToCorrectFace, stToNormal } from './lookup.js'
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
  return getTileShifts(params).map(([s,t])=>{
    let normal = stToNormal(s,t,params.face); 
    let coords = normalToST(normal);
    let {lod} = params;
    let nTile = getTileProps(coords, lod);

    return {face: nTile.face, tile: nTile.tile, lod, textureType:'height'};
  })
}

export function create(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = input.params;

  let thisTileProps = calculateTileProperties(face, lod, tile);
  let normalMap = new Buffer(TextureSize * TextureSize * 3);
  async.map(getTileShifts(params), getTextures(thisTileProps), (err, results)=>{
    if(err) throw err;

    console.log('create normalmap' ,`lod:${lod}, face:${face}, tile:${tile}`);
    let __t = Date.now();
    for(let i = 0; i < TextureSize; ++i){
      for(let j = 0; j < TextureSize; ++j){
        let vs = [[i-1, j],
         [i+1, j],
          

         [i, j-1],
         [i, j+1],
         [i, j]
        ];

        let heights = [];
        let HHH = [];
        for(let d =0; d< vs.length; ++d){
          let [di, dj] = vs[d];
          let {normal, height} = lookupAt(di, dj, thisTileProps, d);
          normal = normalize(normal);
          let mheight = height + 10;
          let hgt = [normal[0] * mheight, normal[1] *mheight, normal[2] *mheight];
          heights.push(hgt);
          HHH.push(height);
        }
        let v1 = [
          heights[0][0] - heights[1][0],
          heights[0][1] - heights[1][1],
          heights[0][2] - heights[1][2]];

        let v2 = [
          heights[2][0] - heights[3][0],
          heights[2][1] - heights[3][1],
          heights[2][2] - heights[3][2]];

        let H = HHH[4];
        let pnormal = normalize(cross(normalize(v1), normalize(v2)));
        let direction = dot(pnormal, heights[4]);

        if(direction < 0){
           pnormal = [pnormal[0] * -1, pnormal[1]*-1, pnormal[2]*-1];
        }
        
        let ix = 3*(j * TextureSize + i);
        let color = COLORS[thisTileProps.face];
        if(!color) console.log(H);


        pnormal = ntc(pnormal);

        if(false){
          if(false){
            normalMap[ix]   = Math.floor(255 *(H +1.0) / 2) ;
            normalMap[ix+1] = Math.floor(255 *(H +1.0) / 2) ;
            normalMap[ix+2] = Math.floor(255 *(H +1.0) / 2) ;
          }else{
            normalMap[ix]   = Math.floor(255 *(pnormal[0] +1.0) / 2) ;
            normalMap[ix+1] = Math.floor(255 *(pnormal[1] +1.0) / 2) ;
            normalMap[ix+2] = Math.floor(255 *(pnormal[2] +1.0) / 2) ;
          }
        }else{
          normalMap[ix]   = pnormal[0];
          normalMap[ix+1] = pnormal[1];
          normalMap[ix+2] = pnormal[2];
        }
      }
    }
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
    return ([s,t], next)=>{
      //------------------------------------------
      let normal = stToNormal(s,t, mainTileProps.face); 
      let coords = normalToST(normal);
      let {lod} = mainTileProps;
      let nTile = getTileProps(coords, lod);
      let fromFace = `|| face:${mainTileProps.face}, c.s=${coords.s}, c.t=${coords.t}`
      let normal__ = `N:[${normal[0]},${normal[1]},${normal[2]}]`;
      console.log(`openFile: lod ${lod}, face:${nTile.face}, tile:${nTile.tile}: ${fromFace}: ${normal__}`);

      let filename = getTextureFilename({
        planetUUID: planet.uuid,
        textureType: 'height',
        lod: mainTileProps.lod, 
        tile: nTile.tile, 
        face: nTile.face
      });

      inflateFromIfExists(filename, (err, content)=>{
        if (err) return next(err);
        filesOpened[filename] = floatArray(content)
        next(null);
      })
    }
  }

  function lookupAt(i, j, tileProps, d){
    let division = Math.pow(2, tileProps.lod);
    let tt = i / TextureSize / division;
    let ts = j / TextureSize / division;
    let {t,s,face} = tileProps;

    let normal = stToNormal(ts+s, tt+t, face);
    let coords = normalToST(normal);
    let newTileProps = getTileProps(coords, tileProps.lod)

    let height = getHeightAt(newTileProps, {...tileProps, coords, normal,ts, tt, i,j});
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

function cross([x,y,z], v){
  return [
    y*v[2] - z * v[1],
    z*v[0] - x * v[2],
    x*v[1] - y * v[0],
  ]
}
