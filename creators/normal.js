import async from 'async';
import {getTileProps, calculateTileProperties, transformToCorrectFace, stToNormal } from './lookup.js'
import {getTextureFilename} from '../utils.js';
import {inflateFromIfExists} from './zipper.js';
import * as zipper from './zipper';

let TextureSize = 512;

const COLORS = [
  [255, 0,0],
  [255, 255,0],

  [0, 255,0],
  [0, 255,255], // 3

  [0, 0,255],
  [255, 0, 255],
]

export function create(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = input.params;
  // if(face != 0) return callback();

  let division = Math.pow(2, lod);
  let size = 1.0 / division;
  let T = Math.floor(tile / division);
  let S = tile % division;

  let s = S / division;
  let t = T / division;

  let filesOpened = {};
  let queryQueue = {};
  let dv = 1e-6

  let thisTileProps = calculateTileProperties(face, lod, tile);
  let sts = [
    [thisTileProps.s-dv, thisTileProps.t],
    [thisTileProps.s + size + dv, thisTileProps.t],

    [thisTileProps.s+dv, thisTileProps.t+dv],

    [thisTileProps.s, thisTileProps.t-dv],
    [thisTileProps.s, thisTileProps.t+size+dv],
  ]
  let normalMap = new Buffer(TextureSize * TextureSize * 3);
  async.map(sts, getTextures(thisTileProps), (err, results)=>{
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
          let {normal, height} = lookupAt(di, dj, thisTileProps);
          normal = normalize(normal);
          height += 1;
          let hgt = [normal[0] * height, normal[1] *height, normal[2] *height];
          heights.push(hgt);
          HHH.push(height);
        }
        // console.log('----\n', heights, '\n-----');
        // let b =
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
        
        // pnormal = normalize(heights[4])
        let ix = 3*(j * TextureSize + i);
        if(i == 0) H = HHH[0];
        if(i == 511) H = HHH[1];
        if(j == 0) H = HHH[2];
        if(j == 511) H = HHH[3];
        let color = COLORS[Math.floor((H + 6)%6)];
        if(!color) console.log(H);

        if(false){
          if(true){
            normalMap[ix]   = Math.floor(255 *(H +1.0) / 2) ;
            normalMap[ix+1] = Math.floor(255 *(H +1.0) / 2) ;
            normalMap[ix+2] = Math.floor(255 *(H +1.0) / 2) ;
          }else{
            normalMap[ix]   = Math.floor(255 *(pnormal[0] +1.0) / 2) ;
            normalMap[ix+1] = Math.floor(255 *(pnormal[1] +1.0) / 2) ;
            normalMap[ix+2] = Math.floor(255 *(pnormal[2] +1.0) / 2) ;
          }
        }else{
          normalMap[ix]   = color[0];
          normalMap[ix+1] = color[1];
          normalMap[ix+2] = color[2];
        }
      }
    }
    console.log("time for preparation:", Date.now() -__t, normalMap);
    let filePath = getTextureFilename({
      planetUUID: planet.uuid,
      textureType: 'normal',
      lod: thisTileProps.lod, 
      tile: thisTileProps.tile, 
      face: thisTileProps.face
    })
    zipper.deflateTo(filePath, normalMap, callback);
  })


  function getTextures(mainTileProps){
    return ([s,t], next)=>{
      let coords = transformToCorrectFace(s,t, mainTileProps.face);
      let nTile = getTileProps(coords, mainTileProps.lod);

      // let tile = nFace.j*mainTileProps.division + nFace.i;
      let filename = getTextureFilename({
        planetUUID: planet.uuid,
        textureType: 'height',
        lod: mainTileProps.lod, 
        tile: nTile.tile, 
        face: nTile.face
      });
      console.log('request', filename, mainTileProps.face);
      inflateFromIfExists(filename, (err, content)=>{
        if (err) {
          console.log('props', mainTileProps, filename, err);
          next(err);
        }
        console.log('one more file', filename);
        filesOpened[filename] = floatArray(content)
        next(null);
      })

    }
  }


  function lookupAt(i, j, tileProps){
    let tt = i / TextureSize / division;
    let ts = j / TextureSize / division;
    let {t,s,face} = tileProps;
    let coords = {...transformToCorrectFace(s+ts, t+tt, face)};

    let normal = stToNormal(coords.s, coords.t, coords.face);

    let height = getHeightAt(getTileProps(coords, lod));
    return {normal, height, i, j}
  }

  function getHeightAt(tileProp){
    let filename = getTextureFilename({
      planetUUID: planet.uuid,
      textureType: 'height',
      lod: tileProp.lod, 
      tile: tileProp.tile, 
      face: tileProp.face
    });

    let data = filesOpened[filename];
    return readHeight(data, tileProp);
  }

  function readHeight(buffer, tileProps){
    let {lod} = tileProps;
    let division = Math.pow(2, lod);
    let size = 1/division;

    let I = Math.floor(tileProps.inTileT * TextureSize)
    let J = Math.floor(tileProps.inTileS * TextureSize);
    let ix = J*TextureSize + I;
    let hv = buffer[ix];
    return   hv
  }
}

function floatArray(buf){
  console.log(buf.length);
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
