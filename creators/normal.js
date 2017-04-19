import async from 'async';
import {getTileProps, calculateTileProperties, transformToCorrectFace, stToNormal} from './lookup.js'
import {getTextureFilename} from '../utils.js';
import {inflateFromIfExists} from './zipper.js';

let TextureSize = 512;

export function create(input, callback){
  let {planet, params} = input;
  let {lod, face, tile} = input.params;

  let division = Math.pow(2, lod);
  let T = Math.floor(tile / division);
  let S = tile % division;

  let s = S / division;
  let t = T / division;

  let filesOpened = {};
  let queryQueue = {};

  let thisTileProps = calculateTileProperties(face, lod, tile);

  let normalMap = new Buffer(TextureSize * TextureSize * 3);
  console.log('create normalmap' ,`lod:${lod}, face:${face}, tile:${tile}`);
  for(let i = 0; i < TextureSize; ++i){
    for(let j = 0; j < TextureSize; ++j){
      let vs = [[i-1, j],
       [i+1, j],

       [i, j-1],
       [i, j+1]]

       async.map(vs, lookupAt(thisTileProps), (err,results)=>{
         console.log("about to write this results", i,j, results);
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

      getHeightAt(getTileProps(coords, lod), height=>{
        next(null, {normal, height, i, j});
      });
    }
  }

  function getHeightAt(tileProp,  callback){
    let filename = getTextureFilename({
      planetUUID: planet.uuid,
      textureType: 'height',
      lod: tileProp.lod, 
      tile: tileProp.tile, 
      face: tileProp.face
    });

    enqueue(filename, (data)=>{
      let buffer = data;
      let height = readHeight(buffer, tileProp);
      return callback(height);
    });

    if(!filesOpened[filename]){

      filesOpened[filename] = {
        state:'queried',
        data:null
      }

      inflateFromIfExists(filename, (err, content)=>{
        if (err) {
          console.log('props', tileProp, filename, err);
          throw err;
        }
        console.log('one more file', filename);
        filesOpened[filename].data = floatArray(content)
        filesOpened[filename].state = 'readed';
        flushQueue(filename);
      });
    }else{
      flushQueue(filename);
    }
  }

  function enqueue(filename, fn){
    if(!queryQueue[filename])
      queryQueue[filename] = [];
    queryQueue[filename].push(fn);
    flushQueue(filename);
  }

  function flushQueue(filename){
    let o = filesOpened[filename];
    if(!o || o.state !== 'readed') {
      //console.log('WTF', o);
      return;
    }
    let data = o.data;
    let queue = queryQueue[filename];
    queryQueue[filename] = [];
    console.log('need to flush', queue.length);

    for(let i = 0; i < queue.length; ++i){
      queue[i](data);
    }
    console.log("flush done");
  }

  function readHeight(buffer, tileProps){
    let {lod} = tileProps;
    let division = Math.pow(2, lod);
    let size = 1/division;

    let I = Math.floor(tileProps.inTileT * TextureSize)
    let J = Math.floor(tileProps.inTileS * TextureSize);
    let ix = J*division + I;
    console.log('-----------------', J,I,  buffer[ix]);
    return buffer[ix];
  }
}

function floatArray(buf){
  console.log(buf.length);
  return new Float32Array(buf.buffer);
}

