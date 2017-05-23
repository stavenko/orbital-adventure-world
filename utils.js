import path from 'path';
import config from './config.json';

const dir = config.rootDir;
const filesDir = path.join(dir, 'data');

export function getTextureFilename({planetUUID, textureType, lod, tile, face}){
  return path.join(filesDir, 'planets',planetUUID, 'textures', textureType, `${lod}`, `${face}`, `${tile}.raw`);
}

export function generateSystemId(){
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    var r = Math.random()*16|0, v = c == 'x' ? r : (r&0x3|0x8);
    return v.toString(16);
  });
}

export function getAtmosphereTextureFileName({planetUUID, resolution, textureType, resMu, resR, resMus, resNu}){

  if(!resolution)
    resolution = [resR, resMu, resMus, resNu].join('x'); 

  return path.join(filesDir, 'planets',planetUUID, 'textures/atmosphere',  resolution, `${textureType}.raw`); }

