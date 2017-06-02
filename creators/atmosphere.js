import async from 'async';
import fs from 'fs';
import {getAtmosphereTextureFileName} from '../utils.js';
import {generateAtmosphere} from './atmospehereGenerators/generator.js';
import * as zipper from './zipper';
export const atmosphereTextures=[
'transmittanceTexture',
'deltaIrradianceTexture',
'scatteringTexture',
'irradianceTexture'
]

export function getRequirements(params){
  return [];
}

export function isExists(forPlanet, params, callback){
  let p = {...params};
  p.resolution = ['resR', 'resMu', 'resMus', 'resNu'].map(k=>params[k]).join('x');
  let filename = getAtmosphereTextureFileName({
    planetUUID: forPlanet.uuid,
    ...p
  });
  fs.access(filename, err=>{
    console.log("exists", filename, err )
    if(!err) return callback(true);
    callback(false);
  })
}

export function create(event, callback){
  let {planet, params} = event;
  let {resolution, textureType} = params;
  isExists(planet, params, (e)=>{
    if(!e) return createAtmosphere(planet, params, callback);
    callback();
  })

}

function createAtmosphere(planet, params, callback){
  let generatedTextures = generateAtmosphere(params, planet)
  let keys = Object.keys(generatedTextures);

  console.log("actually create");
  async.map(keys, zip, (err,v)=>{
    if(err) throw err;
    console.log("err",err,v);
    setTimeout(()=>callback(), 5000);
  })

  function zip(key, next){
    let texture = generatedTextures[key].texture;
    let ps = {...params};
    ps.textureType = key;
    let filename = getAtmosphereTextureFileName({
      planetUUID: planet.uuid,
      ...ps
    })

    let buffer = Buffer.from(texture.buffer);

    console.log('write', filename);
    zipper.deflateTo(filename, buffer, next);
  };
}
