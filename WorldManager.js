import path from 'path';
import fs from 'fs';
import {writeFileInDir} from './fsUtils.js';
import async from 'async'

const dir = '/home/azl/projects/';
const filesDir = path.join(dir, 'data');

export function createWorld(worldCandidate, end){
  let systemId = null;
  

  if(worldCandidate.star){
    systemId = generateSystemId();
    for(let i =0; i < worldCandidate.planets.length; ++i){
      worldCandidate.planets[i].uuid = generateSystemId();
      worldCandidate.planets[i].systemId = systemId;
    }

    let starSystem = path.join(filesDir, '/systems/', `${systemId}.json`);
    writeFileInDir(starSystem, JSON.stringify(worldCandidate), createPlanets(worldCandidate.planets, systemId))
  }

  function createPlanets(planets, systemId){
    let planetsPath = path.join(filesDir, 'planets');
    let pp = p=>path.joun(planetsPath,p);
    return err =>{
      if(err) throw err;
      async.map(planets, _writer, (err, res)=>{
        end(systemId);
      })
    }

    function _writer(planet, cb){
      let file = path.join(planetsPath, planet.uuid, 'description.json');
      writeFileInDir(file, JSON.stringify(planet), cb);
    }
  }
}

export function getWorldList(cb){
  let sd = path.join(filesDir, 'systems');
  fs.readdir(sd, (err, files)=>{
    if(err) return cb('[]');
    async.map(files, (f, next)=>fs.readFile(path.join(sd,f),{encoding:'utf-8'}, next),(err, result)=>{
      result = result.map(s=>JSON.parse(s));
      cb(JSON.stringify(result));
    })
  })
}

export function retrievePlanetDescription(uuid, cb){
  
}

export function getTextureFilename({planetUUID, textureType, lod, tile}){
  return path.join(filesDir, 'planets', planetUUID, textureType, lod, `${tile}.raw`);
}

function generateSystemId(){
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    var r = Math.random()*16|0, v = c == 'x' ? r : (r&0x3|0x8);
    return v.toString(16);
  });
}
