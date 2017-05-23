import async from 'async';
import child_process from 'child_process';
import {ClassicalNoise} from './PerlinNoise.js';
import {generateSystemId} from './utils.js';
import isEqual from 'lodash/isEqual'

import * as HM from './creators/heightmap.js';
import * as NM from './creators/normal.js';
import * as AA from './creators/atmosphere.js';

const CREATORS = {
  'normal': NM,
  'height': HM,
  'irradianceTexture':AA,
  'inscatterTexture':AA,
  'deltaSRTexture':AA,
  'deltaSMTexture':AA,
  'deltaJTexture':AA,
  'deltaETexture':AA,
  'transmittanceTexture':AA
}

const TextureSize = 2048

const childProcess = child_process.fork('./worker.js');
const processState = {};
const processQueue = [];
const runningTasks = {};
let busy = false;



childProcess.on('message', function(event){
  if(event.type == 'completed'){
    busy = false;
    processState[event.uuid] = 'completed';
    delete runningTasks[event.uuid];
    tryPutNextTask();
  }
  if(event.type == 'postpone'){
    let tasklist = event.tasklist;
    processState[event.uuid] = 'queued'
    processQueue.push(...tasklist);
    busy = false;
    tryPutNextTask();
  }
});

export function getProcessStates(processes) {
  let result = {}
  for(let i =0 ;i <processes.length; ++i){
    let pr = processes[i];
    if(!processState[pr]) result[pr] = 'notfound';
    else{
      result[pr] = processState[pr];
      if(processState[pr] === 'completed')
        delete processState[pr]
    }
  }
  return result;
}

export function generateTexture(planetJSON, params, callback){
  let uuid = generateSystemId();
  processState[uuid] = 'queued';
  let planet = JSON.parse(planetJSON);
  //console.log(`generate query: ${params.textureType}, ${params.lod}, ${params.face}, ${params.tile} for planet:${planet.uuid}`);
  //console.log()
  let creator = CREATORS[params.textureType]
  let requirements = creator.getRequirements(params);
  if(requirements.length  !== 0){
    //console.log(`-----requirements for lod:${params.lod}, face:${params.face}, tile:${params.tile}------`);
    //console.log(requirements.map(req=>`lod:${req.lod}, face:${req.face}, tile:${req.tile}`).join('\n'));
    //console.log('------------------------------------------------------------------------------------');
    async.map(requirements, checkAndLaunchReq, (err, result)=>{
      if(err) throw err;
      //console.log("go! ------------------>");
      launch(params, callback);
    });
  }else{ 
    return launch(params, callback);
  }


  function checkAndLaunchReq(reqParams, next){
    creator.isExists(planet, reqParams, isExists=>{
      if(isExists) return next(null, "texture exists");
      generateTexture(planetJSON, reqParams, (err, uuid)=>{
        if(!err) next(null);
      })
    })
  }

  function launch(params, cb){
    let existUUID = alreadyInQueueOrRunning(planet, params)
    if(existUUID) {
      console.log("task already exists");
      return cb(null, existUUID);
    }
    console.log(`enqeued task: ${params.textureType}, lod:${params.lod}, face:${params.face}, tile:${params.tile} for planet:${planet.uuid}`);
    processQueue.push({planet, params, uuid});
    tryPutNextTask();
    cb(null, uuid);
  }

  function alreadyInQueueOrRunning(pl, pr){
    let ix = processQueue.findIndex(finder);
    if(ix !== -1) return processQueue[ix].uuid;
    let currentTasks = [];
    for(let k in runningTasks) currentTasks.push(runningTasks[k]);
    let task = currentTasks.find(finder);
    if(task) return task.uuid;

    function finder({planet, params}){
      let paramsList = []
      switch(params.textureType){
        case 'irradianceTexture':
        case 'inscatterTexture':
        case 'deltaSRTexture':
        case 'deltaSMTexture':
        case 'deltaJTexture':
        case 'deltaETexture':
        case 'transmittanceTexture':
          paramsList = ['textureType', 'resMu', 'resNu', 'resR', 'resMus']
          break;
        default:
          paramsList =['lod', 'face', 'tile', 'textureType'] ;
      }
      let paramsEq = paramsList.map(k=>params[k] == pr[k]).reduce((a,b)=>a&&b, true);
      console.log(planet.uuid === pl.uuid && paramsEq)
      return planet.uuid === pl.uuid && paramsEq
    }
  }

}

function tryPutNextTask(){
  if(busy || processQueue.length == 0) return;
  console.log("tasks left ", processQueue.length);

  let obj = processQueue.shift();
  busy = true;
  processState[obj.uuid] = 'running'
  runningTasks[obj.uuid] = obj;
  childProcess.send(obj);
}

function normalize(n){
  let l = Math.sqrt(n[0] * n[0] +  n[1] * n[1] + n[2] * n[2]);
  return n.map(i=>i/l);
}
