import child_process from 'child_process';
import {ClassicalNoise} from './PerlinNoise.js';
import {generateSystemId} from './utils.js';

const TextureSize = 2048

const childProcess = child_process.fork('./worker.js');
const callbackMap = {}
const processState = {};
const processQueue = [];
let busy = false;



childProcess.on('message', function(event){
  if(event.type == 'completed'){
    busy = false;
    processState[event.uuid] = 'completed';
    tryPutNextTask();
  }
});

export function getProcessStates(processes){
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

export function generateTexture(planet, params, callback){
  let uuid = generateSystemId();
  callbackMap[uuid] = callback;
  processState[uuid] = 'queued';
  planet = JSON.parse(planet);
  processQueue.push({planet, params, uuid});
  tryPutNextTask();
  callback(null, uuid);
}

function tryPutNextTask(){
  if(busy || processQueue.length == 0) return;

  let obj = processQueue.shift();
  busy = true;
  processState[obj.uuid] = 'running'
  childProcess.send(obj);
}

function normalize(n){
  let l = Math.sqrt(n[0] * n[0] +  n[1] * n[1] + n[2] * n[2]);
  return n.map(i=>i/l);
}
