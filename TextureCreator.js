import child_process from 'child_process';
import {ClassicalNoise} from './PerlinNoise.js';
import {generateSystemId} from './utils.js';

const TextureSize = 2048

const childProcess = child_process.fork('./worker.js');
const callbackMap = {}
const processQueue = [];
let busy = false;

childProcess.on('message', function(event){
  // let data = event.data;
  if(event.type == 'completed'){
    busy = false;
    let cb = callbackMap[event.uuid];
    delete callbackMap[event.uuid];
    if(cb) cb(event);
    tryPutNextTask();
  }
});

export function generateTexture(planet, params, callback){
  let uuid = generateSystemId();
  callbackMap[uuid] = callback;
  planet = JSON.parse(planet);
  processQueue.push({planet, params, uuid});
  tryPutNextTask();
}

function tryPutNextTask(){
  if(busy || processQueue.length == 0) return;

  let obj = processQueue.shift();
  busy = true;
  childProcess.send(obj);
}

function normalize(n){
  let l = Math.sqrt(n[0] * n[0] +  n[1] * n[1] + n[2] * n[2]);
  return n.map(i=>i/l);
}
