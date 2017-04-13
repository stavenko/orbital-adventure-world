import {create} from './creators/heightmap.js';

process.on('message', function(event){
  let taskId = event.uuid;
  let {params} = event;
  if(params.textureType != 'height') 
    process.send({uuid: taskId, type:'completed'});

  create(event, err=>{
    if(err) throw err;
    process.send({uuid: taskId, type:'completed'});
  })
})
