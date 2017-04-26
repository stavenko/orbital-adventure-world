import * as heightMap from './creators/heightmap.js';
import * as normalMap from './creators/normal.js';

process.on('message', function(event){
  let taskId = event.uuid;
  let {params} = event;

  if(params.textureType == 'height')
     return heightMap.create(event, err=>{
        if(err) {
          if(err.postpone){
            return process.send({
              uuid: taskId, 
              type:'postpone', 
              tasklist: err.postpone.tasklist
            })
          }
          throw err;
        }
        process.send({uuid: taskId, type:'completed'});
     })

  if(params.textureType == 'normal'){
     return normalMap.create(event, err=>{
        if(err) {
          if(err.postpone){
            console.log('postpone task');
            return process.send({
              uuid: taskId, 
              type:'postpone', 
              tasklist: err.postpone.tasklist
            })
          }
          throw err;
        }
        process.send({uuid: taskId, type:'completed'});
     })

  }

  process.send({uuid: taskId, type:'completed'});
})
