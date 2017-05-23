import express from 'express';
import jpegJs from 'jpeg-js';
import path from 'path';
import fs from 'fs';
import zlib from 'zlib';
import bodyParser from 'body-parser';
import {getAtmosphereTextureFileName, getTextureFilename} from './utils.js';
import * as WorldManager from './WorldManager.js';
import * as TextureGenerator from './TextureCreator.js';

let app = express();


app.use(bodyParser.json()); 
app.use(bodyParser.urlencoded({ extended: false })); 

app.all('/*', (req,res, next)=>{
  res.setHeader("Access-Control-Allow-Origin", "*");
  next();
})

app.get('/get-list',(req,res)=>{
  WorldManager.getWorldList(content=>{
    res.end(content)
  });
})

app.post('/get-task-list-state', (req,res)=>{
  let uuids = JSON.parse(req.body.json)
  let states = TextureGenerator.getProcessStates(uuids);
  res.end(JSON.stringify(states));

})


app.post('/create-world', (req, res)=>{
  let object = JSON.parse(req.body.json);
  WorldManager.createWorld(object, end);
  function end(systemId){
    res.end(systemId);
  }
})


app.post('/generate-texture/',(req,res)=>{
  let object = JSON.parse(req.body.json)
  let {planetUUID} = object;
  WorldManager.retrievePlanetDescription(planetUUID, generateTexture(object))

  function generateTexture(params){
    return planet=>{
      TextureGenerator.generateTexture(planet, params, (err,uuid)=>{
        res.send(uuid)
      })

    }
  }
})

app.get('/texture/:planetUUID/:textureType/:lod/:face/:tile.raw',(req, res)=>{
  let file = getTextureFilename(req.params);
  fs.access(file, err=>{
    if(err) return res.status(404).send('not found');
    res.sendFile(file);
  })
});

app.get('/texture/:planetUUID/atmosphere/:resolution/:textureType.raw',(req, res)=>{
  let file = getAtmosphereTextureFileName(req.params);
  console.log("WTF", file);
  fs.access(file, err=>{
    if(err) return res.status(404).send('not found');
    res.sendFile(file);
  })
});

app.get('/', (req, res)=>{
  res.send("I'm working, it's fine");
})


app.listen(8082,()=>{
  console.log("listening...");
})


function monad(fn){
  return (err, ...args)=>{
    if(err) throw err;
    fn(...args);
  }
}


app.get('/bin-create/', (req,res) =>{
  let resolution = 2048;
  let ab = new Buffer(resolution*resolution * 4);

  for(let i = 0; i < resolution; ++i){
    for(let j = 0; j < resolution; ++j){
      let ix = (i*resolution + j)*4;
      let r = parseInt(i / resolution * 255);
      let g = parseInt(j / resolution * 255);
      ab[ix] = r;
      ab[ix+1] = g;
      ab[ix+2] = 0;
      ab[ix+3] = 0xFF;
    }
  }
  let gz = zlib.createGzip();
  zlib.deflate(ab, {level:9}, (err, buffer)=>{
    if(err) throw err;
    fs.writeFile(binFileName, buffer, err=>{
      if(err) throw err;
      res.send("bin - ok");
    });
  })
})

app.get('/read-test/', (req,res)=>{
  fs.readFile(binFileName, (err,data)=>{
    if(err) throw err;
    zlib.unzip(data, (err,buffer)=>{
      if(err) throw err;

      let image = jpegJs.encode({
        width:2048, 
        height:2048,
        data: buffer
      }, 100);

      fs.writeFile(goodFileName, image.data, err => {
        if(err) throw err;
          res.send('done');
      });

    })
  })
})


app.get('/jpg-create', (req, res)=>{
  let resolution = 8192;
  let ab = new Buffer(resolution*resolution * 4);

  for(let i = 0; i < resolution; ++i){
    for(let j = 0; j < resolution; ++j){
      let ix = (i*resolution + j)*4;
      let r = parseInt(i / resolution * 255);
      let g = parseInt(j / resolution * 255);
      ab[ix] = r;
      ab[ix+1] = g;
      ab[ix+2] = 0;
      ab[ix+3] = 0xFF;
    }
  }
  
  let image = {
    width: resolution,
    height: resolution,
    data: ab
  }
  let data = jpegJs.encode(image, 150);

  fs.writeFile(goodFileName, data.data, err => {
    if(err) throw err;
      res.send('done');
  });
})

app.get('/jpg', (req, res)=>{
  res.sendFile(goodFileName);
})
