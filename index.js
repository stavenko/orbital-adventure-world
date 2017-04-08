import express from 'express';
import jpegJs from 'jpeg-js';
import path from 'path';
import fs from 'fs';
import zlib from 'zlib';

let app = express();
const dir = '/Users/vasilijstavenko/projects/';
const filesDir = path.join(dir, 'data');

const goodFileName = path.join(filesDir, 'synth.jpg'); 
const binFileName = path.join(filesDir, 'synth.raw'); 

app.get('/', (req, res)=>{
  res.send("I'm working, it's fine");
})


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
  // let file = path.join(dir, 'orbital-adventure-data/textures/Earth.png');
  res.sendFile(goodFileName);
})

app.listen(8082,()=>{
  console.log("listening...");
})

