import zlib from 'zlib';
import {writeFileInDir} from '../fsUtils.js';
import fs from 'fs';


export function inflateFromIfExists(path, callback){
  fs.access(path, (err)=>{
    if(err) return callback(err); 

    let gz = zlib.createGzip();
    fs.readFile(path, (err, content)=>{
      if(err) throw err;
      zlib.inflate(content, (err, buffer)=>{
        callback(err, buffer);
      })
    })
  });
}

export function inflateFrom(path, callback){
  let gz = zlib.createGzip();
  fs.readFile(path, (err, content)=>{
    if(err) throw err;
    zlib.inflate(content, (err, buffer)=>{
      callback(buffer);
    })
  })
}

export function deflateTo(filePath, buffer, callback){
  zlib.deflate(buffer, {level:9}, (err, buf)=>{
    if(err) throw err;
    writeFileInDir(filePath, buf, err=>{
      if(err) throw err;
      callback(null);
    });
  })
}
