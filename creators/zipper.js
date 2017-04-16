import zlib from 'zlib';
import {writeFileInDir} from '../fsUtils.js';
import fs from 'fs';


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
  zlib.deflate(buffer, {level:9}, (err, buffer)=>{
    if(err) throw err;
    writeFileInDir(filePath, buffer, err=>{
      if(err) throw err;
      callback(null);
    });
  })
}
