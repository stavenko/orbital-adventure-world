import fs from 'fs';
import path from 'path';
import mkdirp from 'mkdirp'

export function writeFileInDir(filePath, content, cb){
  let dir = path.dirname(filePath);
  fs.access(dir, (err)=>{
    if(err) return mkdirp(dir, writeFile);
    writeFile(err);
  })

  function writeFile(err){
    if(err) throw err;
    fs.writeFile(filePath, content, cb);
  }
}
