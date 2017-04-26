function transformPosX(s,t, face){
    if(t >= 1) return {face: 1, t: s, s: 2-t } // face: -y
    if(t  < 0) {
      // console.log(">>>>", 1+t);
      return {face: 3, t: 1-s, s: 1+t }
    } // face: +y

    if(s >= 1) return {face: 5, t: t, s: s-1 } // face: -z
    if(s  < 0) return {face: 4, t: t, s: 1+s } // face: +z
    return {face, s,t}
}
function transformNegX(s,t, face){
    if(t >= 1) return {face: 1, t: 1-s, s: t-1 } // face: -y
    if(t  < 0) return {face: 3, t: s,   s: -t } // face: +y
    if(s >= 1) return {face: 4, t: t,   s: s-1 } // face: -z
    if(s  < 0) return {face: 5, t: t,   s: 1+s } // face: +z
    return {face, s,t}
}
function transformPosY(s,t, face){
    if(t >= 1) return {face: 4, t: 1-(2-t), s: s} // face: +z

    if(t  < 0) return {face: 5, t: -t, s:1-s } // face: -z
    if(s >= 1) return {face: 2, t: s-1, s: 1-t } // face: +x
    if(s  < 0) return {face: 0, t: -s, s: t  } // face: -x
    return {face, s,t}
}
function transformNegY(s,t, face){
    if(t >= 1) return {face: 5, t: 2-t, s:1-s} // face: -z
    if(t  < 0) return {face: 4, t: 1+t, s:s} // face: +z
    if(s >= 1) return {face: 2, t: 2-s, s:t} // face: +x
    if(s  < 0) return {face: 0, t: s+1, s:1-t} // face: -x
    return {face, s,t}
}
function transformPosZ(s,t, face){
    if(t >= 1) return {face: 1, t: t-1, s:s} // face: -y
    if(t  < 0) return {face: 3, t: 1+t, s:s} // face: +y

    if(s >= 1) return {face: 2, t: t, s:s-1} // face: +x
    if(s  < 0) return {face: 0, t: t, s:1+s} // face: -x
    return {face, s,t}
}

function transformNegZ(s,t, face){
    if(t >= 1) return {face: 1, t: 2-t, s:1-s} // face: -y
    if(t  < 0) return {face: 3, t: -t , s:1-s} // face: +y

    if(s >= 1) return {face: 0, t: t, s: s-1} // face: -x
    if(s  < 0) return {face: 2, t: t, s: 1+s} // face: +x
    return {face, s,t}
}

export function transformToCorrectFace(t, s, face){
  throw "NOOOO";
  if(face == 2){
    return transformPosX(s,t,face)
  }
  if(face == 0){
    return transformNegX(s,t,face)
  }
  if(face == 1){
    return transformNegY(s,t,face)
  }
  if(face == 3){
    return transformPosY(s,t,face)
  }
  if(face == 4){
    return transformPosZ(s,t,face)
  }
  if(face == 5){
    return transformNegZ(s,t,face)
  }
}


export function tileST2Normal(s,t,tileProps){
  

}

export function normal2FaceST(normal, lod){
  let coords = normalToST(normal);
  return getTileProps(coords);

}

export function normalToST(normal){
  let face = faceFromNormal(normal);
  let [s,t] = stFromNormal(normal, face).map(i=>0.5*(i+1));
  return {s,t,face};

}

function stFromNormal(n, f){
  let abs = Math.abs;
  if(f == 2) return [-n[2]/ abs(n[0]), -n[1]/abs(n[0])]; // +x
  if(f == 0) return [n[2]/ abs(n[0]),  -n[1]/abs(n[0])]; 

  if(f == 4) return [n[0]/ abs(n[2]),  -n[1]/abs(n[2])]; // +z
  if(f == 5) return [-n[0]/ abs(n[2]), -n[1]/abs(n[2])]; 

  if(f == 3) return [n[0]/ abs(n[1]),   n[2]/abs(n[1])]; // +y
  if(f == 1) return [n[0]/ abs(n[1]),  -n[2]/abs(n[1])]; 
  throw "undefined face " + f;
}

function faceFromNormal(n){
  let ax = Math.abs(n[0]);
  let ay = Math.abs(n[1]);
  let az = Math.abs(n[2]);

  if(ay >= ax && ay >= az){
    if(n[1] > 0.0) return 3;
    else return 1;
  } 

  if(ax >= ay && ax >= az){
    if(n[0] > 0.0) return 2;
    else return 0;
  }

  if(az >= ay && az >= ax){
    if(n[2] > 0.0) return 4;
    else return 5;
  }
  throw `incorrect data normal: ${n[0]}, ${n[1]}, ${n[2]}`; 
}


export function stToNormal(s,t, face){
  //if(s == 0) s+=1e-6;
  //if(t == 0) t+=1e-6;
  //if(s == 1) s-=1e-6;
  //if(t == 1) t-=1e-6;

  let ss = s * 2 - 1;
  let tt = t * 2 - 1;

  return [
    [-1, -tt, ss], // back

    [ss, -1, -tt], // left
    [1, -tt,-ss], // front
    [ss,  1, tt], // right
    [ss, -tt, 1], // top
    [-ss, -tt, -1] // bottom
  ][face];
}

export function getTileProps(coords, lod){
  let {s,t, face} = coords;
  if(s >= 1) s-=1e-10;
  if(t >= 1) t-=1e-10;
  let division = Math.pow(2,lod);
  let size = 1/division;
  let J = Math.floor(s/size);
  let I = Math.floor(t/size);
  let tile = J*division + I;
  let inTileS = (s-(J*size))/size;
  let inTileT = (t-(I*size))/size;
  return { tile, face, lod, s: J*size, t:I*size, inTileS, inTileT,_coords:coords};
}

export function calculateTileProperties(face, lod, tile){
  let division = Math.pow(2, lod);
  let S = Math.floor(tile / division);
  let T = tile % division;

  let s = S / division;
  let t = T / division;
  return {division, face, s,t, lod, tile, J:S, I:T}
}

export function calculateTile(face, lod, s,t){
  let division = Math.pow(2, lod);
  let S = s*division;
  let T = t*division;
  let tile = S*division + T;
  return {J:S, I:T, tile, lod}
}
