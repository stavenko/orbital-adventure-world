function transformPosX(s,t, face){
    if(t >= 1) return {face: 1, t: s, s: 2-t } // face: -y
    if(t  < 0) return {face: 3, t: s, s: 1+t } // face: +y
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
    if(t >= 1) return {face: 4, t: 2-t, s: s} // face: +z
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
    if(t >= 1) return {face: 1, t: 2-t, s:1-s} // face: -z
    if(t  < 0) return {face: 3, t: -t , s:1-s} // face: +z
    if(s >= 1) return {face: 0, t: t, s: s-1} // face: -x
    if(s  < 0) return {face: 2, t: t, s: 1+s} // face: +x
    return {face, s,t}
}

export function transformToCorrectFace(s,t, face){
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
    return transformNegZ(s,t,face)
  }
  if(face == 5){
    return transformPosZ(s,t,face)
  }
}

export function switchFace(j,i, face, lod){
  let division = Math.pow(2, lod);
  let d = division;
  let d2 = 2*division;

  // J == S
  // I == T

  if(face == 2){
    return switchPosX()
  }
  if(face == 0){
    return switchNegX()
  }
  if(face == 1){
    return switchNegY()
  }
  if(face == 3){
    return switchPosY()
  }
  if(face == 4){
    return switchNegZ()
  }
  if(face == 5){
    return switchPosZ();
  }
  function switchPosX(){
    if(i >= division) return {face: 1, i: j, j: i % division } // face: -y
    if(i  < 0) return {face: 3, i: j, j: d+i } // face: +y
    if(j >= division) return {face: 5, i: i, j: j-d } // face: -z
    if(j  < 0) return {face: 4, i: i, j: d+j } // face: +z
    return {face, j,i}
  }
  function switchNegX(){
    if(i >= 1) return {face: 1, i: d-j, j: i-d } // face: -y
    if(i  < 0) return {face: 3, i: j,   j: d-i } // face: +y
    if(j >= 1) return {face: 4, i: i,   j: j-d } // face: -z
    if(j  < 0) return {face: 5, i: i,   j: d+j } // face: +z
    return {face, j,i}
  }
  function switchPosY(){
    if(i >= 1) return {face: 4, i: i%d, j: j} // face: +z
    if(i  < 0) return {face: 5, i: -i, j:d-j } // face: -z
    if(j >= 1) return {face: 2, i: j-d, j: d-i } // face: +x
    if(j  < 0) return {face: 0, i: -j, j: i  } // face: -x
    return {face, j,i}
  }
  function switchNegY(){
    if(i >= 1) return {face: 5, i: i%d, j:d-j} // face: -z
    if(i  < 0) return {face: 4, i: d+i, j:j} // face: +z
    if(j >= 1) return {face: 2, i: j%d, j:i} // face: +x
    if(j  < 0) return {face: 0, i: j+d, j:d-i} // face: -x
    return {face, j,i}
  }
  function switchPosZ(){
    if(i >= 1) return {face: 1, i: i-d, j:j} // face: -y
    if(i  < 0) return {face: 3, i: d+i, j:j} // face: +y
    if(j >= 1) return {face: 2, i: i, j:j-d} // face: +x
    if(j  < 0) return {face: 0, i: i, j:d+j} // face: -x
    return {face, j,i}
  }
  function switchNegZ(){
    if(i >= 1) return {face: 1, i: i%d, j:d-j} // face: -z
    if(i  < 0) return {face: 3, i: -i , j:d-j} // face: +z
    if(j >= 1) return {face: 0, i: i, j: j-d} // face: -x
    if(j  < 0) return {face: 2, i: i, j: d+j} // face: +x
    return {face, j,i}
  }
}

export function stToNormal(s,t, face){
  let ss = s * 2 - 1;
  let tt = t * 2 - 1;

  return [
    [-1, -tt, ss], //back
    [ss, -1, -tt], // left
    [1, -tt,-ss], // front
    [ss,  1, tt], // right
    [ss, -tt, 1], // top
    [-ss, -tt, -1], // bottom
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
  let inTileS = (s-(J*size))/size
  let inTileT = (t-(I*size))/size
  return { tile, face, lod, s: J*size, t:I*size, inTileS, inTileT,_coords:coords};
}

export function calculateTileProperties(face, lod, tile){
  let division = Math.pow(2, lod);
  let T = Math.floor(tile / division);
  let S = tile % division;

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
