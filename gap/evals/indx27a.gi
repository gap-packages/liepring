
BindGlobal( "IsSquareGF", function( F, a )
    local e;
    for e in Elements(F) do
        if e^2 = a then return e; fi;
    od;
    return false;
end );

#############################################################################
BindGlobal( "GetIndex2", function(P,W,A)
local F, x, y, z, t, u, v;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v <> Zero(F) then
  x:=x+(u^2*z+u*v-2*t*v*z-u*v*y)/(2*v^2);
  y:=y-u*z/v;
  t:=t/v-u^2/(4*v^2);
  if x <> Zero(F) then
    t:=t/x^2;
    return P^5+P^4*IntFFE(y)+P^2*IntFFE(t)+1;
  fi;
  if t = Zero(F) then return P^4*IntFFE(y)+1; fi;
  if IsSquareGF(GF(P),t)<>false then return P^4*IntFFE(y)+P^2+1; fi;
  return P^4*IntFFE(y)+P^2*W+1;
fi;

if z <> Zero(F) and u <> Zero(F) then
  x:=z*x+(t^2*z+t*u-t*u*y)*z/u^2;
  return P^5*IntFFE(x)+P^3+P;
fi;

if z <> Zero(F) and u = Zero(F) then
  if t <> Zero(F) then return P^3+P^2; fi;
  x:=x*z+y/2-y^2/4;
  return P^5*IntFFE(x)+P^3;
fi;

if u <> Zero(F) then
  t:=t-t*y+x*u;
  if t = Zero(F) then return P; fi;
  return P^2+P;
fi;

if t <> Zero(F) then return P^4*IntFFE(y)+P^2; fi;

if y <> One(F) or x = Zero(F) then return P^4*IntFFE(y); fi;

return P^5+P^4;

end );

#############################################################################
BindGlobal( "GetIndex3", function(P,W,A)
local F, x, y, z, t, u, v, c, e;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v <> Zero(F) then
  c:=z/v;
  x:=x-c*t;
  y:=y-c*u;
  t:=t/v;
  u:=u/v;
  if u <> Zero(F) then
    x:=x/u^2;
    y:=y/u;
    t:=t/u^2;
    return P^5*IntFFE(x)+P^4*IntFFE(y)+P^2*IntFFE(t)+P+1;
  fi;
  if y <> Zero(F) then
    x:=x/y^2;
    t:=t/y^2;
    return P^5*IntFFE(x)+P^4+P^2*IntFFE(t)+1;
  fi;
  if t <> Zero(F) then
    e:=t^-1;
    if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
    x:=x*e;
    t:=t*e;
    return P^5*IntFFE(x)+P^2*IntFFE(t)+1;
  fi;
  if x = Zero(F) then return 1; fi;
  if IsSquareGF(GF(P),x)<>false then return P^5+1; fi;
  return P^5*W+1;
fi;

if u <> Zero(F) then
  x:=x-t*y/u;
  t:=t/u;
  if t <> Zero(F) then
    x:=x/t^2;
    return P^5*IntFFE(x)+P^3*IntFFE(z)+P^2+P;
  fi;
  if x = Zero(F) then return P^3*IntFFE(z)+P; fi;
  if IsSquareGF(GF(P),x)<>false then return P^5+P^3*IntFFE(z)+P; fi;
  return P^5*W+P^3*IntFFE(z)+P;
fi;

if t <> Zero(F) then
  if y = Zero(F) then return P^3*IntFFE(z)+P^2; fi;
  return P^4+P^3*IntFFE(z)+P^2;
fi;

if y <> Zero(F) then
  x:=x/y^2;
  return P^5*IntFFE(x)+P^4+P^3*IntFFE(z);
fi;

if x = Zero(F) then return P^3*IntFFE(z); fi;
if IsSquareGF(GF(P),x)<>false then return P^5+P^3*IntFFE(z); fi;
return P^5*W+P^3*IntFFE(z);

end );

#############################################################################
BindGlobal( "GetIndex4", function(P,W,A)
local F, x, y, z, t, u, v, c, e;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v <> Zero(F) then
  x:=v^2*x+u*z^2-t*v*z-v*y*z;
  y:=z-u*z+v*y;
  t:=t*v-z-u*z+z^2;
  u:=u-2*z;
  return P^5*IntFFE(x)+P^4*IntFFE(y)+P^2*IntFFE(t)+P*IntFFE(u)+1;
fi;

if z = One(F) and u = -One(F) then
  if y <> Zero(F) then x:=x/y^2; t:=t/y; y:=One(F); fi;
  if y = Zero(F) and t <> Zero(F) then x:=x/t^2; t:=One(F); fi;
  if y+t <> Zero(F) or x = Zero(F) then 
     return P^4*IntFFE(y)+P^3+P^2*IntFFE(t)+P*IntFFE(u); fi;
  if y = One(F) then return P^5*x+P^4+P^3+P^2*IntFFE(t)+P*IntFFE(u); fi;
  if IsSquareGF(GF(P),x)<>false then return P^5+P^3+P*(P-1); fi;
  return P^5*W+P^3+P*(P-1);
fi;

if z <> One(F) and u = One(F)-2*z then
  x:=(t^2+2*y*t-4*x+4*x*z)/(4*z-4);
  if y <> Zero(F) then x:=x/y^2;
     return P^5*IntFFE(x)+P^4+P^3*IntFFE(z)+P*IntFFE(u); fi;
  if x = Zero(F) then return P^3*IntFFE(z)+P*IntFFE(u); fi;
  if IsSquareGF(GF(P),x)<>false then return P^5+P^3*IntFFE(z)+P*IntFFE(u); fi;
  return P^5*W+P^3*IntFFE(z)+P*IntFFE(u);
fi;

x:=(x*u^2-t*u*y+4*x*u*z-2*x*u-y^2*z+y^2-2*t*y*z+t*y+4*x*z^2-4*x*z+x)
   *(u+2*z-1)^-2;
t:=-(t+y-t*u-2*t*z+u*y)*(u+2*z-1)^-1;
if t <> Zero(F) then x:=x/t^2; 
  return P^5*IntFFE(x)+P^3*IntFFE(z)+P^2+P*IntFFE(u); fi;
if x = Zero(F) then return P^3*IntFFE(z)+P*IntFFE(u); fi;
if IsSquareGF(GF(P),x)<>false then return P^5+P^3*IntFFE(z)+P*IntFFE(u); fi;
return P^5*W+P^3*IntFFE(z)+P*IntFFE(u);

end );

#############################################################################
BindGlobal( "GetIndex5", function(P,W,A, curoots)
local F, x, y, z, t, u, v, c, e, a;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v <> Zero(F) then
  x:=x*v^3;
  y:=y*v^2;
  z:=z*v;
  t:=t*v;
  return P^5+P^4*IntFFE(z)+P^3*IntFFE(t)+P^2*IntFFE(y)+P*IntFFE(x)+IntFFE(u);
fi;
if z <> Zero(F) then
  x:=x/z^3;
  y:=y/z^2;
  t:=t/z;
  return P^4+P^3*IntFFE(t)+P^2*IntFFE(y)+P*IntFFE(x)+IntFFE(u);
fi;
if t <> Zero(F) then
  x:=x/t^3;
  y:=y/t^2;
  return P^3+P^2*IntFFE(y)+P*IntFFE(x)+IntFFE(u);
fi;
if y <> Zero(F) then
  e:=y^-1;
  if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
  a:=IsSquareGF(GF(P),e);
  x:=x*a^3;
  y:=y*e;
  if IntFFE(x) > (P-1)/2 then x:=-x; fi;
  return P^2*IntFFE(y)+P*IntFFE(x)+IntFFE(u);
fi;
if x = Zero(F) then return IntFFE(u); fi;
e:=x;
x:=One(F);
if curoots[IntFFE(e)] = 0 then e:=e/W; x:=W*One(F); fi;
if curoots[IntFFE(e)] = 0 then x:=x^2; fi;
return P*IntFFE(x)+IntFFE(u);
end );

#############################################################################
BindGlobal( "GetIndex6", function(P,W,A, curoots)
local F, x, y, z, t, u, v, c, e, x1, y1, z1, t1, u1, v1, a;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v <> Zero(F) then
  c:=z/v;
  x:=x-c*t;
  y:=y-c*u;
  z:=Zero(F);
  t:=t/v;
  u:=u/v;
  v:=One(F);
  if u <> Zero(F) then
    x:=x/u^3;
    y:=y/u^2;
    t:=t/u^2;
    u:=One(F);
    return P^5*IntFFE(x)+P^4*IntFFE(y)+P^2*IntFFE(t)+P+1;
  fi;
  if y <> Zero(F) then
    e:=y^-1;
    if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
    a:=IsSquareGF(GF(P),e);
    x:=x*a^3;
    y:=y*e;
    t:=t*e;
    if IntFFE(x) > (P-1)/2 then x:=-x; fi;
    return P^5*IntFFE(x)+P^4*IntFFE(y)+P^2*IntFFE(t)+1;
  fi;
  if t <> Zero(F) then
    e:=t^-1;
    if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
    a:=IsSquareGF(GF(P),e);
    x:=x*a^3;
    t:=t*e;
    if IntFFE(x) > (P-1)/2 then x:=-x; fi;
    return P^5*IntFFE(x)+P^2*IntFFE(t)+1;
  fi;
  if x = Zero(F) then return 1; fi;
  e:=x;
  x:=One(F);
  if curoots[IntFFE(e)] = 0 then e:=e/W; x:=W*One(F); fi;
  if curoots[IntFFE(e)] = 0 then x:=x^2; fi;
  return P^5*IntFFE(x)+1;
fi;
if u <> Zero(F) then
  c:=y/u;
  x:=x-c*t;
  y:=Zero(F);
  t:=t/u;
  u:=One(F);
  if z <> Zero(F) then
    x:=x/z^3;
    t:=t/z;
    return P^5*IntFFE(x)+P^3+P^2*IntFFE(t)+P;
  fi;
  if t <> Zero(F) then
    x:=x/t^3;
    return P^5*IntFFE(x)+P^2+P;
  fi;
  if x = Zero(F) then return P; fi;
  e:=x;
  x:=One(F);
  if curoots[IntFFE(e)] = 0 then e:=e/W; x:=W*One(F); fi;
  if curoots[IntFFE(e)] = 0 then x:=x^2; fi;
  return P^5*IntFFE(x)+P;
fi;

if t <> Zero(F) then
  if z <> Zero(F) then
    y:=y/z^2;
    return P^4*IntFFE(y)+P^3+P^2;
  fi;
  if y = Zero(F) then return P^2; fi;
  if IsSquareGF(GF(P),y)<>false then return P^4+P^2; fi;
  return P^4*W+P^2;
fi;

if z <> Zero(F) then
  x:=x/z^3;
  y:=y/z^2;
  return P^5*IntFFE(x)+P^4*IntFFE(y)+P^3;
fi;

if y <> Zero(F) then
  e:=y^-1;
  if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
  a:=IsSquareGF(GF(P),e);
  x:=x*a^3;
  y:=y*e;
  if IntFFE(x) > (P-1)/2 then x:=-x; fi;
  return P^5*IntFFE(x)+P^4*IntFFE(y);
fi;

if x = Zero(F) then return Zero(F); fi;
e:=x;
x:=One(F);
if curoots[IntFFE(e)] = 0 then e:=e/W; x:=W*One(F); fi;
if curoots[IntFFE(e)] = 0 then x:=x^2; fi;
return P^5*IntFFE(x);
end );

#############################################################################
BindGlobal( "GetIndex7", function(P,W,A)
local F, x, y, z, t, u, v, c, e, x1, y1, z1, t1, u1, v1;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v <> Zero(F) then
  c:=z/v;
  x:=x-c*t;
  y:=y-c*u;
  z:=Zero(F);
  t:=t/v;
  u:=u/v;
  v:=One(F);
  u1:=-u;
  x1:=-x;
  if [IntFFE(u1),IntFFE(x1)] < [IntFFE(u),IntFFE(x)] then u:=u1; x:=x1; fi;
  return P^5*IntFFE(x)+P^4*IntFFE(y)+P^2*IntFFE(t)+P*IntFFE(u)+1;
fi;

if u <> Zero(F) then
  c:=y/u;
  x:=x-c*t;
  y:=Zero(F);
  t:=t/u;
  u:=One(F);
  z1:=-z;
  t1:=-t;
  x1:=-x;
  if [IntFFE(z1),IntFFE(t1),IntFFE(x1)] < [IntFFE(z),IntFFE(t),IntFFE(x)] then 
     z:=z1; t:=t1; x:=x1; 
  fi;
  return P^5*IntFFE(x)+P^3*IntFFE(z)+P^2*IntFFE(t)+P;
fi;

if t <> Zero(F) then
  if IntFFE(z) > (P-1)/2 then z:=-z; fi;
  return P^4*IntFFE(y)+P^3*IntFFE(z)+P^2;
fi;

x1:=-x;
z1:=-z;
if [IntFFE(z1),IntFFE(x1)] < [IntFFE(z),IntFFE(x)] then x:=x1; z:=z1; fi;
return P^5*IntFFE(x)+P^4*IntFFE(y)+P^3*IntFFE(z);
end );

#############################################################################
BindGlobal( "GetIndex9", function(P,W,A,curoots)
local F, x, y, z, t, u, v, c, e, a;
F:=GF(P);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if z <> Zero(F) then
  x:=x/z^3;
  y:=y/z^2;
  t:=t/z^2;
  u:=u/z;
  return P^5+P^4*IntFFE(u)+P^3*IntFFE(y)+P^2*IntFFE(t)+P*IntFFE(x)+IntFFE(v);
fi;
if u <> Zero(F) then
  x:=x/u^3;
  y:=y/u^2;
  t:=t/u^2;
  return P^4+P^3*IntFFE(y)+P^2*IntFFE(t)+P*IntFFE(x)+IntFFE(v);
fi;
if y <> Zero(F) then
  e:=y^-1;
  if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
  a:=IsSquareGF(GF(P),e);
  x:=x*a^3;
  y:=y*e;
  t:=t*e;
  if IntFFE(x) > (P-1)/2 then x:=-x; fi;
  return P^3*IntFFE(y)+P^2*IntFFE(t)+P*IntFFE(x)+IntFFE(v);
fi;
if t <> Zero(F) then
  e:=t^-1;
  if IsSquareGF(GF(P),e)=false then e:=e*W; fi;
  a:=IsSquareGF(GF(P),e);
  x:=x*a^3;
  t:=t*e;
  if IntFFE(x) > (P-1)/2 then x:=-x; fi;
  return P^2*IntFFE(t)+P*IntFFE(x)+IntFFE(v);
fi;
if x = Zero(F) then return IntFFE(v); fi;
e:=x;
x:=One(F);
if curoots[IntFFE(e)] = 0 then e:=e/W; x:=W*One(F); fi;
if curoots[IntFFE(e)] = 0 then x:=x^2; fi;
return P*IntFFE(x)+IntFFE(v);
end );

#############################################################################
BindGlobal( "GetIndex10", function(P,W,A)
local x, y, z, t, u, v, c, e, x1, z1, u1;
x:=IntFFE(A[3][1]);
y:=IntFFE(A[3][2]);
z:=IntFFE(A[3][3]);
t:=IntFFE(A[4][1]);
u:=IntFFE(A[4][2]);
v:=IntFFE(A[4][3]);
x1:=(-x) mod P;
z1:=(-z) mod P;
u1:=(-u) mod P;
if [x1,z1,u1] < [x,z,u] then
  x:=x1; z:=z1; u:=u1;
fi;
return P^5*x+P^4*y+P^3*z+P^2*t+P*u+v;
end );

