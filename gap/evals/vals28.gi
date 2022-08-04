
BindGlobal( "ValsPreFunction28", function(P, case)
local W, F, SQ, i, CU, x, y, z, t, u, v, k, ii,  
transversal3, sofar, transversal4, transversal5, transversal6, params, 
xrange, yrange, zrange, urange, vrange, mats, y5, y6, y3, y4, A, y1, y2, 
a, b, c, m, B, D, C, z3, z4, z5, z6, z1, z2, new, index, ind1, AS, stab1, 
stab2, half, bcrange, sol, expect, val, B1, C1, B2, C2, BB1, BB2, lns, 
x1, u1, range, k2, km2, wover2, bc, s, A1, A2, AS1, AS2;

W:=PrimitiveRootMod(P);
F:=GF(P);
SQ:=Set(List([1..(P-1)/2], X -> (X^2) mod P));
for i in [2..P-1] do
  if not i in SQ then lns:=i; break; fi;
od;
CU:=[1];
for x in [2..P-1] do
  if (x^3) mod P = 1 then AddSet(CU,x); fi;
od;

##If p=1 mod 3 compute a transversal for the cube roots of 1
transversal3:=[];
if P mod 3 = 1 then
  transversal3:=[];
  sofar:=[];
  for i in [1..P-1] do
    if i^3 mod P in sofar then continue; fi;
    Add(transversal3,i);
    AddSet(sofar,i^3 mod P);
  od;
fi;

##If p=1 mod 4 then compute a transversal for the fourth roots of 1
transversal4:=[];
if P mod 4 = 1 then
  sofar:=[];
  for i in [1..(P-1)/2] do
    if i^4 mod P in sofar then continue; fi;
    Add(transversal4,i);
    AddSet(sofar,i^4 mod P);
  od;
fi;

##If p=1 mod 5 then compute a transversal for the fifth roots of 1
transversal5:=[];
if P mod 5 = 1 then
  sofar:=[];
  for i in [1..P-1] do
    if i^5 mod P in sofar then continue; fi;
    Add(transversal5,i);
    AddSet(sofar,i^5 mod P);
  od;
fi;

##If p=1 mod 3 then compute a transversal for the sixth roots of 1
transversal6:=[];
if P mod 3 = 1 then
  sofar:=[];
  for i in [1..P-1] do
    if i^6 mod P in sofar then continue; fi;
    Add(transversal6,i);
    AddSet(sofar,i^6 mod P);
  od;
fi;

params := List([1..24], x -> []);

##Case 1: cb=caa=cab=cac=0
if case = true or case = 1 then 

Add(params[1],[0,0,0,0,0,0]);
Add(params[1],[1,0,0,0,0,0]);
Add(params[1],[0,1,0,0,0,0]);
Add(params[1],[0,W,0,0,0,0]);
Add(params[1],[0,0,1,0,0,0]);
Add(params[1],[0,0,W,0,0,0]);
Add(params[1],[0,1,1,0,0,0]);
Add(params[1],[0,1,W,0,0,0]);
Add(params[1],[0,W,1,0,0,0]);
Add(params[1],[0,W,W,0,0,0]);
for y in [0..P-1] do
  Add(params[1],[1,y,1,0,0,0]);
  Add(params[1],[1,y,W,0,0,0]);
od;
Add(params[1],[1,1,0,1,0,0]);
Add(params[1],[1,W,0,1,0,0]);
for x in [0..P-1] do
  Add(params[1],[x,0,0,1,0,0]);
od;
Add(params[1],[0,0,0,0,1,0]);
Add(params[1],[0,1,0,0,1,0]);
Add(params[1],[0,W,0,0,1,0]);
Add(params[1],[0,0,0,1,1,0]);
Add(params[1],[0,1,0,1,1,0]);
Add(params[1],[0,W,0,1,1,0]);
Add(params[1],[1,0,0,0,0,1]);
Add(params[1],[0,0,0,0,0,1]);
Add(params[1],[0,0,1,0,0,1]);
Add(params[1],[0,0,W,0,0,1]);

if Length(params[1]) <> 3*P+22 then Error("params[1] wrong"); fi;
if case = 1 then return params[1]; fi;
fi;


##Case 2: caa=cab=cac=0, cb=baa
if case = true or case = 2 then 

Add(params[2],[0,0,0,0,0,0]);
Add(params[2],[1,0,0,0,0,0]);
Add(params[2],[0,1,0,0,0,0]);
Add(params[2],[0,W,0,0,0,0]);
Add(params[2],[0,0,1,0,0,0]);
Add(params[2],[0,0,W,0,0,0]);
Add(params[2],[0,1,1,0,0,0]);
Add(params[2],[0,1,W,0,0,0]);
Add(params[2],[0,W,1,0,0,0]);
Add(params[2],[0,W,W,0,0,0]);
for y in [0..P-1] do
  Add(params[2],[1,y,1,0,0,0]);
  Add(params[2],[1,y,W,0,0,0]);
od;
Add(params[2],[1,1,0,1,0,0]);
Add(params[2],[1,W,0,1,0,0]);
for x in [0..P-1] do
  Add(params[2],[x,0,0,1,0,0]);
od;
for x in [0..P-1] do
  Add(params[2],[0,x,0,0,1,0]);
  Add(params[2],[0,x,0,1,1,0]);
od;
Add(params[2],[0,0,0,0,0,1]);
Add(params[2],[1,0,0,0,0,1]);
if P mod 3 = 1 then
  Add(params[2],[W,0,0,0,0,1]);
  Add(params[2],[W^2 mod P,0,0,0,0,1]);
fi;
Add(params[2],[0,0,1,0,0,1]);
Add(params[2],[0,0,W,0,0,1]);
if P mod 4 = 1 then
  Add(params[2],[0,0,W^2 mod P,0,0,1]);
  Add(params[2],[0,0,W^3 mod P,0,0,1]);
fi;

if Length(params[2])<>5*P+13+Gcd(P-1,3)+Gcd(P-1,4) 
then Error("params[2] wrong"); fi;
if case = 2 then return params[2]; fi;
fi;

##Case 3: cb=bab=bac=cac=0
if case = true or case = 3 then 

Add(params[3],[0,0,0,0,0,0]);
Add(params[3],[1,0,0,0,0,0]);
Add(params[3],[0,0,0,0,1,0]);
Add(params[3],[0,1,0,0,1,0]);
Add(params[3],[0,0,0,0,1,1]);
Add(params[3],[1,0,0,0,1,1]);
Add(params[3],[0,0,0,0,1,W]);
Add(params[3],[1,0,0,0,1,W]);
Add(params[3],[0,0,1,0,0,1]);
Add(params[3],[0,0,W,0,0,W]);
for i in [1..P-1] do
  Add(params[3],[0,0,0,i,1,1]);
  Add(params[3],[0,0,0,i,1,W]);
od;
Add(params[3],[0,0,0,1,1,0]);
Add(params[3],[0,0,0,W,1,0]);
if P mod 4 = 1 then
  Add(params[3],[0,0,0,W^2 mod P,1,0]);
  Add(params[3],[0,0,0,W^3 mod P,1,0]);
fi;

if Length(params[3]) <> 2*P+8+Gcd(P-1,4) then Error("params[3] wrong"); fi;
if case = 3 then return params[3]; fi;
fi;

##Case 4: bab=bac=cac=0, cb=baa
if case = true or case = 4 then 

Add(params[4],[0,0,0,0,0,0]);
Add(params[4],[1,0,0,0,0,0]);
Add(params[4],[0,1,0,0,0,0]);
if P mod 3 = 1 then
  Add(params[4],[0,W,0,0,0,0]);
  Add(params[4],[0,W^2 mod P,0,0,0,0]);
fi;

Add(params[4],[0,0,0,0,1,0]);
Add(params[4],[0,1,0,0,1,0]);
if P mod 3 = 1 then
  Add(params[4],[0,W,0,0,1,0]);
  Add(params[4],[0,W^2 mod P,0,0,1,0]);
fi;

Add(params[4],[0,0,0,0,0,1]);
Add(params[4],[1,0,0,0,0,1]);
Add(params[4],[0,0,0,0,0,W]);
Add(params[4],[1,0,0,0,0,W]);

for y in [0..P-1] do
  Add(params[4],[0,0,1,0,0,y]);
  Add(params[4],[0,0,W,0,0,y]);
od;
Add(params[4],[0,0,1,0,1,1]);
Add(params[4],[0,0,W,0,1,W]);
for z in [1..(P-1)/2] do
  Add(params[4],[0,z,1,0,0,0]);
  Add(params[4],[0,z,W,0,0,0]);
od;

Add(params[4],[0,0,0,1,0,0]);
Add(params[4],[1,0,0,1,0,0]);
if P mod 5 = 1 then
  Add(params[4],[W,0,0,1,0,0]);
  Add(params[4],[W^2 mod P,0,0,1,0,0]);
  Add(params[4],[W^3 mod P,0,0,1,0,0]);
  Add(params[4],[W^4 mod P,0,0,1,0,0]);
fi;

Add(params[4],[0,0,0,1,1,0]);
Add(params[4],[0,0,0,1,W,0]);
if P mod 4 = 1 then
  Add(params[4],[0,0,0,1,W^2 mod P,0]);
  Add(params[4],[0,0,0,1,W^3 mod P,0]);
fi;
for x in [0..P-1] do
  Add(params[4],[0,0,0,1,x,1]);
  Add(params[4],[0,0,0,1,x,W]);
od;
for z in [1..(P-1)/2] do
  Add(params[4],[z,0,0,1,0,1]);
  Add(params[4],[z,0,0,1,0,W]);
od;

if Length(params[4]) <> 6*P+8+2*Gcd(P-1,3)+Gcd(P-1,4)+Gcd(P-1,5) then 
  Error("params[4] wrong"); 
fi;
if case = 4 then return params[4]; fi;
fi;

##Case 5: cb=bac=cac=0, caa=bab
if case = true or case = 5 then

Add(params[5],[0,0,0,0,0,0]);
Add(params[5],[1,0,0,0,0,0]);
Add(params[5],[0,1,0,0,0,0]);
Add(params[5],[0,W,0,0,0,0]);

Add(params[5],[0,0,1,0,0,0]);
Add(params[5],[0,0,W,0,0,0]);
Add(params[5],[0,1,1,0,0,0]);
Add(params[5],[0,1,W,0,0,0]);
Add(params[5],[0,W,1,0,0,0]);
Add(params[5],[0,W,W,0,0,0]);

for x in [0..P-1] do
  Add(params[5],[x,0,0,1,0,0]);
od;

Add(params[5],[0,0,0,0,1,0]);
Add(params[5],[0,1,0,0,1,0]);
Add(params[5],[0,W,0,0,1,0]);
if P mod 3 = 1 then
  Add(params[5],[0,W^2 mod P,0,0,1,0]);
  Add(params[5],[0,W^3 mod P,0,0,1,0]);
  Add(params[5],[0,W^4 mod P,0,0,1,0]);
  Add(params[5],[0,W^5 mod P,0,0,1,0]);
fi;

if P mod 4 = 1 then
for x in [0..(P-1)/2] do
  Add(params[5],[0,x,0,1,1,0]);
  Add(params[5],[0,x,0,W,1,0]);
  Add(params[5],[0,x,0,W^2 mod P,1,0]);
  Add(params[5],[0,x,0,W^3 mod P,1,0]);
od;
fi;

if P mod 4 = 3 then
for x in [0..P-1] do
  Add(params[5],[0,x,0,1,1,0]);
  Add(params[5],[0,x,0,W,1,0]);
od;
fi;

for x in [0..P-1] do
Add(params[5],[0,0,x,0,0,1]);
Add(params[5],[0,0,x,0,0,W]);
od;
Add(params[5],[1,0,0,0,0,1]);
Add(params[5],[1,0,0,0,0,W]);
Add(params[5],[1,0,1,0,0,1]);
Add(params[5],[1,0,W,0,0,W]);

if Length(params[5]) <> 5*P+13+2*Gcd(P-1,3)+Gcd(P-1,4) then 
  Error("params[5] wrong"); 
fi;
if case = 5 then return params[5]; fi;
fi;


##Case 6: bac=cac=0, caa=bab, cb=baa
if case = true or case = 6 then 

Add(params[6],[0,0,0,0,0,0]);
Add(params[6],[1,0,0,0,0,0]);
if P mod 5 = 1 then
  Add(params[6],[W,0,0,0,0,0]);
  Add(params[6],[W^2 mod P,0,0,0,0,0]);
  Add(params[6],[W^3 mod P,0,0,0,0,0]);
  Add(params[6],[W^4 mod P,0,0,0,0,0]);
fi;
Add(params[6],[0,1,0,0,0,0]);
Add(params[6],[0,W,0,0,0,0]);
if P mod 3 = 1 then
  Add(params[6],[0,W^2 mod P,0,0,0,0]);
  Add(params[6],[0,W^3 mod P,0,0,0,0]);
  Add(params[6],[0,W^4 mod P,0,0,0,0]);
  Add(params[6],[0,W^5 mod P,0,0,0,0]);
fi;
if P mod 4 = 3 then
  for x in [0..P-1] do
    Add(params[6],[0,x,1,0,0,0]);
    Add(params[6],[0,x,W,0,0,0]);
  od;
fi;
if P mod 4 = 1 then
  for x in [0..(P-1)/2] do
    Add(params[6],[0,x,1,0,0,0]);
    Add(params[6],[0,x,W,0,0,0]);
    Add(params[6],[0,x,W^2 mod P,0,0,0]);
    Add(params[6],[0,x,W^3 mod P,0,0,0]);
  od;
fi;
for x in [0..P-1] do
  Add(params[6],[x,0,0,1,0,0]);
  if P mod 5 = 1 then
    Add(params[6],[x,0,0,W,0,0]);
    Add(params[6],[x,0,0,W^2 mod P,0,0]);
    Add(params[6],[x,0,0,W^3 mod P,0,0]);
    Add(params[6],[x,0,0,W^4 mod P,0,0]);
  fi;
od;
if P mod 3 = 2 then
  for x in [0..P-1] do
  for y in [0..P-1] do
    Add(params[6],[0,y,0,x,1,0]);
  od;
  od;
fi;
if P mod 3 = 1 then
  xrange:=ShallowCopy(transversal3);
  Add(xrange,0);
  for x in xrange do
  for y in [0..P-1] do
    Add(params[6],[0,y,0,x,1,0]);
    Add(params[6],[0,y,0,x,W,0]);
    Add(params[6],[0,y,0,x,W^2 mod P,0]);
  od;
  od;
fi;
if P mod 4 = 3 then
  for y in [0..(P-1)/2] do
    Add(params[6],[y,0,0,0,0,1]);
    Add(params[6],[y,0,0,0,0,W]);
  od;
  for x in [1..P-1] do
    Add(params[6],[0,0,x,0,0,1]);
    Add(params[6],[0,0,x,0,0,W]);
  od;
  for y in [1..(P-1)/2] do
    Add(params[6],[y,0,1,0,0,1]);
    Add(params[6],[y,0,W,0,0,W]);
  od;
fi;
if P mod 4 = 1 then
  yrange:=ShallowCopy(transversal4);
  Add(yrange,0);
  for y in yrange do
    Add(params[6],[y,0,0,0,0,1]);
    Add(params[6],[y,0,0,0,0,W]);
    Add(params[6],[y,0,0,0,0,W^2 mod P]);
    Add(params[6],[y,0,0,0,0,W^3 mod P]);
  od;
  for x in [1..P-1] do
    Add(params[6],[0,0,x,0,0,1]);
    Add(params[6],[0,0,x,0,0,W]);
    Add(params[6],[0,0,x,0,0,W^2 mod P]);
    Add(params[6],[0,0,x,0,0,W^3 mod P]);
  od;
  for y in transversal4 do
    Add(params[6],[y,0,1,0,0,1]);
    Add(params[6],[y,0,W,0,0,W]);
    Add(params[6],[y,0,W^2 mod P,0,0,W^2 mod P]);
    Add(params[6],[y,0,W^3 mod P,0,0,W^3 mod P]);
  od;
fi;

if Length(params[6])<>P^2+3*P-3+(P+2)*Gcd(P-1,3)+(P+1)*Gcd(P-1,4)+(P+1)*Gcd(P-1,5)
then 
  Error("params[6] wrong"); 
fi;
if case = 6 then return params[6]; fi;
fi;


##Case 7: cb=baa=bac=cac=0
if case = true or case = 7 then 

for u in [0,1,W] do
for t in [0,1] do
for x in [0,1] do
  Add(params[7],[u,0,t,x,0,0]);
od;
od;
od;

for u in [0,1,W] do
for x in [0,1] do
  Add(params[7],[u,1,0,x,0,0]);
od;
od;

for x in [0,1,W] do
  Add(params[7],[0,1,1,x,0,0]);
od;

for u in [1,W] do
for x in [0..P-1] do
  Add(params[7],[u,1,1,x,0,0]);
od;
od;

for u in [0,1,W] do
for x in [0,1] do
for z in [1,W] do
  Add(params[7],[u,0,0,x,0,z]);
od;
od;
od;

for u in [0..P-1] do
for x in [0,1] do
for z in [1,W] do
  Add(params[7],[u,0,1,x,0,z]);
od;
od;
od;

for v in [0,1,W] do
for z in [0,1,W] do
  Add(params[7],[0,v,0,0,1,z]);
od;
od;

vrange:=[0,1,W];
if P mod 4 = 1 then vrange:=[0,1,W,W^2 mod P,W^3 mod P]; fi;
for v in vrange do
  Add(params[7],[0,v,0,1,1,0]);
od;

for v in [0..P-1] do
for z in [1,W] do
  Add(params[7],[0,v,0,1,1,z]);
od;
od;

for v in [0..P-1] do
for z in [0,1,W] do
  Add(params[7],[0,v,1,0,1,z]);
od;
od;

for v in [0..P-1] do
for x in [1,W] do
for z in [0..P-1] do
  Add(params[7],[0,v,1,x,1,z]);
od;
od;
od;

if Length(params[7]) <> 2*P^2+11*P+43+Gcd(P-1,4) then Error("params[7] wrong"); fi;
if case = 7 then return params[7]; fi;
fi;

##Case 8: cb=caa, baa=bac=cac=0
if case = true or case = 8 then 

urange:=[0,1,W];
if P mod 4 = 1 then urange:=[0,1,W,W^2 mod P,W^3 mod P]; fi;
for u in urange do
  Add(params[8],[u,0,0,0,0,0]);
  Add(params[8],[u,0,0,1,0,0]);
  Add(params[8],[u,1,0,0,0,0]);
od;

if P mod 3 = 1 then
  urange:=ShallowCopy(transversal3);
  Add(urange,0);
else;
  urange:=[0..P-1];
fi;
for u in urange do
  Add(params[8],[u,0,1,0,0,0]);
  Add(params[8],[u,0,1,1,0,0]);
  Add(params[8],[u,1,1,0,0,0]);
  if P mod 3 = 1 then
    Add(params[8],[u,0,W,0,0,0]);
    Add(params[8],[u,0,W,1,0,0]);
    Add(params[8],[u,1,W,0,0,0]);
    Add(params[8],[u,0,W^2 mod P,0,0,0]);
    Add(params[8],[u,0,W^2 mod P,1,0,0]);
    Add(params[8],[u,1,W^2 mod P,0,0,0]);
  fi;
od;

for u in [0..P-1] do
for t in [0..P-1] do
  Add(params[8],[u,1,t,1,0,0]);
od;
od;

for u in [0..P-1] do
for t in [0..(P-1)/2] do
for x in [0,1] do
for z in [1,W] do
  Add(params[8],[u,0,t,x,0,z]);
od;
od;
od;
od;

vrange:=[0,1,W];
if P mod 3 = 1 then vrange:=[0,1,W,W^2,W^3,W^4,W^5] mod P; fi;
for v in vrange do
  Add(params[8],[0,v,0,0,1,0]);
od;

if P mod 5 <> 1 then
  for v in [0..P-1] do
    Add(params[8],[0,v,0,1,1,0]);
  od;
fi;

if P mod 5 = 1 then
  vrange:=ShallowCopy(transversal5);
  Add(vrange,0);
  for v in vrange do
  for x in ([1,W,W^2,W^3,W^4] mod P) do
    Add(params[8],[0,v,0,x,1,0]);
  od;
  od;
fi;

if P mod 3 = 2 then
  for v in [0..P-1] do
  for x in [0..P-1] do
    Add(params[8],[0,v,1,x,1,0]);
  od;
  od;
fi;

if P mod 3 = 1 then
  xrange:=ShallowCopy(transversal3);
  Add(xrange,0);
  for v in [0..P-1] do
  for x in xrange do
    Add(params[8],[0,v,1,x,1,0]);
    Add(params[8],[0,v,W,x,1,0]);
    Add(params[8],[0,v,W^2 mod P,x,1,0]);
  od;
  od;
fi;

for v in [0..P-1] do
for x in [0..(P-1)/2] do
for z in [1,W] do
  Add(params[8],[0,v,0,x,1,z]);
od;
od;
od;

for v in [0..P-1] do
for t in [1..(P-1)/2] do
for x in [0..P-1] do
for z in [1,W] do
  Add(params[8],[0,v,t,x,1,z]);
od;
od;
od;
od;

if Length(params[8]) <> P^3+4*P^2+6*P+(P+5)*Gcd(P-1,3)+3*Gcd(P-1,4)+Gcd(P-1,5) 
then Error("params[8] wrong"); fi;
if case = 8 then return params[8]; fi;
fi;

#Cases 9 & 10: cb=bac=caa=0, cac=kbab (k=1,W)
if case = true or case = 9 or case = 10 then

Add(params[9],[0,0,0,0,0,0]);
Add(params[9],[1,0,0,0,0,0]);
Add(params[9],[0,1,0,0,0,0]);
Add(params[9],[0,W,0,0,0,0]);

for y in [0,1,W] do
  Add(params[9],[0,y,0,0,1,0]);
  if P mod 4 = 1 then
    Add(params[9],[0,y,0,0,W,0]);
   fi;
od;
yrange:=[0..P-1];
if P mod 4 = 1 then yrange:=[0..(P-1)/2]; fi;
for y in yrange do
  Add(params[9],[1,y,0,0,1,0]);
  if P mod 4 = 1 then
    Add(params[9],[1,y,0,0,W,0]);
  fi;
od;

for x in [0..(P-1)/2] do
  Add(params[9],[x,0,0,0,0,1]);
od;
Add(params[9],[0,1,0,0,0,1]);
Add(params[9],[0,W,0,0,0,1]);

for y in [0,1,W] do
for z in [0..(P-1)/2] do
  Add(params[9],[0,y,1,0,z,0]);
  Add(params[9],[0,y,W,0,z,0]);
od;
od;

for y in [0..P-1] do
for z in [0..(P-1)/2] do
  Add(params[9],[0,y,1,0,z,1]);
  Add(params[9],[0,y,W,0,z,1]);
od;
od;

for y in [0..P-1] do
for t in [0..(P-1)/2] do
  Add(params[9],[1,y,1,0,0,t]);
  Add(params[9],[1,y,W,0,0,t]);
od;
od;

for y in [0..P-1] do
for z in [1..(P-1)/2] do
for t in [0..P-1] do
  Add(params[9],[1,y,1,0,z,t]);
  Add(params[9],[1,y,W,0,z,t]);
od;
od;
od;

for x in [0..(P-1)/2] do
for y in [0..P-1] do
  Add(params[9],[y,0,0,1,0,x]);
od;
od;
for x in [0..(P-1)/2] do
  Add(params[9],[1,1,0,1,0,x]);
  Add(params[9],[1,W,0,1,0,x]);
od;

if P mod 4 = 1 then
  for x in [0..P-1] do
  for y in [0..(P-1)/2] do
    Add(params[9],[x,y,0,1,1,0]);
    Add(params[9],[x,y,0,1,W,0]);
  od;
  od;
fi;

if P mod 4 = 3 then
  for x in [0..P-1] do
  for y in [0..P-1] do
    Add(params[9],[x,y,0,1,1,0]);
  od;
  od;
fi;

if Length(params[9]) <> P^3+5*P^2/2+7*P+19/2+(P+4)*Gcd(P-1,4)/2
then Error("params[9] wrong"); fi;
if case = 9 then return params[9]; fi;
if case = 10 then return params[9]; fi;
if case = true then params[10] := params[9]; fi;
fi;

##Cases 11 & 12: bac=caa=0, cb=baa, cac=kbab (k=1,W)
if case = true or case = 11 or case = 12 then

Add(params[11],[0,0,0,0,0,0]);
Add(params[11],[1,0,0,0,0,0]);
if P mod 3 = 1 then
  Add(params[11],[W,0,0,0,0,0]);
  Add(params[11],[W^2 mod P,0,0,0,0,0]);
fi;
Add(params[11],[0,1,0,0,0,0]);
Add(params[11],[0,W,0,0,0,0]);
if P mod 4 = 1 then
  Add(params[11],[0,W^2 mod P,0,0,0,0]);
  Add(params[11],[0,W^3 mod P,0,0,0,0]);
fi;

if P mod 4 = 1 then
  xrange:=ShallowCopy(transversal4);
  Add(xrange,0);
  for x in xrange do
  for y in [0..P-1] do
    Add(params[11],[x,y,0,0,1,0]);
    Add(params[11],[x,y,0,0,W,0]);
  od;
  od;
fi;
if P mod 4 = 3 then
  for x in [0..(P-1)/2] do
  for y in [0..P-1] do
    Add(params[11],[x,y,0,0,1,0]);
  od;
  od;
fi;

if P mod 3 = 1 then
  for x in [0..(P-1)/2] do
    Add(params[11],[x,0,0,0,0,1]);
    Add(params[11],[x,0,0,0,0,W]);
    Add(params[11],[x,0,0,0,0,W^2 mod P]);
  od;
  for x in transversal3 do
    Add(params[11],[0,x,0,0,0,1]);
    Add(params[11],[0,x,0,0,0,W]);
    Add(params[11],[0,x,0,0,0,W^2 mod P]);
  od;
fi;
if P mod 3 = 2 then
  for x in [0..(P-1)/2] do
    Add(params[11],[x,0,0,0,0,1]);
  od;
  for x in [1..P-1] do
    Add(params[11],[0,x,0,0,0,1]);
  od;
fi;

for y in [0..P-1] do
for z in [0..(P-1)/2] do
for t in [0..(P-1)/2] do
  Add(params[11],[0,y,1,0,z,t]);
  Add(params[11],[0,y,W,0,z,t]);
od;
od;
od;
for y in [0..P-1] do
for x in [1..(P-1)/2] do
for z in [0..(P-1)/2] do
  Add(params[11],[x,y,1,0,z,0]);
  Add(params[11],[x,y,W,0,z,0]);
od;
od;
od;
for y in [0..P-1] do
for x in [1..(P-1)/2] do
for t in [1..(P-1)/2] do
for z in [0..P-1] do
  Add(params[11],[x,y,1,0,z,t]);
  Add(params[11],[x,y,W,0,z,t]);
od;
od;
od;
od;


if P mod 3 = 1 then
  for x in [0..P-1] do
    Add(params[11],[x,0,0,1,0,0]);
    Add(params[11],[x,0,0,W,0,0]);
    Add(params[11],[x,0,0,W^2 mod P,0,0]);
  od;
  for y in transversal3 do
    Add(params[11],[1,y,0,1,0,0]);
    Add(params[11],[W,y,0,W,0,0]);
    Add(params[11],[W^2 mod P,y,0,W^2 mod P,0,0]);
  od;
  for x in transversal6 do
  for y in [0..P-1] do
  for z in [0..P-1] do
    Add(params[11],[y,z,0,1,x,0]);
    Add(params[11],[y,z,0,W,x,0]);
    Add(params[11],[y,z,0,W^2 mod P,x,0]);
  od;
  od;
  od;
  for x in [0..P-1] do
  for z in [1..(P-1)/2] do
    Add(params[11],[x,0,0,1,0,z]);
    Add(params[11],[x,0,0,W,0,z]);
    Add(params[11],[x,0,0,W^2 mod P,0,z]);
  od;
  od;
  for y in transversal3 do
  for z in [1..(P-1)/2] do
    Add(params[11],[1,y,0,1,0,z]);
    Add(params[11],[W,y,0,W,0,z]);
    Add(params[11],[W^2 mod P,y,0,W^2 mod P,0,z]);
  od;
  od;
fi;

if P mod 3 = 2 then
  for x in [0..P-1] do
    Add(params[11],[x,0,0,1,0,0]);
  od;
  for y in [1..P-1] do
    Add(params[11],[1,y,0,1,0,0]);
  od;
  for x in [1..(P-1)/2] do
  for y in [0..P-1] do
  for z in [0..P-1] do
    Add(params[11],[y,z,0,1,x,0]);
  od;
  od;
  od;
  for x in [0..P-1] do
  for z in [1..(P-1)/2] do
    Add(params[11],[x,0,0,1,0,z]);
  od;
  od;
  for y in [1..P-1] do
  for z in [1..(P-1)/2] do
    Add(params[11],[1,y,0,1,0,z]);
  od;
  od;
fi;

if Length(params[11]) <> 
(P^4+P^3+4*P^2+P-1+(P^2+2*P+3)*Gcd(P-1,3)+(P+2)*Gcd(P-1,4))/2
then Error("params[11] wrong"); fi;
if case = 11 then return params[11]; fi;
if case = 12 then return params[11]; fi;
if case = true then params[12] := params[11]; fi;
fi;


##Case 13, cb=bac=0, caa=baa, cac=-bab
if case = true or case = 13 then

##pb=pc=0
Add(params[13],[0,0,0,0,0,0]);
Add(params[13],[1,0,0,0,0,0]);
Add(params[13],[0,1,0,0,0,0]);

##pb=0, pc=baa or wbaa
Add(params[13],[0,0,0,0,1,0]);
Add(params[13],[0,0,0,0,W,0]);
Add(params[13],[0,1,0,0,1,0]);
Add(params[13],[0,1,0,0,W,0]);
Add(params[13],[0,W,0,0,1,0]);
Add(params[13],[0,W,0,0,W,0]);
for y in [0..P-1] do
  Add(params[13],[1,y,0,0,1,0]);
  Add(params[13],[1,y,0,0,W,0]);
od;

##pb=0, pc=bab
for x in [0..P-1] do
  Add(params[13],[x,0,0,0,0,1]);
od;
Add(params[13],[-1,1,0,0,0,1] mod P);
Add(params[13],[-1,W,0,0,0,1] mod P);

##pb=pc=bab
Add(params[13],[0,0,0,1,0,1]);
Add(params[13],[1,0,0,1,0,1]);
Add(params[13],[0,1,0,1,0,1]);

##pb=bab, pc=-bab
for x in [0..P-1] do
  Add(params[13],[x,0,0,1,0,-1] mod P);
od;
Add(params[13],[2,1,0,1,0,-1] mod P);

##pb=bab, pc=baa or wbaa
for x in [0..P-1] do
for y in [0..P-1] do
  Add(params[13],[x,y,0,1,1,0]);
  Add(params[13],[x,y,0,1,W,0]);
od;
od;

##pb=pc=baa
Add(params[13],[0,0,1,0,1,0]);
Add(params[13],[1,0,1,0,1,0]);
Add(params[13],[0,1,1,0,1,0]);
Add(params[13],[1,1,1,0,1,0]);

##pb=baa, pc=baa+bab
for x in [0..P-2] do ##Note that -1 is excluded
  Add(params[13],[x,0,1,0,1,1]);
  y:=((x+1)/2)*One(F); y:=IntFFE(y);
  Add(params[13],[x,y,1,0,1,1]);
od;
Add(params[13],[-1,0,1,0,1,1] mod P);
Add(params[13],[-1,1,1,0,1,1] mod P);

##pb=baa, pc=-baa
Add(params[13],[0,0,1,0,-1,0] mod P);
Add(params[13],[0,1,1,0,-1,0] mod P);
Add(params[13],[0,W,1,0,-1,0] mod P);
for y in [0..P-1] do
  Add(params[13],[1,y,1,0,-1,0] mod P);
od;

##pb=baa, pc=-baa+bab
half:=(P+1)/2;
Add(params[13],[1,-half,1,0,-1,1] mod P);
Add(params[13],[1,half,1,0,-1,1] mod P);
Add(params[13],[1,W-half,1,0,-1,1] mod P);
if P mod 4 = 1 then
  Add(params[13],[1,W^2-half,1,0,-1,1] mod P);
  Add(params[13],[1,W^3-half,1,0,-1,1] mod P);
fi;
for y in [0..P-1] do
  Add(params[13],[0,y,1,0,-1,1] mod P);
  Add(params[13],[1-W,y,1,0,-1,1] mod P);
od;

##pb=pc=wbaaa
Add(params[13],[0,0,W,0,W,0]);
Add(params[13],[1,0,W,0,W,0]);
Add(params[13],[0,1,W,0,W,0]);
Add(params[13],[1,1,W,0,W,0]);

##pb=wbaa, pc=wbaa+bab
for x in [0..P-2] do ##Note that -1 is excluded
  Add(params[13],[x,0,W,0,W,1]);
  y:=One(F)*((x+1)/(2*W)); y:=IntFFE(y);
  Add(params[13],[x,y,W,0,W,1]);
od;
Add(params[13],[-1,0,W,0,W,1] mod P);
Add(params[13],[-1,1,W,0,W,1] mod P);

if Length(params[13]) <> 2*P^2+11*P+27+Gcd(P-1,4)
then Error("params[13] wrong"); fi;
if case = 13 then return params[13]; fi;
fi;

##Case 14, bac=0, cb=caa=baa, cac=-bab
if case = true or case = 14 then 

##1. v <> +/- y 
for y in [1,W] do
for u in [0..P-1] do
for x in [0..(P-1)/2] do
  zrange:=[0..P-1];
  if x = 0 then zrange:=[0..(P-1)/2]; fi;
  for z in zrange do
    Add(params[14],[0,u,0,x,y,z]);
  od;
od;
od;
od;

for u in [1,W] do
for x in [0..(P-1)/2] do
  zrange:=[0..P-1];
  if x = 0 then zrange:=[0..(P-1)/2]; fi;
  for z in zrange do
    Add(params[14],[0,u,1,x,-1,z] mod P);
  od;
od;
od;

for x in [0..P-1] do
for z in [0..P-1] do
  new:=true;
  for a in [1..P-1] do
    x1:=One(F)*((x+z+x*a^2-z*a^2)/(2*a^3)); x1:=IntFFE(x1);
    z1:=One(F)*((x+z-x*a^2+z*a^2)/(2*a^3)); z1:=IntFFE(z1);
    if [x1,z1] < [x,z] then new:=false; break; fi;
  od;
  if new then Add(params[14],[0,0,1,x,-1,z] mod P); fi;
od;
od;



for v in [1,W] do
for x in [0..P-1] do
for z in [0..P-1] do
  if x = (z+1) mod P then continue; fi;
  x1:=(-x) mod P; z1:=(-z) mod P;
  if [z1,x1] < [x,z] then continue; fi;
  Add(params[14],[1,0,v,x,v,z]);
od;
od;
od;

for x in [0..P-1] do
for z in [0..P-1] do
  if x = (z+1) mod P then continue; fi;
  new:=true;
  for a in [1..P-1] do
    x1:=One(F)*((x+z+x*a^3-z*a^3)/(2*a^3)); x1:=IntFFE(x1);
    z1:=One(F)*((x+z-x*a^3+z*a^3)/(2*a^3)); z1:=IntFFE(z1);
    if [x1,z1] < [x,z] then new:=false; break; fi;
  od;
  if new then Add(params[14],[1,0,0,x,0,z]); fi;
od;
od;

for v in [1,W] do
for u in [0..(P-1)/2] do
  Add(params[14],[1,u,v,0,v,-1] mod P);
od;
od;

for z in [0..P-1] do
  Add(params[14],[1,1,0,z+1,0,z] mod P);
od;

for z in [0..P-1] do
  new:=true;
  for a in [1..P-1] do
    z1:=One(F)*((-a^3+2*z+1)/(2*a^3)); z1:=IntFFE(z1);
    if z1 < z then new:=false; break; fi;
  od;
  if new then Add(params[14],[1,0,0,z+1,0,z] mod P); fi;
od;

for v in [1,W] do
for z in [1..(P-1)/2] do
  Add(params[14],[0,0,v,0,v,z]);
od;
od;

Add(params[14],[0,0,0,0,0,1]);
if P mod 3 = 1 then
  Add(params[14],[0,0,0,0,0,W]);
  Add(params[14],[0,0,0,0,0,W^2 mod P]);
fi;

Add(params[14],[0,0,0,1,0,-1] mod P);
Add(params[14],[0,0,W,1,W,-1] mod P);
Add(params[14],[0,0,W^2,1,W^2,-1] mod P);

for v in [1,W] do
for u in [0,1] do
  Add(params[14],[0,u,v,0,v,0]);
od;
od;

xrange:=[0,1];
if P mod 3 = 1 then xrange:=[0,1,W,W^2 mod P]; fi;
for u in [0,1] do
for x in xrange do
  Add(params[14],[0,u,0,x,0,x]);
od;
od;

if Length(params[14]) <> P^3+2*P^2+6*P+10+(P+4)*Gcd(P-1,3)
then Error("params[14] wrong"); fi;
if case = 14 then return params[14]; fi;
fi;

##Case 15, cb=baa=bac=caa=0
if case = true or case = 15 then

Add(params[15],[0,0,0,0,0,0]);
Add(params[15],[0,0,0,0,1,0]);
Add(params[15],[0,0,0,1,1,0]);
if P mod 3 = 1 then
  Add(params[15],[0,0,0,1,W,0]);
fi;

Add(params[15],[0,1,0,0,0,0]);
Add(params[15],[0,W,0,0,0,0]);
Add(params[15],[0,1,0,0,1,0]);
Add(params[15],[0,W,0,0,1,0]);
Add(params[15],[0,1,0,1,0,0]);
Add(params[15],[0,W,0,1,0,0]);
Add(params[15],[0,1,0,1,1,0]);
Add(params[15],[0,W,0,1,1,0]);
if  P mod 3 = 1 then
  Add(params[15],[0,1,0,1,W,0]);
  Add(params[15],[0,W,0,1,W,0]);
  Add(params[15],[0,1,0,1,W^2 mod P,0]);
  Add(params[15],[0,W,0,1,W^2 mod P,0]);
fi;

Add(params[15],[1,1,0,0,0,0]);
Add(params[15],[W,W,0,0,0,0]);
Add(params[15],[1,1,0,0,1,0]);
Add(params[15],[W,W,0,0,1,0]);
for u in [1..(P-1)/2] do
  u1:=One(F)*u^-1; u1:=IntFFE(u1);
  if u1 < u or (P-u1) < u then continue; fi;
  Add(params[15],[1,1,0,1,u,0]);
  Add(params[15],[W,W,0,1,u,0]);
od;

Add(params[15],[1,W,0,0,0,0]);
Add(params[15],[1,W,0,0,1,0]);
Add(params[15],[1,W,0,1,0,0]);
for u in [1..(P-1)/2] do
  Add(params[15],[1,W,0,1,u,0]);
od;

Add(params[15],[0,0,0,0,0,1]);
Add(params[15],[0,0,0,0,1,1]);
Add(params[15],[0,0,0,0,W,1]);
Add(params[15],[0,1,0,0,0,1]);
Add(params[15],[0,1,0,0,1,1]);
Add(params[15],[0,1,0,0,W,1]);
Add(params[15],[0,W,0,0,0,1]);
Add(params[15],[0,W,0,0,1,1]);
Add(params[15],[0,W,0,0,W,1]);
Add(params[15],[1,0,0,0,0,1]);
Add(params[15],[1,0,0,0,1,1]);
Add(params[15],[1,0,0,0,W,1]);
Add(params[15],[W,0,0,0,0,1]);
Add(params[15],[W,0,0,0,1,1]);
Add(params[15],[W,0,0,0,W,1]);
for x in [1,W] do
for y in [1,W] do
for u in [0..P-1] do
  Add(params[15],[x,y,0,0,u,1]);
od;
od;
od;

for x in [1,W] do
for y in [0..P-1] do
for u in [0..P-1] do
  Add(params[15],[x,y,0,1,u,1]);
od;
od;
od;
for y in [0,1,W] do
for u in [0..P-1] do
  Add(params[15],[0,y,0,1,u,1]);
od;
od;

for t in [0..P-2] do
for u in [t+1..P-1] do
  Add(params[15],[0,0,1,t,u,1]);
  Add(params[15],[0,1,1,t,u,1]);
  Add(params[15],[0,W,1,t,u,1]);
  for y in [0..P-1] do
    Add(params[15],[1,y,1,t,u,1]);
    Add(params[15],[W,y,1,t,u,1]);
  od;
od;
od;

range:=[];
for k in [1..(P-1)/2] do
  k2:=k^2 mod P;
  km2:=One(F)*(k2^-1); km2:=IntFFE(km2);
  if k2 <= km2 then AddSet(range,k2); fi;
od;
for t in [0..P-1] do
   Add(params[15],[0,0,1,t,t,1]);
   Add(params[15],[0,1,1,t,t,1]);
   Add(params[15],[0,W,1,t,t,1]);
   for k2 in range do
      Add(params[15],[1,k2,1,t,t,1]);
      Add(params[15],[W,W*k2,1,t,t,1] mod P);
   od;
   for k in [1..(P-1)/2] do
      Add(params[15],[1,W*k^2,1,t,t,1] mod P);
   od;   
od;

if Length(params[15]) <> P^3+(7*P^2+17*P+59+5*Gcd(P-1,3)+(P+1)*Gcd(P-1,4))/2
then Error("params[15] wrong"); fi;
if case = 15 then return params[15]; fi;
fi;

##Case 16, cb=bac=caa=0, baa=cac
if case = true or case = 16 then 

for y in [0,1,W] do
  Add(params[16],[0,y,0,0,0,0]);
  Add(params[16],[0,y,0,1,0,0]);
  Add(params[16],[0,y,0,W,0,0]);
  for t in [0..P-1] do
    Add(params[16],[1,y,0,t,0,0]);
    Add(params[16],[W,y,0,t,0,0]);
  od;
od;

Add(params[16],[0,0,1,0,0,0]);
Add(params[16],[0,0,W,0,0,0]);
Add(params[16],[0,1,1,0,0,0]);
Add(params[16],[0,1,W,0,0,0]);
Add(params[16],[0,W,1,0,0,0]);
Add(params[16],[0,W,W,0,0,0]);
if P mod 4 = 1 then
  Add(params[16],[0,W^2 mod P,1,0,0,0]);
  Add(params[16],[0,W^2 mod P,W,0,0,0]);
  Add(params[16],[0,W^3 mod P,1,0,0,0]);
  Add(params[16],[0,W^3 mod P,W,0,0,0]);
fi;
for y in [0..P-1] do
  Add(params[16],[0,y,1,1,0,0]);
  Add(params[16],[0,y,W,1,0,0]);
  Add(params[16],[0,y,1,W,0,0]);
  Add(params[16],[0,y,W,W,0,0]);
od;
for y in [0..P-1] do
for t in [0..P-1] do
  Add(params[16],[1,y,1,t,0,0]);
  Add(params[16],[1,y,W,t,0,0]);
  Add(params[16],[W,y,1,t,0,0]);
  Add(params[16],[W,y,W,t,0,0]);
od;
od;

Add(params[16],[0,0,0,0,0,1]);
Add(params[16],[0,0,1,0,0,1]);
Add(params[16],[0,0,W,0,0,1]);
if P mod 3 = 1 then
  Add(params[16],[0,0,W^2 mod P,0,0,1]);
  Add(params[16],[0,0,W^3 mod P,0,0,1]);
  Add(params[16],[0,0,W^4 mod P,0,0,1]);
  Add(params[16],[0,0,W^5 mod P,0,0,1]);
fi;
for z in [0..P-1] do
  Add(params[16],[0,0,z,1,0,1]);
  Add(params[16],[0,0,z,W,0,1]);
od;
for z in [0..P-1] do
for t in [0..P-1] do
  Add(params[16],[0,1,z,t,0,1]);
  Add(params[16],[0,W,z,t,0,1]);
od;
od;
for y in [0..P-1] do
for z in [0..P-1] do
for t in [0..P-1] do
  Add(params[16],[1,y,z,t,0,1]);
  Add(params[16],[W,y,z,t,0,1]);
od;
od;
od;

Add(params[16],[0,0,0,0,1,0]);
Add(params[16],[0,1,0,0,1,0]);
Add(params[16],[0,W,0,0,1,0]);
if P mod 3 = 1 then
  Add(params[16],[0,W^2 mod P,0,0,1,0]);
  Add(params[16],[0,W^3 mod P,0,0,1,0]);
  Add(params[16],[0,W^4 mod P,0,0,1,0]);
  Add(params[16],[0,W^5 mod P,0,0,1,0]);
fi;
if P mod 4 = 1 then
  for v in [1,W,W^2 mod P,W^3 mod P] do
  for y in [0..(P-1)/2] do  
    Add(params[16],[0,y,0,0,1,v]);
  od;
  od;
fi;
if P mod 4 = 3 then
  for v in [1,W] do
  for y in [0..P-1] do  
    Add(params[16],[0,y,0,0,1,v]);
  od;
  od;
fi;
for v in [0..P-1] do
for y in [0..P-1] do
  Add(params[16],[0,y,0,1,1,v]);
  Add(params[16],[0,y,0,W,1,v]);
od;
od;
for t in [0..P-1] do
for v in [0..P-1] do
for y in [0..P-1] do
  Add(params[16],[0,y,1,t,1,v]);
  Add(params[16],[0,y,W,t,1,v]);
od;
od;
od;
for z in [0..P-1] do
for t in [0..P-1] do
for v in [0..P-1] do
for y in [0..P-1] do
  Add(params[16],[1,y,z,t,1,v]);
  Add(params[16],[W,y,z,t,1,v]);
od;
od;
od;
od;

if Length(params[16]) <> 2*P^4+4*P^3+8*P^2+14*P+11+4*Gcd(P-1,3)+3*Gcd(P-1,4)
then Error("params[16] wrong"); fi;
if case = 16 then return params[16]; fi;
fi;

##Case 17, cb=bac=0, cac=baa, caa=bab
if case = true or case = 17 then

mats:=[];

for y5 in [0,1,lns] do
for y6 in [0..P-1] do
for y3 in [0..P-1] do
for y4 in [0..P-1] do

A:=One(F)*MyCutVector([y3,y4,y5,y6], 2);
if A = 0*A then
  Add(mats,A);
  continue;
fi;

new:=1;
index:=P^3*y5+P^2*y6+P*y3+y4;

for a in [1..P-1] do
for x in CU do
c:=a*x;

B:=One(F)*MyCutVector([a^-1*c^2,0,0,c], 2);
C:=One(F)*MyCutVector([a*c^2,0,0,a^2*c], 2);

D:=B*A*C^-1;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);


ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=0; fi;

B:=One(F)*MyCutVector([0,a^-1*c^2,c,0], 2);
C:=One(F)*MyCutVector([0,a*c^2,a^2*c,0], 2);

D:=B*A*C^-1;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);


ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=0; fi;


if new = 0 then break; fi;
od;
if new = 0 then break; fi;
od;

if new = 1 then
  Add(mats,A);
fi;

od;
od;
od;
od;


for AS in mats do
stab1:=[];
stab2:=[];
for a in [1..P-1] do
for x in CU do
c:=a*x;

B:=One(F)*MyCutVector([a^-1*c^2,0,0,c], 2);
C:=One(F)*MyCutVector([a*c^2,0,0,a^2*c], 2);

D:=B*AS*C^-1;

if D = AS then
  Add(stab1,One(F)*MyCutVector([a,0,0,0,a^-1*c^2,0,0,0,c], 3));
  Add(stab2,C);
fi;

B:=One(F)*MyCutVector([0,a^-1*c^2,c,0], 2);
C:=One(F)*MyCutVector([0,a*c^2,a^2*c,0], 2);

D:=B*AS*C^-1;

if D = AS then
  Add(stab1,One(F)*MyCutVector([a,0,0,0,0,a^-1*c^2,0,c,0], 3));
  Add(stab2,C);
fi;

od;
od;

y3:=IntFFE(AS[1][1]); 
y4:=IntFFE(AS[1][2]); 
y5:=IntFFE(AS[2][1]); 
y6:=IntFFE(AS[2][2]);

for y1 in [0..P-1] do
for y2 in [0..P-1] do

A:=One(F)*MyCutVector([y1,y2,y3,y4,y5,y6], 3);
if A = 0*A then
  Add(params[17],[0,0,0,0,0,0]);
  continue;
fi;

new:=1;
index:=P*y1+y2;

for ii in [1..Length(stab1)] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);

ind1:=P*z1+z2;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;


if new = 1 then
  Add(params[17],[y1,y2,y3,y4,y5,y6]);
fi;

od;
od;

od;

if Length(params[17]) <> 
(P^4+2*P^3+3*P^2+4*P+2)*(P-1)/Gcd(P-1,3)+3*P+4+(P^2+P+1)*Gcd(P-1,4)/2
then Error("params[17] wrong"); fi;
if case = 17 then return params[17]; fi;
fi;

##Case 18, cb=bac=0, cac=baa, caa=wbab (p=1mod3 only)
if case = true or case = 18 then 

if P mod 3 = 1 then

mats:=[];

for y5 in [0,1,lns] do
for y6 in [0..P-1] do
for y3 in [0..P-1] do
for y4 in [0..P-1] do

A:=One(F)*MyCutVector([y3,y4,y5,y6], 2);
if A = 0*A then
  Add(mats,A);
  continue;
fi;

new:=1;
index:=P^3*y5+P^2*y6+P*y3+y4;

for a in [1..P-1] do
for x in CU do
c:=a*x;

B:=One(F)*MyCutVector([a^-1*c^2,0,0,c], 2);
C:=One(F)*MyCutVector([a*c^2,0,0,a^2*c], 2);

D:=B*A*C^-1;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);


ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=0; fi;


if new = 0 then break; fi;
od;
if new = 0 then break; fi;
od;

if new = 1 then
  Add(mats,A);
fi;

od;
od;
od;
od;

for AS in mats do
stab1:=[];
stab2:=[];
for a in [1..P-1] do
for x in CU do
c:=a*x;

B:=One(F)*MyCutVector([a^-1*c^2,0,0,c], 2);
C:=One(F)*MyCutVector([a*c^2,0,0,a^2*c], 2);

D:=B*AS*C^-1;

if D = AS then
  Add(stab1,One(F)*MyCutVector([a,0,0,0,a^-1*c^2,0,0,0,c], 3));
  Add(stab2,C);
fi;

od;
od;

y3:=IntFFE(AS[1][1]); 
y4:=IntFFE(AS[1][2]); 
y5:=IntFFE(AS[2][1]); 
y6:=IntFFE(AS[2][2]);

for y1 in [0..P-1] do
for y2 in [0..P-1] do

A:=One(F)*MyCutVector([y1,y2,y3,y4,y5,y6], 3);
if A = 0*A then
  Add(params[18],[0,0,0,0,0,0]);
  continue;
fi;

new:=1;
index:=P*y1+y2;

for ii in [1..Length(stab1)] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);

ind1:=P*z1+z2;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;


if new = 1 then
  Add(params[18],[y1,y2,y3,y4,y5,y6]);
fi;

od;
od;

od;

if Length(params[18]) <> (2*P^5+2*P^4+2*P^3+2*P^2+14*P+17)/3
then Error("params[18] wrong"); fi;
fi;
if case = 18 then return params[18]; fi;
fi; 

##Case 19, cb=baa=caa=cac=0
if case = true or case = 19 then

Add(params[19],[0,0,0,0,0,0]);
Add(params[19],[1,0,0,0,0,0]);
Add(params[19],[W,0,0,0,0,0]);
Add(params[19],[0,1,0,0,0,0]);

for x in [0,1,W] do
for y in [0,1] do
Add(params[19],[x,y,1,0,0,0]);
od;
od;

Add(params[19],[0,0,0,1,0,0]);
Add(params[19],[1,0,0,1,0,0]);
Add(params[19],[W,0,0,1,0,0]);
Add(params[19],[0,1,0,1,0,0]);

for t in [0..P-1] do
  Add(params[19],[0,0,1,t,1,0]);
  Add(params[19],[0,1,1,t,1,0]);
  Add(params[19],[0,W,1,t,1,0]);
  for y in [0..P-1] do
    Add(params[19],[1,y,1,t,1,0]);
    Add(params[19],[W,y,1,t,1,0]);
  od;
od;

for x in [0,1,W] do
for t in [0,1,W] do
  Add(params[19],[x,0,0,t,1,0]);
od;
od;

for x in [1,W] do
for t in [0..P-1] do
  Add(params[19],[x,1,0,t,1,0]);
od;
od;

Add(params[19],[0,1,0,0,1,0]);
Add(params[19],[0,1,0,1,1,0]);
Add(params[19],[0,1,0,W,1,0]);
if P mod 4 = 1 then
  Add(params[19],[0,1,0,W^2 mod P,1,0]);
  Add(params[19],[0,1,0,W^3 mod P,1,0]);
fi;

half:=(P+1)/2;
for x in [0..P-1] do
  if x = half then continue; fi;
  for z in [0,1,W] do
  for t in [0,1] do
    Add(params[19],[z,t,x,0,0,1]);
  od;
  od;
od;

Add(params[19],[1,0,half,0,0,1]);
Add(params[19],[W,0,half,0,0,1]);
Add(params[19],[1,0,half,1,0,1]);
Add(params[19],[W,0,half,1,0,1]);
Add(params[19],[0,0,half,0,0,1]);
Add(params[19],[0,0,half,1,0,1]);
Add(params[19],[0,1,half,0,0,1]);
Add(params[19],[0,1,half,1,0,1]);
Add(params[19],[0,W,half,1,0,1]);

if Length(params[19]) <> 2*P^2+11*P+27+Gcd(P-1,4)
then Error("params[19] wrong"); fi;
if case = 19 then return params[19]; fi;
fi;

##Case 20, cb=baa=cac, caa=bab
if case = true or case = 20 then 

for u in [0,1,W] do
  Add(params[20],[0,0,0,0,u,0]);
  Add(params[20],[0,0,0,1,u,0]);
  Add(params[20],[0,0,0,W,u,0]);
  for t in [0..P-1] do
    Add(params[20],[1,0,0,t,u,0]);
    Add(params[20],[W,0,0,t,u,0]);
  od;
od;

Add(params[20],[0,1,0,0,0,0]);
Add(params[20],[0,1,0,0,1,0]);
Add(params[20],[0,1,0,0,W,0]);
if P mod 3 = 1 then
  Add(params[20],[0,1,0,0,W^2 mod P,0]);
  Add(params[20],[0,1,0,0,W^3 mod P,0]);
  Add(params[20],[0,1,0,0,W^4 mod P,0]);
  Add(params[20],[0,1,0,0,W^5 mod P,0]);
fi;

for u in [0..P-1] do
  Add(params[20],[0,1,0,1,u,0]);
  Add(params[20],[0,1,0,W,u,0]);
od;

for t in [0..P-1] do
for u in [0..P-1] do
  Add(params[20],[1,1,0,t,u,0]);
  Add(params[20],[W,1,0,t,u,0]);
od;
od;

Add(params[20],[0,0,0,0,0,1]);
Add(params[20],[0,1,0,0,0,1]);
Add(params[20],[0,W,0,0,0,1]);
if P mod 4 = 1 then
  Add(params[20],[0,W^2 mod P,0,0,0,1]);
  Add(params[20],[0,W^3 mod P,0,0,0,1]);
fi;

for y in [0..P-1] do
  Add(params[20],[0,y,0,0,1,1]);
  Add(params[20],[0,y,0,0,W,1]);
od;

for y in [0..P-1] do
for u in [0..P-1] do
  Add(params[20],[0,y,0,1,u,1]);
  Add(params[20],[0,y,0,W,u,1]);
od;
od;

for y in [0..P-1] do
for t in [0..P-1] do
for u in [0..P-1] do
  Add(params[20],[1,y,0,t,u,1]);
  Add(params[20],[W,y,0,t,u,1]);
od;
od;
od;

for v in [0..P-1] do

Add(params[20],[0,0,1,0,0,v]);
Add(params[20],[0,1,1,0,0,v]);
Add(params[20],[0,W,1,0,0,v]);
if P mod 4 = 1 then
  Add(params[20],[0,W^2 mod P,1,0,0,v]);
  Add(params[20],[0,W^3 mod P,1,0,0,v]);
fi;

for y in [0..P-1] do
  Add(params[20],[0,y,1,0,1,v]);
  Add(params[20],[0,y,1,0,W,v]);
od;

for y in [0..P-1] do
for u in [0..P-1] do
  Add(params[20],[0,y,1,1,u,v]);
  Add(params[20],[0,y,1,W,u,v]);
od;
od;

for y in [0..P-1] do
for t in [0..P-1] do
for u in [0..P-1] do
  Add(params[20],[1,y,1,t,u,v]);
  Add(params[20],[W,y,1,t,u,v]);
od;
od;
od;

od;

if Length(params[20]) <> 2*P^4+4*P^3+6*P^2+11*P+11+2*Gcd(P-1,3)+(P+1)*Gcd(P-1,4)
then Error("params[20] wrong"); fi;
if case = 20 then return params[20]; fi;
fi;


##Case 21, cb=caa=cac=0, bab=baa
if case = true or case = 21 then 

Add(params[21],[0,0,0,0,0,0]);
Add(params[21],[1,0,0,0,0,0]);
Add(params[21],[W,0,0,0,0,0]);
Add(params[21],[0,1,0,0,0,0]);

for x in [0..P-1] do
  Add(params[21],[x,0,1,0,0,0]);
  Add(params[21],[x,0,W,0,0,0]);
  Add(params[21],[x,1,1,0,0,0]);
  Add(params[21],[x,1,W,0,0,0]);
od;

Add(params[21],[1,0,0,1,0,0]);
Add(params[21],[W,0,0,1,0,0]);
for x in [0..P-1] do
  Add(params[21],[0,x,0,1,0,0]);
od;

for y in [0..P-1] do
for z in [0..P-1] do
for t in [0..P-1] do
  Add(params[21],[1,y,z,t,1,0]);
  Add(params[21],[W,y,z,t,1,0]);
od;
od;
od;

for y in [0..P-1] do
for t in [0..P-1] do
  Add(params[21],[0,y,1,t,1,0]);
  Add(params[21],[0,y,W,t,1,0]);
od;
od;

for t in [0..P-1] do
  Add(params[21],[0,1,0,t,1,0]);
  Add(params[21],[0,W,0,t,1,0]);
  if P mod 4 = 1 then
    Add(params[21],[0,W^2 mod P,0,t,1,0]);
    Add(params[21],[0,W^3 mod P,0,t,1,0]);
  fi;
od;

Add(params[21],[0,0,0,0,1,0]);
Add(params[21],[0,0,0,1,1,0]);
Add(params[21],[0,0,0,W,1,0]);
if P mod 4 = 1 then
  Add(params[21],[0,0,0,W^2 mod P,1,0]);
  Add(params[21],[0,0,0,W^3 mod P,1,0]);
fi;

for x in [0..P-1] do
if x = 1 then continue; fi;
for z in [0..P-1] do
for t in [0,1] do
Add(params[21],[x,0,z,t,0,1]);
od;
od;
od;

for x in [0..P-1] do
if x = W then continue; fi;
for z in [0..P-1] do
for t in [0,1] do
Add(params[21],[x,0,z,t,0,W]);
od;
od;
od;

half:=(P+1)/2;
wover2:=One(F)*(W/2); wover2:=IntFFE(wover2);
for z in [0..P-1] do
  if z = half then continue; fi;
  Add(params[21],[1,0,z,0,0,1]);
  Add(params[21],[1,1,z,0,0,1]);
od;
for z in [0..P-1] do
  if z = wover2 then continue; fi;
  Add(params[21],[W,0,z,0,0,W]);
  Add(params[21],[W,1,z,0,0,W]);
od;

Add(params[21],[1,0,half,0,0,1]);
Add(params[21],[W,0,wover2,0,0,W]);
Add(params[21],[1,0,half,1,0,1]);
Add(params[21],[W,0,wover2,1,0,W]);
for t in [0..P-1] do
  Add(params[21],[1,1,half,t,0,1]);
  Add(params[21],[W,1,wover2,t,0,W]);
od;

if Length(params[21]) <> 2*P^3+6*P^2+7*P+7+(P+1)*Gcd(P-1,4)
then Error("params[21] wrong"); fi;
if case = 21 then return params[21]; fi;
fi;

##Case 22, cb=baa=caa=0, cac=wbab
if case = true or case = 22 then 

mats:=[];

bcrange:=[[0,1]];
for i in [0..P-1] do
  Add(bcrange,[1,i]);
od;

for y3 in [0,1] do
for y4 in [0..P-1] do
for y5 in [0..P-1] do
for y6 in [0..P-1] do

A:=One(F)*MyCutVector([y3,y4,y5,y6], 2);
if A = 0*A then
  Add(mats,A);
  continue;
fi;

new:=1;
index:=P^3*y3+P^2*y4+P*y5+y6;

for bc in bcrange do
b:=bc[1]; c:=bc[2];
for m in [-1,1] do

B:=One(F)*MyCutVector([W*b,m*c,W*c,W*m*b], 2);
C:=One(F)*MyCutVector([W*(W*b^2+c^2),m*2*W*b*c,2*W^2*b*c,m*W*(W*b^2+c^2)], 2);

D:=B*A*C^-1;
if D[1][1] <> Zero(F) then D:=D[1][1]^-1*D; fi;
if D[1][1] = Zero(F) and D[1][2] <> Zero(F) 
then D:=D[1][2]^-1*D; fi;
if D[1][1] = Zero(F) and D[1][2] = Zero(F) and D[2][1] <> Zero(F) 
then D:=D[2][1]^-1*D; fi;
if D[1][1]=Zero(F) and D[1][2]=Zero(F) and D[2][1]=Zero(F) and D[2][2]<>Zero(F) 
then D:=D[2][2]^-1*D; fi;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);


ind1:=P^3*z3+P^2*z4+P*z5+z6;

if ind1 < index then new:=0; fi;


if new = 0 then break; fi;
od;
if new = 0 then break; fi;
od;

if new = 1 then Add(mats,A); fi;

od;
od;
od;
od;

for AS in mats do
stab1:=[];
stab2:=[];
for b in [0..P-1] do
for c in [0..P-1] do
for m in [-1,1] do
if b+c = 0 then continue; fi;

B:=One(F)*MyCutVector([W*b,m*c,W*c,W*m*b], 2);
C:=One(F)*MyCutVector([W*(W*b^2+c^2),m*2*W*b*c,2*W^2*b*c,m*W*(W*b^2+c^2)], 2);

D:=B*AS*C^-1;
if D[1][1] <> Zero(F) 
then D:=D[1][1]^-1*D; fi;
if D[1][1] = Zero(F) and D[1][2] <> Zero(F) 
then D:=D[1][2]^-1*D; fi;
if D[1][1] = Zero(F) and D[1][2] = Zero(F) and D[2][1] <> Zero(F) 
then D:=D[2][1]^-1*D; fi;
if D[1][1]=Zero(F) and D[1][2]=Zero(F) and D[2][1]=Zero(F) and D[2][2]<>Zero(F) 
then D:=D[2][2]^-1*D; fi;

if D = AS then
  Add(stab1,One(F)*MyCutVector([1,0,0,0,W*b,m*c,0,W*c,W*m*b], 3));
  Add(stab2,C);
fi;

od;
od;
od;

y3:=IntFFE(AS[1][1]); 
y4:=IntFFE(AS[1][2]); 
y5:=IntFFE(AS[2][1]); 
y6:=IntFFE(AS[2][2]);

for y1 in [0..P-1] do
for y2 in [0..P-1] do

A:=One(F)*MyCutVector([y1,y2,y3,y4,y5,y6], 3);
if A = 0*A then
  Add(params[22],[0,0,0,0,0,0]);
  continue;
fi;

new:=1;
index:=P*y1+y2;

for ii in [1..Length(stab1)] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);

ind1:=P*z1+z2;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;


if new = 1 then
  Add(params[22],[y1,y2,y3,y4,y5,y6]);
fi;

od;
od;

od;

if Length(params[22]) <> (2*P^3+3*P^2+3*P+13-Gcd(P-1,3)+(P+1)*Gcd(P-1,4))/2
then Error("params[22] wrong"); fi;
if case = 22 then return params[22]; fi;
fi;

## Case 23, cb=baa=0, caa=bac, cac=wbab
if case = true or case = 23 then 

sol:=0;
if P mod 3 = 2 then
  for x in [1..P-1] do
    if One(F)*(12*W*x^2+1) = Zero(F) then sol:=x; break; fi;
  od;
fi;

mats:=[];

for y5 in [0,1,lns] do
for y6 in [0..P-1] do
for y3 in [0..P-1] do
for y4 in [0..P-1] do

A:=One(F)*MyCutVector([y3,y4,y5,y6], 2);
if A = 0*A then
  Add(mats,A);
  continue;
fi;

new:=1;
index:=P^3*y5+P^2*y6+P*y3+y4;

for a in [1..P-1] do
for m in [-1,1] do

B:=One(F)*MyCutVector([a,0,0,m*a], 2);
C:=One(F)*MyCutVector([a^3,0,0,m*a^3], 2);

D:=B*A*C^-1;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);


ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=0; fi;

if P mod 3 = 2 and new = 1 then
for b in [sol,-sol] do

B:=One(F)*MyCutVector([-2*W*a*b,a,m*W*a,-2*m*W*a*b], 2);
C:=One(F)*MyCutVector([8/3*W^2*a^3*b,4/3*W*a^3,4/3*m*W^2*a^3,8/3*m*W^2*a^3*b], 2);

D:=B*A*C^-1;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);

ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;

fi;

if new = 0 then break; fi;
od;
if new = 0 then break; fi;
od;

if new = 1 then
  Add(mats,A);
fi;

od;
od;
od;
od;

for AS in mats do
stab1:=[];
stab2:=[];
for a in [1..P-1] do
for m in [-1,1] do

B:=One(F)*MyCutVector([a,0,0,m*a], 2);
C:=One(F)*MyCutVector([a^3,0,0,m*a^3], 2);

D:=B*AS*C^-1;

if D = AS then
  Add(stab1,One(F)*MyCutVector([a,0,0,0,a,0,0,0,m*a], 3));
  Add(stab2,C);
fi;

if P mod 3 = 2 then
for b in [sol,-sol] do

B:=One(F)*MyCutVector([-2*W*a*b,a,m*W*a,-2*m*W*a*b], 2);
C:=One(F)*MyCutVector([8/3*W^2*a^3*b,4/3*W*a^3,4/3*m*W^2*a^3,8/3*m*W^2*a^3*b], 2);

D:=B*AS*C^-1;

if D = AS then
  Add(stab1,
  One(F)*MyCutVector([4*W*a*b,-3*W*a*b,a/2,0,-2*W*a*b,a,0,m*W*a,-2*m*W*a*b], 3));
  Add(stab2,C);
fi;

od;
fi;

od;
od;

y3:=IntFFE(AS[1][1]); 
y4:=IntFFE(AS[1][2]); 
y5:=IntFFE(AS[2][1]); 
y6:=IntFFE(AS[2][2]);

for y1 in [0..P-1] do
for y2 in [0..P-1] do

A:=One(F) * MyCutVector([y1,y2,y3,y4,y5,y6], 3);
if A = 0*A then
  Add(params[23],[0,0,0,0,0,0]);
  continue;
fi;

new:=1;
index:=P*y1+y2;

for ii in [1..Length(stab1)] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);

ind1:=P*z1+z2;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;


if new = 1 then
  Add(params[23],[y1,y2,y3,y4,y5,y6]);
fi;

od;
od;

od;

if P mod 3 = 1 then
  expect:=P^5+P^4+P^3+P^2+P+2+(P^2+P+1)*Gcd(P-1,4)/2;
else
  expect:=P^5/3+P^4/3+P^3/3+P^2/3+P+2+(P^2+P+1)*Gcd(P-1,4)/2;
fi;

if Length(params[23]) <> expect then Error("params[23] wrong"); fi;
if case = 23 then return params[23]; fi;
fi;

## Case 24, cb=baa=0, caa=xbab+bac, cac=wbab
## where x is not a value of y(y^2+3w)/(3y^2+w)
if case = true or case = 24 then 

if P mod 3 = 2 then

val:=[];
for x in [0..P-1] do
  if One(F)*(3*x^2+W) <> Zero(F) then
    a:=One(F)*(x*(x^2+3*W))*(One(F)*(3*x^2+W))^-1;
    AddSet(val,IntFFE(a));
  fi;
od;
for a in [1..P-1] do
  if not a in val then s:=a; break; fi;
od;

for y in [1..P-1] do
  if One(F)*(W*y^2+3) = Zero(F) then b:=y; break; fi;
od;

B1:=One(F)*MyCutVector([2,2*b,2*W*b,2], 2);
C1:=One(F)*MyCutVector([32,-32*b,-32*W*b,32], 2);
B2:=One(F)*MyCutVector([2,-2*b,-2*W*b,2], 2);
C2:=One(F)*MyCutVector([32,32*b,32*W*b,32], 2);
BB1:=One(F)*MyCutVector([-4,s*b+3,3*s*W^-1+b,0,2,2*b,0,2*W*b,2], 3);
BB2:=One(F)*MyCutVector([-4,-s*b+3,3*s*W^-1-b,0,2,-2*b,0,-2*W*b,2], 3);

mats:=[];

for y5 in [0,1,lns] do
for y6 in [0..P-1] do
for y3 in [0..P-1] do
for y4 in [0..P-1] do

A:=One(F)*MyCutVector([y3,y4,y5,y6], 2);
if A = 0*A then
  Add(mats,A);
  continue;
fi;

new:=true;
index:=P^3*y5+P^2*y6+P*y3+y4;

A1:=B1*A*C1^-1;
A2:=B2*A*C2^-1;

for a in [1..(P-1)/2] do

D:=(One(F)*a)^-2*A;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);


ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=false; break; fi;

D:=(One(F)*a)^-2*A1;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);

ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=false; break; fi;

D:=(One(F)*a)^-2*A2;

z3:=IntFFE(D[1][1]);
z4:=IntFFE(D[1][2]);
z5:=IntFFE(D[2][1]);
z6:=IntFFE(D[2][2]);

ind1:=P^3*z5+P^2*z6+P*z3+z4;

if ind1 < index then new:=false; break; fi;

od;

if new then
  Add(mats,A);
fi;

od;
od;
od;
od;


for AS in mats do
stab1:=[];
stab2:=[];

AS1:=B1*AS*C1^-1;
AS2:=B2*AS*C2^-1;

for a in [1..(P-1)/2] do

D:=(One(F)*a)^-2*AS;

if D = AS then
  Add(stab1,One(F)*MyCutVector([1,0,0,0,1,0,0,0,1], 3));
  Add(stab2,One(F)*MyCutVector([a^2,0,0,a^2], 2));
fi;

D:=(One(F)*a)^-2*AS1;

if D = AS then
  Add(stab1,BB1);
  Add(stab2,a^2*C1);
fi;

D:=(One(F)*a)^-2*AS2;

if D = AS then
  Add(stab1,BB2);
  Add(stab2,a^2*C2);
fi;

od;

y3:=IntFFE(AS[1][1]); 
y4:=IntFFE(AS[1][2]); 
y5:=IntFFE(AS[2][1]); 
y6:=IntFFE(AS[2][2]);

for y1 in [0..P-1] do
for y2 in [0..P-1] do

A:=One(F)*MyCutVector([y1,y2,y3,y4,y5,y6], 3); 
if A = 0*A then
  Add(params[24],[0,0,0,0,0,0,s]);
  continue;
fi;

new:=1;
index:=P*y1+y2;

for ii in [1..Length(stab1)] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);

ind1:=P*z1+z2;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;


if new = 1 then
  Add(params[24],[y1,y2,y3,y4,y5,y6,s]);
fi;

od;
od;

od;

if Length(params[24]) <> 2*(P^5+P^4+P^3+P^2)/3+2*P+3 then 
Error("params[24] wrong"); fi;
fi;
if case = 24 then return params[24]; fi;
fi; 

if Sum(List(params, Length))<> 2*P^5+7*P^4+19*P^3+49*P^2+128*P+256
  +(P^2+7*P+29)*Gcd(P-1,3)+ (P^2+7*P+24)*Gcd(P-1,4)+(P+3)*Gcd(P-1,5) then 
    Error("sum is wrong"); 
fi;

for i in [1..Length(params)] do
    if not Length(params[i]) = Length(Set(params[i])) then 
        Error("params ",i," has duplicates");
    fi;
    if not IsSubset( [0..P-1], Set(Flat(params[i]))) then 
        Error("params ",i," is out of range");
    fi;
od;

return params;
end );

BindGlobal( "ValsFunction28", function(P, case)
    local vals, i;
    vals := ValsPreFunction28(P, case);
    if IsInt(case) and case < 24 then 
       return List(vals, x -> [x[4], x[5], x[6], x[1], x[2], x[3]]);
    elif IsInt(case) and case = 24 then 
       return List(vals, x -> [x[4], x[5], x[6], x[1], x[7], x[2], x[3]]);
    else
       for i in [1..23] do
           vals[i] := List(vals[i], x->[x[4],x[5],x[6],x[1],x[2],x[3]]);
       od;
       vals[24] := List(vals[24], x->[x[4],x[5],x[6],x[1],x[7],x[2],x[3]]);
       return vals;
    fi;
end );

