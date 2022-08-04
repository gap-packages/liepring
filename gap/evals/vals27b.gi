
BindGlobal( "Leading27b", function(A)
  if A[1][1] <> 0*A[1][1] then return A[1][1]; fi;
  if A[1][2] <> 0*A[1][1] then return A[1][2]; fi;
  if A[2][1] <> 0*A[1][1] then return A[2][1]; fi;
  return A[2][2];
end );

BindGlobal( "Linearise27b", function(list)
   return list{[10,11,12,7,1,2,3,4,5,6,8,9]};
end );

BindGlobal( "ValsFunction27b", function(P)
local W, F, SQ, x, y, i, range, y1, y2, y3, y4, y4range, y3range, new, index,
      A, r, a, d, b, c, B, C, k, D, z1, z2, z3, z4, ind1, l, 
      y5, y6, A2, D2, n, D3, params, mats, k2, k1, s, c1, b1, u, lns, v,
      lastt, t, z; 

W:=PrimitiveRootMod(P);
F:=GF(P);
params:=[];

SQ:=[];
for x in [1..((P-1)/2)] do
  y:=x^2 mod P;
  Add(SQ,y);
od;
for i in [2..P-1] do
  if not i in SQ then lns:=i; break; fi;
od;

range:=[[0,1]];
for i in [0..P-1] do
  Add(range,[1,i]);
od;

Add(params,Linearise27b([0,0,0,0,0,0,0,0,1,0,0,0]));
Add(params,Linearise27b([0,0,0,0,1,0,0,0,1,0,0,0]));

y1:=0;
for y2 in [0,1] do
y3range:=[0..P-1];
if y2 = 0 then y3range:=[0,1]; fi;
for y3 in y3range do
y4range:=[0..P-1];
if y1 = 0 and y2 = 0 then y4range:=[0,1]; fi;
for y4 in y4range do

new:=1;
index:=P^3*y1+P^2*y2+P*y3+y4;

A:=MyCutVector([y1,y2,y3,y4], 2)*One(F);

if A = 0*A then
  Add(params,[0,0,0,0,0,0,0,0,0,0,0,0]);
  continue;
fi;

for r in range do
a:=r[1];
d:=r[2];
for s in range do
b:=s[1];
c:=s[2];

B:=MyCutVector([c,W*b,b,c], 2)*One(F);
C:=MyCutVector([a*c-W*b*d,W*a*b-W*c*d, a*b-c*d,a*c-W*b*d], 2)*One(F);

D:=B*A*C^-1;

k:=Leading27b(D);
D:=k^-1*D;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);
z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P^3*z1+P^2*z2+P*z3+z4;

if ind1 < index then new:=0; break; fi;

B[1]:=-B[1];
C[1]:=-C[1];

D:=B*A*C^-1;

k:=Leading27b(D);
D:=k^-1*D;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);
z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P^3*z1+P^2*z2+P*z3+z4;

if ind1 < index then new:=0; break; fi;

od;
if new = 0 then break; fi;
od;

if new = 1 then
  Add(params,Linearise27b([0,0,0,0,y2,0,y3,y4,0,0,0,0]));
fi;

od;
od;
od;

for z in [1..((P-1)/2)] do
for u in [0,1] do
for t in [0..P-1] do
  Add(params,Linearise27b([0,0,0,0,0,z,u,0,t,0,1,0]));
od;
od;
od;

for t in [1..P-1] do
for x in [0,1] do
  Add(params,Linearise27b([0,0,0,x,0,0,0,0,t,0,1,0]));
od;
od;

Add(params,Linearise27b([0,0,0,0,0,0,0,0,0,0,1,0]));

Add(params,Linearise27b([0,0,0,0,0,0,1,0,0,0,1,0]));

for u in [0..((P-1)/2)] do
  Add(params,Linearise27b([0,0,0,1,0,0,u,0,0,0,1,0]));
od;

for y1 in [0,1,lns] do
for y2 in [0..(P-1)/2] do
for y3 in [0..P-1] do
for y4 in [0..P-1] do

A:=MyCutVector([y1,y2,y3,y4], 2)*One(F);

if A = 0*A then
  Add(params,Linearise27b([0,0,0,0,0,0,0,0,0,0,0,1]));
  continue;
fi;

new:=1;
index:=P^3*y1+P^2*y2+P*y3+y4;

for r in range do
b:=r[1]; c:=r[2];

B:=MyCutVector([c,W*b,b,c], 2)*One(F);
C:=MyCutVector([c*(c^2-W*b^2),W*b*(c^2-W*b^2), 
              b*(c^2-W*b^2),c*(c^2-W*b^2)], 2)*One(F);

D:=B*A*C^-1;

k:=Leading27b(D);
D:=k^-1*D;
if IsSquareGF(F, k)=false then D:=lns*D; fi;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);
z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P^3*z1+P^2*z2+P*z3+z4;

if ind1 < index then new:=0; break; fi;

B:=MyCutVector([-c,-W*b,b,c], 2)*One(F);
C:=MyCutVector([-c*(c^2-W*b^2),-W*b*(c^2-W*b^2),
               b*(c^2-W*b^2),c*(c^2-W*b^2)], 2)*One(F);

D:=B*A*C^-1;

k:=Leading27b(D);
D:=k^-1*D;
if IsSquareGF(F,k)=false then D:=lns*D; fi;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);
z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P^3*z1+P^2*z2+P*z3+z4;

if ind1 < index then new:=0; break; fi;

od;

if new = 1 then
  Add(params,Linearise27b([0,0,0,y1,y2,0,y3,y4,0,0,0,1]));
fi;

od;
od;
od;
od;

for x in [0..P-1] do
for u in [0..((P-1)/2)] do
for v in [0..P-1] do
lastt:=P-1;
if u = 0 then lastt:=(P-1)/2; fi;
for t in [0..lastt] do
  Add(params,Linearise27b([0,1,0,x,0,0,u,v,t,0,0,1]));
od;
od;
od;
od;

y1:=0;
y2:=1;
for y3 in [1..P-1] do
for y4 in [0..(P-1)/2] do

new:=1;
index:=P*y3+y4;

A:=MyCutVector([y1,y2,y3,y4], 2)*One(F);

for r in range do
a:=r[1];
d:=r[2];
for s in range do
b1:=s[1];
c1:=s[2];
if (b1*W*y3*d^2+b1*a*y4*d+b1*a^2-y4*c1*d^2-d*a*c1-d*y3*a*c1) mod P <> 0 then 
  continue; 
fi;
k1:=(-b1*W*d*(y4*d+a+a*y3)+c1*(W*y3*d^2+y4*a*d+a^2)) * One(F);
k2:=((a^2-W*d^2)*(c1^2-W*b1^2)) * One(F);
k:=IntFFE(k1*k2^-1);
b:=k*b1;
c:=k*c1;

B:=MyCutVector([a,d,W*d,a], 2)*One(F);
C:=MyCutVector([a*c-W*b*d,W*a*b-W*c*d,
        a*b-c*d,a*c-W*b*d], 2)*One(F);

D:=B*A*C^-1;

z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P*z3+z4;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;
if new = 0 then break; fi;
od;

if new = 1 then
  Add(params,Linearise27b([0,1,0,0,0,0,0,0,0,y3,y4,0]));
fi;

od;
od;

mats:=[];

for y1 in [0..P-1] do
for y2 in [0..((P-1)/2)] do
for y3 in [0..P-1] do
for y4 in [0..P-1] do

new:=1;
index:=P^3*y1+P^2*y2+P*y3+y4;

A:=MyCutVector([y1,y2,y3,y4], 2)*One(F);
if Rank(A) <> 2 then continue; fi;

for r in range do
a:=r[1];
d:=r[2];

B:=MyCutVector([a,d,W*d,a], 2)*One(F);
C:=MyCutVector([a,-W*d,-d,a], 2)*One(F);

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);
z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P^3*z1+P^2*z2+P*z3+z4;

if ind1 < index then new:=0; fi;

B[2]:=-B[2];
C[2]:=-C[2];

D:=B*A*C^-1;

z1:=IntFFE(D[1][1]);
z2:=IntFFE(D[1][2]);
z3:=IntFFE(D[2][1]);
z4:=IntFFE(D[2][2]);

ind1:=P^3*z1+P^2*z2+P*z3+z4;

if ind1 < index then new:=0; fi;

if new = 0 then break; fi;
od;

if new = 1 then
  Add(params,Linearise27b([y1,y2,0,0,0,0,0,0,1,y3,y4,0]));
  Add(mats,A);
fi;


od;
od;
od;
od;

for A in mats do
y1:=IntFFE(A[1][1]);
y2:=IntFFE(A[1][2]);
y5:=IntFFE(A[2][1]);
y6:=IntFFE(A[2][2]);

y3:=0; y4:=1;

A2:=MyCutVector([y1,y2,0,0,0,0,y3,y4,1,y5,y6,0], 4)*One(F);
new:=1;

for r in range do
a:=r[1];
d:=r[2];

B:=MyCutVector([a,d,W*d,a], 2)*One(F);
C:=MyCutVector([a,-W*d,-d,a], 2)*One(F);
D:=B*A*C^-1;

B[2]:=-B[2];
C[2]:=-C[2];
D2:=B*A*C^-1;

if D = A or D2 = A then

for n in [0..P-1] do
for x in [0..P-1] do

if D = A then

B:=MyCutVector([a,0,0,d, 0,1,0,0, n,0,1,x, W*d,0,0,a], 4)*One(F);
C:=MyCutVector([a,-W*d,0, -d,a,0, -n,W*x,1], 3)*One(F);
D3:=B*A2*C^-1;

z1:=IntFFE(D3[3][1]);
z2:=IntFFE(D3[3][2]);

if z1+z2 = 0 then new:=0; break; fi;

fi;

if D2 = A then

B:=MyCutVector([a,0,0,d, 0,1,0,0, n,0,-1,x, -W*d,0,0,-a], 4)*One(F);
C:=MyCutVector([a,-W*d,0, d,-a,0, -n,W*x,-1], 3)*One(F);

D3:=B*A2*C^-1;

z1:=IntFFE(D3[3][1]);
z2:=IntFFE(D3[3][2]);

if z1+z2 = 0 then new:=0; break; fi;

fi;

od;
if new = 0 then break; fi;
od;

fi;

if new = 0 then break; fi;
od;

if new = 1 then
  Add(params,Linearise27b([y1,y2,0,0,0,0,0,1,1,y5,y6,0]));
fi;

od;

l := (P^4+P^3+6*P^2+9*P+13)/2;
if Length(params) <> l then Error("wrong number"); fi;
return params;
end );

