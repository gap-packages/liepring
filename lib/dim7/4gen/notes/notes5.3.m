//Descendants of 5.3 of order p^7

readi p,"Input the prime p";
//Get a primitive element
w:=0;
F:=FiniteField(p);
for i in [2..p-1] do
a:=F!i;
S:={a};
for j in [2..p-1] do
  Include(~S,a^j);
end for;
if #S eq p-1 then
  w:=i;
  break;
end if;
end for;

print "p equals",p;
print "w equals",w;

tt:=Cputime();

/*
Case 4 in the descendants of 5.3

We have a 2x2 matrix A representing pa,pb and we
map it to (det P)^-1.PAP^-1 where

  P =   [a ,b]
     +/-[wb,a]

Get orbits of rank two matrices A

Store orbit representatives in params4
*/

params4:=[];
count:=0;

Z:=Integers();
V2:=VectorSpace(F,2);
H22:=Hom(V2,V2);

range:={[0,1]};
for i in [0..p-1] do
  Include(~range,[1,i]);
end for;

SQ:={};
for i in [1..((p-1) div 2)] do
  Include(~SQ,F!(i^2));
end for;

for i in [2..p-1] do
  if F!i notin SQ then
    lnsq:=F!i;
    zlnsq:=i;
    break;
  end if;
end for;

for y1 in [0,1,zlnsq] do
for y2 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

//A represents pa,pb
A:=H22![y1,y2,y3,y4];

if Rank(A) eq 2 then

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

//Find first non-zero entry in A
nza:=1;
u:=A[1][1];
if u eq 0 then
  nza:=2;
  u:=A[1][2];
end if;

if u eq 1 or u eq lnsq then

for r in range do
a:=r[1]; b:=r[2];

P:=H22![a,b,w*b,a];
c:=F!(a^2-w*b^2);

//D is the image of A under the action of the group element
D:=c^-1*P*A*P^-1;

//Find first non-zero entry in D
nzd:=1;
u:=D[1][1];
if u eq 0 then
  nzd:=2;
  u:=D[1][2];
end if;

if nza lt nzd then new:=0; end if;

if nza eq nzd then

D:=u^-1*D;
if u notin SQ then D:=lnsq*D; end if;

//Get entries in image of A
z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

P[2]:=-P[2];
c:=-c;

//D is the image of A under the action of the group element
D:=c^-1*P*A*P^-1;

//Find first non-zero entry in D
nzd:=1;
u:=D[1][1];
if u eq 0 then
  nzd:=2;
  u:=D[1][2];
end if;

if nza lt nzd then new:=0; end if;

if nza eq nzd then

D:=u^-1*D;
if u notin SQ then D:=lnsq*D; end if;

//Get entries in image of A
z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  //We have a new orbit representative
  count:=count+1;
  Append(~params4,[y1,y2,y3,y4]);
end if;

end if;
end if;

end for;
end for;
end for;
end for;

print count,#params4;
print "p^2+(p+1-gcd(p-1,4))/2 =",p^2+(p+1-GCD(p-1,4))/2;
print "";

/*
This program computes the orbits of the 4x2 matrices
representing a^p,b^p,c^p,d^p needed for Case 6 in the
descendants of 5.3

[c,a]=[b,a,b], [c,b]=[b,a,a], [d,a]=1, [d,b]=[b,a,b], [d,c]=1

We have a 4x2 matrix A with rows representing pa,pb,pc,pd
and we premultiply by

     [a,-b,c,d]
  +/-[b,a,l,m]
     [0,0,a^2-b^2,-4ab]
  +/-[0,0,ab,a^2-b^2]

and post multiply by the inverse of

  (a^2+b^2). +/-[a,-b]
                [b,a]

We store the orbit representatives in params6, as vectors [y1,y2,y3,y4,y5,y6,y7,y8]
where row 1 of A is [y1,y2], row 2 is [y3,y4] etc.
*/

count:=0;
params6:=[];

Z:=Integers();
V2:=VectorSpace(F,2);
H22:=Hom(V2,V2);

range:={[0,1]};
for i in [0..p-1] do
  if (1+i^2) mod p ne 0 then
    Include(~range,[1,i]);
  end if;
end for;

nz:=[[1,1],[1,2],[2,1],[2,2]];

SQ:={};
for i in [1..((p-1) div 2)] do
  Include(~SQ,(F!i)^2);
end for;

for i in [2..p-1] do
  j:=F!i;
  if j notin SQ then
    lnsq:=j; zlnsq:=i; break;
  end if;
end for;

mats:={};
L:=[]; Nmr := 0;

//First get non-zero possibilities for pc,pd.

for y1 in [0..1] do
for y2 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y1,y2,y3,y4];
//Don't keep zero matrix
if A eq 0 then continue; end if;

//Get first non-zero entry in A
for i in [1..4] do
  spot:=nz[i];
  if A[spot[1],spot[2]] ne 0 then break; end if;
end for;

//We only need A with leading entry 1
if A[spot[1],spot[2]] ne 1 then continue; end if;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a^2-b^2,-4*a*b,a*b,a^2-b^2];
C:=F!(a^2+b^2)*H22![a,-b,b,a];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then; break; end if;
end for;

//Normalize D
u:=D[spot[1],spot[2]];
D:=u^-1*D;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

B[2]:=-B[2];
C[1]:=-C[1];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then; break; end if;
end for;

//Normalize D
u:=D[spot[1],spot[2]];
D:=u^-1*D;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  Include(~mats,A);
end if;

end for;
end for;
end for;
end for;

count:=0;

//Get pa,pb when pc=pd=0
for y1 in [0,1,zlnsq] do
for y2 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A:=H22![y1,y2,y3,y4];
if A eq 0 then
  count:=count+1;
  Append(~params6,[0,0,0,0,0,0,0,0]);
  continue;
end if;

//Get first non-zero entry in A
for i in [1..4] do
  spot:=nz[i];
  if A[spot[1],spot[2]] ne 0 then break; end if;
end for;

if (A[spot[1],spot[2]] ne 1) and (A[spot[1],spot[2]] ne zlnsq) then continue; end if;

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a,-b,b,a];
C:=F!(a^2+b^2)*H22![a,-b,b,a];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then break; end if;
end for;

u:=D[spot[1],spot[2]];
D:=u^-1*D;
if u notin SQ then D:=lnsq*D; end if;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

B[2]:=-B[2];
C[1]:=-C[1];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then break; end if;
end for;

u:=D[spot[1],spot[2]];
D:=u^-1*D;
if u notin SQ then D:=lnsq*D; end if;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  Append(~params6,[y1,y2,y3,y4,0,0,0,0]);
end if;

end for;
end for;
end for;
end for;

//Now get rest

for A in mats do

t1:=Z!(A[1][1]);
t2:=Z!(A[1][2]);
t3:=Z!(A[2][1]);
t4:=Z!(A[2][2]);

if Rank(A) eq 2 then
  count:=count+1;
  Append(~params6,[0,0,0,0,t1,t2,t3,t4]);
else;

if t1 eq 0 and t3 eq 0 then

if t2 ne 0 then
  v:=A[1];
else;
  v:=A[2];
end if;
u:=v[2];
v:=u^-1*v;

for y1 in [0..p-1] do
y2:=0;
for y3 in [0..p-1] do
y4:=0;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A2:=H22![y1,y2,y3,y4];

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a^2-b^2,-4*a*b,a*b,a^2-b^2];
C:=F!(a^2+b^2)*H22![a,-b,b,a];
D1:=B*A*C^-1;

//Get first non-zero entry in D1
for i in [1..4] do
  spot:=nz[i];
  if D1[spot[1],spot[2]] ne 0 then break; end if;
end for;

u1:=D1[spot[1],spot[2]];
D1:=u1^-1*D1;

B[2]:=-B[2];
C[1]:=-C[1];
D2:=B*A*C^-1;

//Get first non-zero entry in D2
for i in [1..4] do
  spot:=nz[i];
  if D2[spot[1],spot[2]] ne 0 then break; end if;
end for;

u2:=D2[spot[1],spot[2]];
D2:=u2^-1*D2;

if D1 eq A or D2 eq A then

B:=H22![a,-b,b,a];
C:=F!(a^2+b^2)*H22![a,-b,b,a];

if D1 eq A then

C1:=u1^2*C;
D:=B*A2*C1^-1;

z2:=D[1][2];
D[1]:=D[1]-z2*v;
z4:=D[2][2];
D[2]:=D[2]-z4*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

if D2 eq A then

B[2]:=-B[2];
C[1]:=-C[1];
C2:=u2^2*C;
D:=B*A2*C2^-1;

z2:=D[1][2];
D[1]:=D[1]-z2*v;
z4:=D[2][2];
D[2]:=D[2]-z4*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;
end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  Append(~params6,[y1,0,y3,0,0,t2,0,t4]);
end if;

end for;
end for;

else;
//one of t1,t3 is non-zero

if t1 ne 0 then
  v:=A[1];
else;
  v:=A[2];
end if;
u:=v[1];
v:=u^-1*v;

y1:=0;
for y2 in [0..p-1] do
y3:=0;
for y4 in [0..p-1] do

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A2:=H22![y1,y2,y3,y4];

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a^2-b^2,-4*a*b,a*b,a^2-b^2];
C:=F!(a^2+b^2)*H22![a,-b,b,a];
D1:=B*A*C^-1;

//Get first non-zero entry in D1
for i in [1..4] do
  spot:=nz[i];
  if D1[spot[1],spot[2]] ne 0 then break; end if;
end for;

u1:=D1[spot[1],spot[2]];
D1:=u1^-1*D1;

B[2]:=-B[2];
C[1]:=-C[1];
D2:=B*A*C^-1;

//Get first non-zero entry in D2
for i in [1..4] do
  spot:=nz[i];
  if D2[spot[1],spot[2]] ne 0 then break; end if;
end for;

u2:=D2[spot[1],spot[2]];
D2:=u2^-1*D2;

if D1 eq A or D2 eq A then

B:=H22![a,-b,b,a];
C:=F!(a^2+b^2)*H22![a,-b,b,a];

if D1 eq A then

C1:=u1^2*C;
D:=B*A2*C1^-1;

z1:=D[1][1];
D[1]:=D[1]-z1*v;
z3:=D[2][1];
D[2]:=D[2]-z3*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

if D2 eq A then

B[2]:=-B[2];
C[1]:=-C[1];
C2:=u2^2*C;

D:=B*A2*C2^-1;

z1:=D[1][1];
D[1]:=D[1]-z1*v;
z3:=D[2][1];
D[2]:=D[2]-z3*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  Append(~params6,[0,y2,0,y4,t1,t2,t3,t4]);
end if;

end for;
end for;

end if;
end if;
end for;

print "There are",count,"four generator p-class three rings of order";
print "p^7 satisfying ca=bab, cb=baa, da=0, db=bab, dc=0";

print "params6 has size",#params6;

if p mod 12 eq 1 then print "Total number of Lie rings is 3p^2+(19/2)p+15+p^3/2 =",3*p^2+(19/2)*p+15+p^3/2; end if;
if p mod 12 eq 5 then print "Total number of Lie rings is 3p^2+(19/2)p+12+p^3/2 =",3*p^2+(19/2)*p+12+p^3/2; end if;
if p mod 12 eq 7 then print "Total number of Lie rings is 2p^2+(5/2)p+2+p^3/2 =",2*p^2+(5/2)*p+2+p^3/2; end if;
if p mod 12 eq 11 then print "Total number of Lie rings is 2p^2+(5/2)p+3+p^3/2 =",2*p^2+(5/2)*p+3+p^3/2; end if;

print "";


/*
This program computes the orbits of the 4x2 matrices representing
a^p,b^p,c^p,d^p needed for Case 7 in the descendants of 5.3

[c,a]=[b,a,b], [c,b]=[b,a,a]^w, [d,a]=1, [d,b]=[b,a,b], [d,c]=1

We have a 4x2 matrix A with rows representing pa,pb,pc,pd
and we premultiply by

     [a,b,c,d]
  +/-[-wb,a,l,m]
     [0,0,a^2-wb^2,4wab]
  +/-[0,0,-ab,a^2-wb^2]

and post multiply by the inverse of

  (a^2+wb^2). +/-[a,b]
                 [-wb,a]

We store a set of orbit representatives in params7
*/

count:=0;
params7:=[];

Z:=Integers();
V2:=VectorSpace(F,2);
H22:=Hom(V2,V2);

range:={[0,1]};
for i in [0..p-1] do
  if (1+w*i^2) mod p ne 0 then
    Include(~range,[1,i]);
  end if;
end for;

nz:=[[1,1],[1,2],[2,1],[2,2]];

SQ:={};
for i in [1..((p-1) div 2)] do
  Include(~SQ,(F!i)^2);
end for;

for i in [2..p-1] do
  j:=F!i;
  if j notin SQ then
    lnsq:=j; zlnsq:=i; break;
  end if;
end for;

mats:={};

//First get non-zero possibilities for pc,pd.

for y1 in [0..1] do
for y2 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y1,y2,y3,y4];
//Don't store zero matrix
if A eq 0 then continue; end if;

//Get first non-zero entry in A
for i in [1..4] do
  spot:=nz[i];
  if A[spot[1],spot[2]] ne 0 then break; end if;
end for;

//We only need A with leading entry 1
if A[spot[1],spot[2]] ne 1 then continue; end if;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a^2-w*b^2,4*w*a*b,-a*b,a^2-w*b^2];
C:=F!(a^2+w*b^2)*H22![a,b,-w*b,a];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then; break; end if;
end for;

//Normalize D
u:=D[spot[1],spot[2]];
D:=u^-1*D;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

B[2]:=-B[2];
C[1]:=-C[1];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then; break; end if;
end for;

//Normalize D
u:=D[spot[1],spot[2]];
D:=u^-1*D;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  Include(~mats,A);
end if;

end for;
end for;
end for;
end for;

count:=0;
L:=[]; Nmr := 0;

//Get pa,pb when pc=pd=0
for y1 in [0,1,zlnsq] do
for y2 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A:=H22![y1,y2,y3,y4];
if A eq 0 then
  count:=count+1;
  Append(~params7,[0,0,0,0,0,0,0,0]);
  continue;
end if;

//Get first non-zero entry in A
for i in [1..4] do
  spot:=nz[i];
  if A[spot[1],spot[2]] ne 0 then break; end if;
end for;

if (A[spot[1],spot[2]] ne 1) and (A[spot[1],spot[2]] ne zlnsq) then continue; end if;

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a,b,-w*b,a];
C:=F!(a^2+w*b^2)*H22![a,b,-w*b,a];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then break; end if;
end for;

u:=D[spot[1],spot[2]];
D:=u^-1*D;
if u notin SQ then D:=lnsq*D; end if;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

B[2]:=-B[2];
C[1]:=-C[1];
D:=B*A*C^-1;

//Get first non-zero entry in D
for i in [1..4] do
  spot:=nz[i];
  if D[spot[1],spot[2]] ne 0 then break; end if;
end for;

u:=D[spot[1],spot[2]];
D:=u^-1*D;
if u notin SQ then D:=lnsq*D; end if;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  Append(~params7,[y1,y2,y3,y4,0,0,0,0]);
end if;

end for;
end for;
end for;
end for;

//Now get rest

for A in mats do

t1:=Z!(A[1][1]);
t2:=Z!(A[1][2]);
t3:=Z!(A[2][1]);
t4:=Z!(A[2][2]);

if Rank(A) eq 2 then
  count:=count+1;
  Append(~params7,[0,0,0,0,t1,t2,t3,t4]);
else;

if t1 eq 0 and t3 eq 0 then

if t2 ne 0 then
  v:=A[1];
else;
  v:=A[2];
end if;
u:=v[2];
v:=u^-1*v;

for y1 in [0..p-1] do
y2:=0;
for y3 in [0..p-1] do
y4:=0;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A2:=H22![y1,y2,y3,y4];

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a^2-w*b^2,4*w*a*b,-a*b,a^2-w*b^2];
C:=F!(a^2+w*b^2)*H22![a,b,-w*b,a];
D1:=B*A*C^-1;

//Get first non-zero entry in D1
for i in [1..4] do
  spot:=nz[i];
  if D1[spot[1],spot[2]] ne 0 then break; end if;
end for;

u1:=D1[spot[1],spot[2]];
D1:=u1^-1*D1;

B[2]:=-B[2];
C[1]:=-C[1];
D2:=B*A*C^-1;

//Get first non-zero entry in D2
for i in [1..4] do
  spot:=nz[i];
  if D2[spot[1],spot[2]] ne 0 then break; end if;
end for;

u2:=D2[spot[1],spot[2]];
D2:=u2^-1*D2;

if D1 eq A or D2 eq A then

B:=H22![a,b,-w*b,a];
C:=F!(a^2+w*b^2)*H22![a,b,-w*b,a];

if D1 eq A then

C1:=u1^2*C;
D:=B*A2*C1^-1;

z2:=D[1][2];
D[1]:=D[1]-z2*v;
z4:=D[2][2];
D[2]:=D[2]-z4*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

if D2 eq A then

B[2]:=-B[2];
C[1]:=-C[1];
C2:=u2^2*C;
D:=B*A2*C2^-1;

z2:=D[1][2];
D[1]:=D[1]-z2*v;
z4:=D[2][2];
D[2]:=D[2]-z4*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;
end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  Append(~params7,[y1,0,y3,0,0,t2,0,t4]);
end if;

end for;
end for;

else;
//one of t1,t3 is non-zero

if t1 ne 0 then
  v:=A[1];
else;
  v:=A[2];
end if;
u:=v[1];
v:=u^-1*v;

y1:=0;
for y2 in [0..p-1] do
y3:=0;
for y4 in [0..p-1] do

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A2:=H22![y1,y2,y3,y4];

for r in range do
a:=r[1]; b:=r[2];

B:=H22![a^2-w*b^2,4*w*a*b,-a*b,a^2-w*b^2];
C:=F!(a^2+w*b^2)*H22![a,b,-w*b,a];
D1:=B*A*C^-1;

//Get first non-zero entry in D1
for i in [1..4] do
  spot:=nz[i];
  if D1[spot[1],spot[2]] ne 0 then break; end if;
end for;

u1:=D1[spot[1],spot[2]];
D1:=u1^-1*D1;

B[2]:=-B[2];
C[1]:=-C[1];
D2:=B*A*C^-1;

//Get first non-zero entry in D2
for i in [1..4] do
  spot:=nz[i];
  if D2[spot[1],spot[2]] ne 0 then break; end if;
end for;

u2:=D2[spot[1],spot[2]];
D2:=u2^-1*D2;

if D1 eq A or D2 eq A then

B:=H22![a,b,-w*b,a];
C:=F!(a^2+w*b^2)*H22![a,b,-w*b,a];

if D1 eq A then

C1:=u1^2*C;
D:=B*A2*C1^-1;

z1:=D[1][1];
D[1]:=D[1]-z1*v;
z3:=D[2][1];
D[2]:=D[2]-z3*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

if D2 eq A then

B[2]:=-B[2];
C[1]:=-C[1];
C2:=u2^2*C;

D:=B*A2*C2^-1;

z1:=D[1][1];
D[1]:=D[1]-z1*v;
z3:=D[2][1];
D[2]:=D[2]-z3*v;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  Append(~params7,[0,y2,0,y4,t1,t2,t3,t4]);
end if;

end for;
end for;

end if;
end if;
end for;

print "There are",count,"four generator p-class three Lie rings of order";
print "p^7 satisfying ca=bab, cb=wbaa da=0, db=bab, dc=0";

print "params7 has size",#params7;

if p mod 12 eq 1 then print "Total number of Lie rings is 2p^2+(5/2)p+2+p^3/2 =",2*p^2+(5/2)*p+2+p^3/2; end if;
if p mod 12 eq 5 then print "Total number of Lie rings is 2p^2+(5/2)p+3+p^3/2 =",2*p^2+(5/2)*p+3+p^3/2; end if;
if p mod 12 eq 7 then print "Total number of Lie rings is 3p^2+(19/2)p+13+p^3/2 =",3*p^2+(19/2)*p+13+p^3/2; end if;
if p mod 12 eq 11 then print "Total number of Lie rings is 3p^2+(19/2)p+10+p^3/2 =",3*p^2+(19/2)*p+10+p^3/2; end if;

print "";
print Cputime(tt);
