//This is the new p^5 version
//This program computes presentations for Lie rings from
//Case 6 in the descendants of 4.1 with L^2 of order p^3.
// da=0, db=wca, dc=ba

//Orbit representatives for pa,pb,pc,pd are stored in paramsCase6
//as 12 vectors, with the first three entries representing pa, etc.

paramsCase6:=[];

readi p,"Input the prime p";
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

//Get the squares, and the least non-square lns
SQ:={};
for x in [1..((p-1) div 2)] do
  y:=x^2 mod p;
  Include(~SQ,y);
end for;
for i in [2..p-1] do
  if i notin SQ then lns:=i; break; end if;
end for;

//Get the leading entry in a matrix
leading:=function(A)
  if A[1][1] ne 0 then return A[1][1]; end if;
  if A[1][2] ne 0 then return A[1][2]; end if;
  if A[2][1] ne 0 then return A[2][1]; end if;
  return A[2][2];
end function;

tt:=Cputime();

Z:=Integers();
V2:=VectorSpace(F,2);
V3:=VectorSpace(F,3);
V4:=VectorSpace(F,4);
H22:=Hom(V2,V2);
H43:=Hom(V4,V3);
H33:=Hom(V3,V3);
H44:=Hom(V4,V4);

range:={[0,1]};
for i in [0..p-1] do
  Include(~range,[1,i]);
end for;

GR:=[];

//Get pb,pc when pa=pd=0.

//If pb,pc notin <ba,ca> then we can assume that
//pc=cb and pb=0 or ca
//GR[1]:=Group<a,b,c,d|(d,a),(d,b)=(c,a)^w,(d,c)=(b,a),
  //           a^p,b^p,c^p=(c,b),d^p>;
Append(~paramsCase6,[0,0,0,0,0,0,0,0,1,0,0,0]);
//GR[2]:=Group<a,b,c,d|(d,a),(d,b)=(c,a)^w,(d,c)=(b,a),
  //           a^p,b^p=(c,a),c^p=(c,b),d^p>;
Append(~paramsCase6,[0,0,0,0,1,0,0,0,1,0,0,0]);

//Get pb,pc when pa=pd=0 and pb,pc in <ba,ca>
//Can assume that pb=0 or ca
y1:=0;
for y2 in [0,1] do
y3range:=[0..p-1];
if y2 eq 0 then y3range:=[0,1]; end if;
for y3 in y3range do
y4range:=[0..p-1];
if y1 eq 0 and y2 eq 0 then y4range:=[0,1]; end if;
for y4 in y4range do

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A:=H22![y1,y2,y3,y4];

if A eq 0 then
  Append(~paramsCase6,[0,0,0,0,0,0,0,0,0,0,0,0]);
  continue;
end if;

for r in range do
a:=r[1];
d:=r[2];
for s in range do
b:=s[1];
c:=s[2];

B:=H22![c,w*b,b,c];
C:=H22![a*c-w*b*d,w*a*b-w*c*d,
        a*b-c*d,a*c-w*b*d];

D:=B*A*C^-1;

k:=leading(D);
D:=k^-1*D;


z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; break; end if;

B[1]:=-B[1];
C[1]:=-C[1];

D:=B*A*C^-1;

k:=leading(D);
D:=k^-1*D;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~paramsCase6,[0,0,0,0,y2,0,y3,y4,0,0,0,0]);
end if;

end for;
end for;
end for;


//Get pb,pc when pa=0, pd=ca

for z in [1..((p-1) div 2)] do
for u in [0,1] do
for t in [0..p-1] do
  Append(~paramsCase6,[0,0,0,0,0,z,u,0,t,0,1,0]);
end for;
end for;
end for;

for t in [1..p-1] do
for x in [0,1] do
  Append(~paramsCase6,[0,0,0,x,0,0,0,0,t,0,1,0]);
end for;
end for;

Append(~paramsCase6,[0,0,0,0,0,0,0,0,0,0,1,0]);

Append(~paramsCase6,[0,0,0,0,0,0,1,0,0,0,1,0]);

for u in [0..((p-1) div 2)] do
  Append(~paramsCase6,[0,0,0,1,0,0,u,0,0,0,1,0]);
end for;


//Get pb,pc when pa=0, pd=cb

for y1 in [0,1,lns] do
for y2 in [0..(p-1) div 2] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y1,y2,y3,y4];

if A eq 0 then
  Append(~paramsCase6,[0,0,0,0,0,0,0,0,0,0,0,1]);
  continue;
end if;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

for r in range do
b:=r[1]; c:=r[2];

B:=H22![c,w*b,b,c];

C:=H22![c*(c^2-w*b^2),w*b*(c^2-w*b^2),
        b*(c^2-w*b^2),c*(c^2-w*b^2)];

D:=B*A*C^-1;

k:=leading(D);
D:=k^-1*D;
if not IsSquare(k) then D:=lns*D; end if;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; break; end if;

B:=H22![-c,-w*b,b,c];

C:=H22![-c*(c^2-w*b^2),-w*b*(c^2-w*b^2),
        b*(c^2-w*b^2),c*(c^2-w*b^2)];

D:=B*A*C^-1;

k:=leading(D);
D:=k^-1*D;
if not IsSquare(k) then D:=lns*D; end if;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; break; end if;

end for;

if new eq 1 then
  Append(~paramsCase6,[0,0,0,y1,y2,0,y3,y4,0,0,0,1]);
end if;

end for;
end for;
end for;
end for;


//Get pb,pc when pa=ca, pd=cb

for x in [0..p-1] do
for u in [0..((p-1) div 2)] do
for v in [0..p-1] do
lastt:=p-1;
if u eq 0 then lastt:=(p-1) div 2; end if;
for t in [0..lastt] do
  Append(~paramsCase6,[0,1,0,x,0,0,u,v,t,0,0,1]);
end for;
end for;
end for;
end for;


//Get pa,pd when pa,pd span <ba,ca> and pb=pc=0

y1:=0;
y2:=1;
for y3 in [1..p-1] do
for y4 in [0..(p-1) div 2] do

new:=1;
index:=p*y3+y4;

A:=H22![y1,y2,y3,y4];

for r in range do
a:=r[1];
d:=r[2];
for s in range do
b1:=s[1];
c1:=s[2];
if F!(b1*w*y3*d^2+b1*a*y4*d+b1*a^2-y4*c1*d^2-d*a*c1-d*y3*a*c1) ne 0 then continue; end if;
k1:=F!(-b1*w*d*(y4*d+a+a*y3)+c1*(w*y3*d^2+y4*a*d+a^2));
k2:=F!((a^2-w*d^2)*(c1^2-w*b1^2));
k:=Z!(k1*k2^-1);
b:=k*b1;
c:=k*c1;

B:=H22![a,d,w*d,a];
C:=H22![a*c-w*b*d,w*a*b-w*c*d,
        a*b-c*d,a*c-w*b*d];

D:=B*A*C^-1;

z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~paramsCase6,[0,1,0,0,0,0,0,0,0,y3,y4,0]);
end if;

end for;
end for;


//Get pa,pd when pa,pd span <ba,ca>, pb=0 and pc notin <ba,ca>

//First get possibilities for pa,pd with pb=0, pc=cb
mats:={};

for y1 in [0..p-1] do
for y2 in [0..((p-1) div 2)] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A:=H22![y1,y2,y3,y4];
if Rank(A) ne 2 then continue; end if;

for r in range do
a:=r[1];
d:=r[2];

B:=H22![a,d,w*d,a];
C:=H22![a,-w*d,-d,a];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

B[2]:=-B[2];
C[2]:=-C[2];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~paramsCase6,[y1,y2,0,0,0,0,0,0,1,y3,y4,0]);
  Include(~mats,A);
end if;


end for;
end for;
end for;
end for;

//For the pa,pd just found, see which go with pb=0, pc=ca+cb
for A in mats do
y1:=Z!(A[1][1]);
y2:=Z!(A[1][2]);
y5:=Z!(A[2][1]);
y6:=Z!(A[2][2]);

y3:=0; y4:=1;

A2:=H43![y1,y2,0,0,0,0,y3,y4,1,y5,y6,0];
new:=1;

for r in range do
a:=r[1];
d:=r[2];

B:=H22![a,d,w*d,a];
C:=H22![a,-w*d,-d,a];
D:=B*A*C^-1;

B[2]:=-B[2];
C[2]:=-C[2];
D2:=B*A*C^-1;

if D eq A or D2 eq A then

for n in [0..p-1] do
for x in [0..p-1] do

if D eq A then

B:=H44![a,0,0,d,
        0,1,0,0,
        n,0,1,x,
        w*d,0,0,a];

C:=H33![a,-w*d,0,
        -d,a,0,
        -n,w*x,1];
D3:=B*A2*C^-1;

z1:=Z!(D3[3][1]);
z2:=Z!(D3[3][2]);

if z1+z2 eq 0 then new:=0; break; end if;

end if;

if D2 eq A then

B:=H44![a,0,0,d,
        0,1,0,0,
        n,0,-1,x,
        -w*d,0,0,-a];

C:=H33![a,-w*d,0,
        d,-a,0,
       -n,w*x,-1];

D3:=B*A2*C^-1;

z1:=Z!(D3[3][1]);
z2:=Z!(D3[3][2]);

if z1+z2 eq 0 then new:=0; break; end if;

end if;

end for;
if new eq 0 then break; end if;
end for;

end if;

if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~paramsCase6,[y1,y2,0,0,0,0,0,1,1,y5,y6,0]);
end if;

end for;


print "There are",#paramsCase6,"four generator p-class 2 Lie rings with";
print " da=0, db=wca, dc=ba";
print "(p^4 + p^3 + 6p^2 + 9p + 13)/2 =",(p^4+p^3+6*p^2+9*p+13)/2;

print Cputime(tt);
