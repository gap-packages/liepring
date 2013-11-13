readi p,"Input prime";

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

SQ:={};
for x in [0..p-1] do
  y:=x^2 mod p;
  Include(~SQ,y);
end for;

//Get least non-square
ns:=0;
for x in [2..p-1] do
  if x in SQ then continue; end if;
  ns:=x;
  break;
end for;

w2:=w^2 mod p;

Z:=Integers();
V2:=VectorSpace(F,2);
H22:=Hom(V2,V2);

//ca=bab, cb=baa
mats1:=[];
range1:=[0,1,ns];
range2:=[0,1];
y1range:=range1;
if p mod 4 ne 1 then y1range:=range2; end if;

for y1 in y1range do
y2range:=[0..p-1];
if y1 eq 0 then y2range:=range1; end if;
for y2 in y2range do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y1,y2,y3,y4];
if A eq 0 then Append(~mats1,[0,0,0,0]); continue; end if;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

for a in [0..p-1] do
for b in [0..p-1] do

e:=F!(a^2-b^2);
if e eq 0 then continue; end if;

P:=H22![a,b,b,a];

D:=e^-1*P*A*P^-1;

z1:=Z!(D[1,1]);
z2:=Z!(D[1,2]);
z3:=Z!(D[2,1]);
z4:=Z!(D[2,2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then
  new:=0;
end if;

P:=H22![a,b,-b,-a];

D:=-e^-1*P*A*P^-1;

z1:=Z!(D[1,1]);
z2:=Z!(D[1,2]);
z3:=Z!(D[2,1]);
z4:=Z!(D[2,2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then
  new:=0;
end if;

if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then 
  Append(~mats1,[y1,y2,y3,y4]);
end if;

end for;
end for;
end for;
end for;

print #mats1,p^2+(7*p+15)/2;

//ca=wbab, cb=baa
mats2:=[];

expect:=p^2+(3*p+5)/2;

for y1 in [0,1] do
y2range:=[0..p-1];
if y1 eq 0 then y2range:=[0,1,ns]; end if;
for y2 in y2range do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y1,y2,y3,y4];
if A eq 0 then Append(~mats2,[0,0,0,0]); continue; end if;

new:=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

for a in [0..p-1] do
for b in [0..p-1] do

e:=F!(a^2-w*b^2);
if e eq 0 then continue; end if;

P:=H22![a,w*b,b,a];

D:=e^-1*P*A*P^-1;

z1:=Z!(D[1,1]);
z2:=Z!(D[1,2]);
z3:=Z!(D[2,1]);
z4:=Z!(D[2,2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then
  new:=0;
end if;

P:=H22![a,w*b,-b,-a];

D:=-e^-1*P*A*P^-1;

z1:=Z!(D[1,1]);
z2:=Z!(D[1,2]);
z3:=Z!(D[2,1]);
z4:=Z!(D[2,2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then
  new:=0;
end if;

if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then 
  Append(~mats2,[y1,y2,y3,y4]);
end if;

end for;
end for;
end for;
end for;

print #mats2,p^2+(3*p+5)/2;

