//Descendants of 6.150, Cases 3,4,5,8
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

Z:=Integers();
V1:=VectorSpace(F,1);
V3:=VectorSpace(F,3);
H31:=Hom(V3,V1);
H33:=Hom(V3,V3);

SQ:={};
for i in [1..((p-1) div 2)] do
  Include(~SQ,i^2 mod p);
end for;
for i in [2..p-1] do
  if i notin SQ then lns:=i; break; end if;
end for;

r:=(p+1) div 2;

gtotal:=0;

GR:=[];
//These store representative parameter sets (x,y,z) in cases 3,4,5
mats3:=[];
mats4:=[];
mats5:=[];
//This stores representative parameter sets (x,y,z,t) in case 8
mats8:=[];

/*
Descendants of 6.150

Case 3: baaa=baab=babb, ca-baa=cb=0
*/

count:=0;

for y1 in [0..p-1] do
for y2 in [0..p-1] do
range:=[0,1,lns];
if y1 gt 0 then range:=[0]; end if;
for y3 in range do

A:=H31![y1,y2,y3];

if A eq 0 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-1, (b,a,b,b)*(b,a,a,a)^-1, 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a), (c,b), a^p, b^p, c^p>;
  Append(~mats3,[0,0,0]);
  continue;
end if;

new:=1;
index:=p^2*y3+p*y1+y2;

for a in [1..p-1] do
for c in [0..p-1] do

B:=H33![a,0,c,
        0,a,-c,
        0,0,a^2];

C:=F!(a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z3+p*z1+z2;

if ind1 lt index then new:=0; break; end if;

B:=H33![0,a,c,
        a,0,-c,
        0,0,a^2];

C:=F!(-a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z3+p*z1+z2;

if ind1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-1, (b,a,b,b)*(b,a,a,a)^-1, 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a), (c,b), 
                                   a^p=(b,a,a,a)^y1, b^p=(b,a,a,a)^y2, c^p=(b,a,a,a)^y3>;
  Append(~mats3,[y1,y2,y3]);
end if;

end for;
end for;
end for;

gtotal:=gtotal+count;
print count,#mats3,(p+1+(p+3)*Gcd(p-1,3)+Gcd(p-1,4))/2;

/*
Descendants of 6.150

Case 4: babb+baaa=baab=ca-baa=cb=0
*/

count:=0;

for y1 in [0..p-1] do
for y2 in [0..p-1] do
range:=[0,1,lns];
if y1+y2 gt 0 then range:=[0]; end if;
for y3 in range do

A:=H31![y1,y2,y3];

if A eq 0 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,b,b)*(b,a,a,a), (b,a,a,b), 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a)^r, (c,b), a^p, b^p, c^p>;
  Append(~mats4,[0,0,0]);
  continue;
end if;

new:=1;
index:=p^2*y1+p*y2+y3;

for a in [1..p-1] do
for b in [1,-1] do
c:=0; e:=0;

B:=H33![a*b,0,c,
        0,a,e,
        0,0,a^2];

C:=F!(a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

B:=H33![0,a*b,c,
        a,0,e,
        0,0,a^2];

C:=F!(a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,b,b)*(b,a,a,a), (b,a,a,b), 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a)^r, (c,b), 
                                   a^p=(b,a,a,a)^y1, b^p=(b,a,a,a)^y2, c^p=(b,a,a,a)^y3>;
  Append(~mats4,[y1,y2,y3]);
end if;

end for;
end for;
end for;

gtotal:=gtotal+count;
print count,#mats4,3+Gcd(p-1,3)*(p+3+Gcd(p-1,4))/4;

/*
Descendants of 6.150

Case 5: babb+wbaaa=baab=ca-baa=cb=0
*/

count:=0;

for y1 in [0..p-1] do
for y2 in [0..p-1] do
range:=[0,1];
if y1+y2 gt 0 then range:=[0]; end if;
for y3 in range do

A:=H31![y1,y2,y3];

if A eq 0 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,b,b)*(b,a,a,a)^w, (b,a,a,b), 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a)^r, (c,b), a^p, b^p, c^p>;
  Append(~mats5,[0,0,0]);
  continue;
end if;

new:=1;
index:=p^2*y1+p*y2+y3;

for a in [1..p-1] do
for b in [1,-1] do
c:=0; e:=0;

B:=H33![a*b,0,c,
        0,a,e,
        0,0,a^2];

C:=F!(a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

B:=H33![0,a*b,c,
        w*a,0,e,
        0,0,w*a^2];

C:=F!(w^2*a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,b,b)*(b,a,a,a)^w, (b,a,a,b), 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a)^r, (c,b),
                                   a^p=(b,a,a,a)^y1, b^p=(b,a,a,a)^y2, c^p=(b,a,a,a)^y3>;
  Append(~mats5,[y1,y2,y3]);
end if;

end for;
end for;
end for;

gtotal:=gtotal+count;
print count,#mats5,2+Gcd(p-1,3)*(p+7-Gcd(p-1,4))/4;


/*
Descendants of 6.150

Case 8: baaa=baab, babb=xbaaa, ca-baa=cb=0 (2<=x<p)
*/

count:=0;

for x in [2..p-1] do

for y1 in [0..p-1] do
for y2 in [0..p-1] do
range:=[0,1,lns];
if y1+y2 gt 0 then range:=[0]; end if;
for y3 in range do

A:=H31![y1,y2,y3];

if A eq 0 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-1, (b,a,b,b)*(b,a,a,a)^-x, 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a), (c,b), a^p, b^p, c^p>;
  Append(~mats8,[x,0,0,0]);
  continue;
end if;

new:=1;
index:=p^2*y1+p*y2+y3;

for a in [1..p-1] do

B:=H33![a,0,0,
        0,a,0,
        0,0,a^2];

C:=F!(a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

B:=H33![0,a,0,
        x*a,0,0,
        0,0,-x*a^2];

C:=F!(x^2*a^4);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

end for;

if new eq 1 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-1, (b,a,b,b)*(b,a,a,a)^-x, 
                                  (c,a)*(b,a,a)^-1*(b,a,a,a), (c,b), 
                                   a^p=(b,a,a,a)^y1, b^p=(b,a,a,a)^y2, c^p=(b,a,a,a)^y3>;
  Append(~mats8,[x,y1,y2,y3]);
end if;

end for;
end for;
end for;

end for;

gtotal:=gtotal+count;
print count,#mats8,(5*p-7+(p^2-5)*Gcd(p-1,3)-Gcd(p-1,4))/2;

