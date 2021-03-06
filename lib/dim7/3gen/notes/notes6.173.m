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

SQ:={};
for i in [0..((p-1) div 2)] do
  Include(~SQ,i^2 mod p);
end for;

range:=[];
for i in [0..p-1] do
  Append(~range,[1,i]);
end for;
Append(~range,[0,1]);

//Descendants of 6.173
//<a,b,c|baab-ubaaa,babb-vbaaa,ca-bab,cb-wbaa,pa-ybaaa,pb-zbaaa,pc-tbaaa,class=4>
//A representative set of parameters [u,v,y,z,t] is stored in the following list.
parameters:=[];

Z:=Integers();
V1:=VectorSpace(F,1);
V3:=VectorSpace(F,3);
H31:=Hom(V3,V1);
H33:=Hom(V3,V3);

GR:=[];
gtotal:=0;

/*
Output from Baker-Campbell-Hausdorff
GR[]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-x, (b,a,b,b)*(b,a,a,a)^-y, 
(c,a)*(b,a,b)^-1*(b,a,a,a)^y, (c,b)*(b,a,a)^-w*(b,a,a,a)^w, a^p, b^p, c^p>;
*/


/*
Descendants of 6.173
Case 1: l=m=0
*/

count:=0;

for y3 in [0..p-1] do
for y1 in [0..p-1] do
for y2 in [0..p-1] do
if y3 gt 0 and y2 gt 0 then continue; end if;

A:=H31![y1,y2,y3];

if A eq 0 then
  count:=count+1;
  GR[count]:=Group<a,b,c | (b,a,a,b), (b,a,b,b), (c,a)*(b,a,b)^-1,
                           (c,b)*(b,a,a)^-w*(b,a,a,a)^w, a^p, b^p, c^p>;
  Append(~parameters,[0,0,0,0,0]);
  continue;
end if;

new:=1;
index:=p^2*y3+p*y1+y2;

for a in [1..p-1] do
for b in [1,-1] do
//for e in [0..p-1] do
e:=0;

  B:=H33![a,0,0,
          0,a*b,e,
          0,0,a^2];

  C:=F!(a^4*b);

  D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z3+p*z1+z2;

if ind1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;
//if new eq 0 then break; end if;
//end for;

if new eq 1 then
  count:=count+1;
  GR[count]:=Group<a,b,c | (b,a,a,b), (b,a,b,b), (c,a)*(b,a,b)^-1,
    (c,b)*(b,a,a)^-w*(b,a,a,a)^w, a^p=(b,a,a,a)^y1, b^p=(b,a,a,a)^y2, c^p=(b,a,a,a)^y3>;
  Append(~parameters,[0,0,y1,y2,y3]);
end if;

end for;
end for;
end for;

//print count,(p+1+(p+3)*Gcd(p-1,3)+Gcd(p-1,4))/2;
gtotal:=count;

/*
Descendants of 6.173

Cases other than l=m=0
*/

//first get p+1 orbits for x,y

mats:=[];

count:=0;
for l in [0..p-1] do
for m in [0..p-1] do
if l^2 mod p eq m then continue; end if;
//The orbit of l,l^2 contains 0,0

new:=1;
index:=p*l+m;
for r in range do
a:=r[1]; b:=r[2];
x:=F!(a^2+2*a*b*l+b^2*m);
if x eq 0 then continue; end if;

y:=F!(w*a*b+a^2*l+w*b^2*l+a*b*m);
z:=F!(w^2*b^2+2*w*a*b*l+a^2*m);
l1:=y*x^-1; l2:=Z!l1; m1:=z*x^-1; m1:=Z!m1;

ind1:=p*l2+m1;
if ind1 lt index then new:=0; break; end if;

l2:=Z!(-l1);

ind1:=p*l2+m1;
if ind1 lt index then new:=0; break; end if;

end for;

if new eq 1 then
  count:=count+1;
  //print l,m;
  Append(~mats,[l,m]);  
end if;


end for;
end for;

//print count,p+1;

count:=0;

for xy in mats do
x:=xy[1]; y:=xy[2];

//y1,y2,y3 represent pa,pb,pc

for y1 in [0..p-1] do
for y2 in [0..p-1] do
yrange:=[0..p-1];
if y1+y2 gt 0 then yrange:=[0]; end if;
for y3 in yrange do

new:=1;
A:=H31![y1,y2,y3];
index:=p^2*y1+p*y2+y3;

if A eq 0 then
  count:=count+1;
  GR[gtotal+count]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-x, (b,a,b,b)*(b,a,a,a)^-y, 
    (c,a)*(b,a,b)^-1*(b,a,a,a)^y, (c,b)*(b,a,a)^-w*(b,a,a,a)^w, a^p, b^p, c^p>;
  Append(~parameters,[x,y,0,0,0]);
  continue;
end if;

//Look for matrices of + type which fix x,y

for r in range do
a:=r[1]; b:=r[2];

z:=F!(a^2+2*a*b*x+b^2*y);
if z eq 0 then continue; end if;

x1:=F!(w*a*b+a^2*x+w*b^2*x+a*b*y)*z^-1; x1:=Z!x1;
if x ne x1 then continue; end if;
ynew:=F!(w^2*b^2+2*w*a*b*x+a^2*y)*z^-1; ynew:=Z!ynew;
if y ne ynew then continue; end if;

for c in [1..p-1] do

B:=H33![a,b,0,
        w*b,a,0,
        0,0,(a^2-w*b^2)*c];

C:=F!((a^2-w*b^2)*(a^2+2*a*b*x+b^2*y)*c^3);

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[2][1]);
z3:=Z!(D[3][1]);

ind1:=p^2*z1+p*z2+z3;

if ind1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 0 then continue; end if;

//Look for matrices of - type which fix x,y

for r in range do
a:=r[1]; b:=r[2];

z:=F!(a^2+2*a*b*x+b^2*y);
if z eq 0 then continue; end if;

x1:=-F!(w*a*b+a^2*x+w*b^2*x+a*b*y)*z^-1; x1:=Z!x1;
if x ne x1 then continue; end if;
ynew:=F!(w^2*b^2+2*w*a*b*x+a^2*y)*z^-1; ynew:=Z!ynew;
if y ne ynew then continue; end if;

for c in [1..p-1] do

B:=H33![a,b,0,
        -w*b,-a,0,
        0,0,(a^2-w*b^2)*c];

C:=-F!((a^2-w*b^2)*(a^2+2*a*b*x+b^2*y)*c^3);

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
  GR[gtotal+count]:=Group<a,b,c | (b,a,a,b)*(b,a,a,a)^-x, (b,a,b,b)*(b,a,a,a)^-y, 
                    (c,a)*(b,a,b)^-1*(b,a,a,a)^y, (c,b)*(b,a,a)^-w*(b,a,a,a)^w,
                    a^p=(b,a,a,a)^y1, b^p=(b,a,a,a)^y2, c^p=(b,a,a,a)^y3>;
  Append(~parameters,[x,y,y1,y2,y3]);end if;

end for;
end for;
end for;

end for;

//print count,(5*p+5+(p^2+p)*Gcd(p-1,3)-Gcd(p-1,4))/2;

gtotal:=gtotal+count;

print "Algebra 6.173 has",gtotal,"descendants of order p^7 and class 4";
print "3p + 3 + (p^2+2p+3)gcd(p-1,3)/2 =",3*p+3+(p^2+2*p+3)*Gcd(p-1,3)/2;
print "List parameters has size",#parameters;
