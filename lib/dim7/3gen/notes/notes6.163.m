//Descendants of 6.163

//Modified version 10/8/2013 with p^5 work

//mats1 stores representative sets of parameters [u,x,y,z,t] for algebras
//<a,b,c|ca-baa,cb,pa-baa-ubab-ybaaa,pb+xbaa+bab-zbaaa,pc-tbaaa,class=4>

//mats2 stores representative sets of parameters [u,x,y,z,t] for algebras
//<a,b,c|ca-baa,cb,pa-baa-ubab-ybaaa,pb+xbaa+(u*x)bab-zbaaa,pc-tbaaa,class=4>

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
for i in [1..((p-1) div 2)] do
  Include(~SQ,i^2 mod p);
end for;
for i in [2..p-1] do
  if i notin SQ then lns:=i; break; end if;
end for;

Z:=Integers();

tt:=Cputime();

mats1:=[];
mats2:=[];


/*
Descendants of 6.163
<a,b,c|ca-baa,cb,pa-baa-ubab,pb+xbaa+bab,pc>
*/

mats:=[];


//First get possible u,x for isomorphism classes of daddys
// [1 u]
// [x 1]
for u in [1,lns] do
for x in [1..p-1] do

new:=1;
index:=p*u+x;

//Transformations of the first kind can only change matrix
//to one with higher index, so ignore them.
//Transformations of the second kind swap u and x, but you
//can divide one by a square and multiply the other by the
//same square.

if IsSquare(F!x) then
  u1:=1;
  x1:=(u*x) mod p;
else;
  u1:=lns;
  x1:=F!(u*x/lns);
  x1:=Z!x1;
end if;

if p*u1+x1 lt index then new:=0; end if;

if new eq 1 then
  //print u,x;
  Append(~mats,[u,x]);
end if;

end for;
end for;

//print #mats,3*(p-1)/2;

for A in mats do
u:=A[1]; x:=A[2];

yrange:=[0..(p-1) div 2];
zrange:=[0..p-1];
if u*x mod p ne 1 then yrange:=[0]; zrange:=[0..(p-1) div 2]; end if;
for t in [0..p-1] do
for y in yrange do
for z in zrange do

new:=1;
index:=p^2*t+p*y+z;

//transformations of type (*)
for a in [1,p-1] do
for c in [0..p-1] do
e:=-F!(c*u^-1+c*t);
if u*x mod p eq 1 then e:=F!(c*u^-1); end if;
y1:=F!(e+c*u^-1+a*y+c*t); y1:=Z!y1;
z1:=F!(x*e-x*c*u^-1-2*e*u^-1+z*a+e*t); z1:=Z!z1;
if p^2*t+p*y1+z1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 0 then continue; end if;

//transformations of type (**)
for a in [1..p-1] do
if u ne (a^2*x) mod p then continue; end if;
for c in [0..p-1] do
e:=F!(x*c-2*c*u^-1+z*a+c*t);
if u*x mod p eq 1 then e:=F!(c*u^-1); end if;

y1:=-F!(x*c-e-2*c*u^-1+z*a+c*t); y1:=Z!y1;
z1:=-F!(x*c*u^-1+e*u^-1+a^-1*y+e*t); z1:=Z!z1;
t1:=F!(-t+u^-1-x); t1:=Z!t1;
if p^2*t1+p*y1+z1 lt index then new:=0; break; end if;


end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~mats1,[u,x,y,z,t]);
end if;

end for;
end for;
end for;
end for;

print #mats1,2*p^2-(5*p-1)/2;

/*
Descendants of 6.163
second type
<a,b,c|ca-baa,cb,pa-baa-ubab,pb+xbaa+(u*x)bab,pc> u*x ne 1
*/

//First get possible u,x
mats:=[];

for u in [1,lns] do
for x in [1..p-1] do

//  [1  u ]
//  [x u*x]

if (u*x) mod p eq 1 then continue; end if;

new:=1;
index:=p*u+x;

//Transformations of the first kind can only change matrix
//to one with higher index, so ignore them.
//Transformations of the second kind swap u and x, but you
//can divide one by a square and multiply the other by the
//same square.

if IsSquare(F!x) then
  u1:=1;
  x1:=F!(u^-1*x^-1);
  x1:=Z!x1;
else;
  u1:=lns;
  x1:=F!(u^-1*x^-1*lns^-1);
  x1:=Z!x1;
end if;

if p*u1+x1 lt index then new:=0; end if;

if new eq 1 then
  Append(~mats,[u,x]);
  //print u,x;
end if;

end for;
end for;

//print #mats,p-3+Gcd(p-1,4)/2;


for A in mats do
u:=A[1]; x:=A[2];

for t in [0..p-1] do
for y in [0..(p-1) div 2] do
for z in [0..p-1] do

new:=1;
index:=p*y+z;

for a in [1,-1] do
for e in [0..p-1] do

y1:=F!(2*e+a*y+u*e*t); y1:=Z!y1;
z1:=F!(-2*x*e+a*z+e*t); z1:=Z!z1;
if p*y1+z1 lt index then new:=0; break; end if;

end for;
if new eq 0 then break; end if;
end for;

if new eq 1 and (u*x+1) mod p eq 0 and p mod 4 eq 1 then

  for a in [2..p-2] do
  if (a^2+1) mod p ne 0 then continue; end if;
  for e in [0..p-1] do

  y1:=-F!(2*e+z*a*u+u*e*t); y1:=Z!y1;
  z1:=F!(x*(2*e+a*y+u*e*t)); z1:=Z!z1;
  if p*y1+z1 lt index then new:=0; break; end if;

  end for;
  if new eq 0 then break; end if;
  end for;

end if;

if new eq 1 then
    Append(~mats2,[u,x,y,z,t]);
end if;

end for;
end for;
end for;
end for;

print #mats2,(p^3-5*p+p*Gcd(p-1,4))/2;
gtotal:=#mats1+#mats2;

print "Algebras 6.163 - 6.167 have",gtotal,"descendants of order p^7 and class 4";
print "(p^3 + 4p^2 - 10p + 1 + pgcd(p-1,4))/2 =",
(p^3+4*p^2-10*p+1+p*Gcd(p-1,4))/2;

print Cputime(tt);
