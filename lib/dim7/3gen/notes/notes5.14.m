//Immediate descendants of order p^7 of algebra 5.14
//Modified to cut down complexity 20 August 2013

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

//Get the squares, and the least non-square lns
SQ:={};
for x in [1..((p-1) div 2)] do
  y:=x^2 mod p;
  Include(~SQ,y);
end for;
for i in [2..p-1] do
  if i notin SQ then lns:=i; break; end if;
end for;

//Get cube roots of 1
CU:={1};
for x in [2..p-1] do
  if x^3 mod p eq 1 then Include(~CU,x); end if;
end for;

//If p=1 mod 3 compute a transversal for the cube roots of 1
transversal3:=[];
if p mod 3 eq 1 then
  transversal3:=[];
  sofar:={};
  for i in [1..p-1] do
    if i^3 mod p in sofar then continue; end if;
    Append(~transversal3,i);
    Include(~sofar,i^3 mod p);
  end for;
end if;

//If p=1 mod 4 then compute a transversal for the fourth roots of 1
transversal4:=[];
if p mod 4 eq 1 then
  sofar:={};
  for i in [1..(p-1) div 2] do
    if i^4 mod p in sofar then continue; end if;
    Append(~transversal4,i);
    Include(~sofar,i^4 mod p);
  end for;
end if;

//If p=1 mod 5 then compute a transversal for the fifth roots of 1
transversal5:=[];
if p mod 5 eq 1 then
  sofar:={};
  for i in [1..p-1] do
    if i^5 mod p in sofar then continue; end if;
    Append(~transversal5,i);
    Include(~sofar,i^5 mod p);
  end for;
end if;

//If p=1 mod 3 then compute a transversal for the sixth roots of 1
transversal6:=[];
if p mod 3 eq 1 then
  sofar:={};
  for i in [1..p-1] do
    if i^6 mod p in sofar then continue; end if;
    Append(~transversal6,i);
    Include(~sofar,i^6 mod p);
  end for;
end if;

Z:=Integers();
V1:=VectorSpace(F,1);
V2:=VectorSpace(F,2);
V3:=VectorSpace(F,3);
H12:=Hom(V1,V2);
H22:=Hom(V2,V2);
H32:=Hom(V3,V2);
H33:=Hom(V3,V3);

gtotal:=0;
GR:=[];
//These store the parameters in each of the 24 cases
//Currently they are stored as integer sequences of length 6,
//but an alternative (which takes less space) is to store them
//as vectors or matrices over GF(p).  Two have size p^5+.
params1:=[];
params2:=[];
params3:=[];
params4:=[];
params5:=[];
params6:=[];
params7:=[];
params8:=[];
params9:=[];
params10:=[];
params11:=[];
params12:=[];
params13:=[];
params14:=[];
params15:=[];
params16:=[];
params17:=[];
params18:=[];
params19:=[];
params20:=[];
params21:=[];
params22:=[];
params23:=[];
params24:=[];

tt:=Cputime();

/*
Descendants of 5.14.

Case 1: cb=caa=cab=cac=0
*/

Append(~params1,[0,0,0,0,0,0]);
Append(~params1,[1,0,0,0,0,0]);
Append(~params1,[0,1,0,0,0,0]);
Append(~params1,[0,w,0,0,0,0]);

Append(~params1,[0,0,1,0,0,0]);
Append(~params1,[0,0,w,0,0,0]);
Append(~params1,[0,1,1,0,0,0]);
Append(~params1,[0,1,w,0,0,0]);
Append(~params1,[0,w,1,0,0,0]);
Append(~params1,[0,w,w,0,0,0]);
for y in [0..p-1] do
  Append(~params1,[1,y,1,0,0,0]);
  Append(~params1,[1,y,w,0,0,0]);
end for;

Append(~params1,[1,1,0,1,0,0]);
Append(~params1,[1,w,0,1,0,0]);
for x in [0..p-1] do
  Append(~params1,[x,0,0,1,0,0]);
end for;

Append(~params1,[0,0,0,0,1,0]);
Append(~params1,[0,1,0,0,1,0]);
Append(~params1,[0,w,0,0,1,0]);
Append(~params1,[0,0,0,1,1,0]);
Append(~params1,[0,1,0,1,1,0]);
Append(~params1,[0,w,0,1,1,0]);

Append(~params1,[1,0,0,0,0,1]);
Append(~params1,[0,0,0,0,0,1]);
Append(~params1,[0,0,1,0,0,1]);
Append(~params1,[0,0,w,0,0,1]);

print #params1,3*p+22;
gtotal:=gtotal+#params1;


/*
Descendants of 5.14.

Case 2: caa=cab=cac=0, cb=baa
*/

Append(~params2,[0,0,0,0,0,0]);
Append(~params2,[1,0,0,0,0,0]);
Append(~params2,[0,1,0,0,0,0]);
Append(~params2,[0,w,0,0,0,0]);

Append(~params2,[0,0,1,0,0,0]);
Append(~params2,[0,0,w,0,0,0]);
Append(~params2,[0,1,1,0,0,0]);
Append(~params2,[0,1,w,0,0,0]);
Append(~params2,[0,w,1,0,0,0]);
Append(~params2,[0,w,w,0,0,0]);
for y in [0..p-1] do
  Append(~params2,[1,y,1,0,0,0]);
  Append(~params2,[1,y,w,0,0,0]);
end for;

Append(~params2,[1,1,0,1,0,0]);
Append(~params2,[1,w,0,1,0,0]);
for x in [0..p-1] do
  Append(~params2,[x,0,0,1,0,0]);
end for;

for x in [0..p-1] do
  Append(~params2,[0,x,0,0,1,0]);
  Append(~params2,[0,x,0,1,1,0]);
end for;

Append(~params2,[0,0,0,0,0,1]);
Append(~params2,[1,0,0,0,0,1]);
if p mod 3 eq 1 then
  Append(~params2,[w,0,0,0,0,1]);
  Append(~params2,[w^2,0,0,0,0,1]);
end if;
Append(~params2,[0,0,1,0,0,1]);
Append(~params2,[0,0,w,0,0,1]);
if p mod 4 eq 1 then
  Append(~params2,[0,0,w^2,0,0,1]);
  Append(~params2,[0,0,w^3,0,0,1]);
end if;

print #params2,5*p+13+Gcd(p-1,3)+Gcd(p-1,4);
gtotal:=gtotal+#params2;

/*
Descendants of 5.14.

Case 3: cb=bab=bac=cac=0
*/

Append(~params3,[0,0,0,0,0,0]);
Append(~params3,[1,0,0,0,0,0]);
Append(~params3,[0,0,0,0,1,0]);
Append(~params3,[0,1,0,0,1,0]);
Append(~params3,[0,0,0,0,1,1]);
Append(~params3,[1,0,0,0,1,1]);
Append(~params3,[0,0,0,0,1,w]);
Append(~params3,[1,0,0,0,1,w]);
Append(~params3,[0,0,1,0,0,1]);
Append(~params3,[0,0,w,0,0,w]);
for i in [1..p-1] do
  Append(~params3,[0,0,0,i,1,1]);
  Append(~params3,[0,0,0,i,1,w]);
end for;
Append(~params3,[0,0,0,1,1,0]);
Append(~params3,[0,0,0,w,1,0]);
if p mod 4 eq 1 then
  Append(~params3,[0,0,0,w^2 mod p,1,0]);
  Append(~params3,[0,0,0,w^3 mod p,1,0]);
end if;

print #params3,2*p+8+Gcd(p-1,4);
gtotal:=gtotal+#params3;

/*
Descendants of 5.14.

Case 4: bab=bac=cac=0, cb=baa
*/

Append(~params4,[0,0,0,0,0,0]);
Append(~params4,[1,0,0,0,0,0]);
Append(~params4,[0,1,0,0,0,0]);
if p mod 3 eq 1 then
  Append(~params4,[0,w,0,0,0,0]);
  Append(~params4,[0,w^2,0,0,0,0]);
end if;

Append(~params4,[0,0,0,0,1,0]);
Append(~params4,[0,1,0,0,1,0]);
if p mod 3 eq 1 then
  Append(~params4,[0,w,0,0,1,0]);
  Append(~params4,[0,w^2,0,0,1,0]);
end if;

Append(~params4,[0,0,0,0,0,1]);
Append(~params4,[1,0,0,0,0,1]);
Append(~params4,[0,0,0,0,0,w]);
Append(~params4,[1,0,0,0,0,w]);

for y in [0..p-1] do
  Append(~params4,[0,0,1,0,0,y]);
  Append(~params4,[0,0,w,0,0,y]);
end for;
Append(~params4,[0,0,1,0,1,1]);
Append(~params4,[0,0,w,0,1,w]);
for z in [1..(p-1) div 2] do
  Append(~params4,[0,z,1,0,0,0]);
  Append(~params4,[0,z,w,0,0,0]);
end for;

Append(~params4,[0,0,0,1,0,0]);
Append(~params4,[1,0,0,1,0,0]);
if p mod 5 eq 1 then
  Append(~params4,[w,0,0,1,0,0]);
  Append(~params4,[w^2,0,0,1,0,0]);
  Append(~params4,[w^3,0,0,1,0,0]);
  Append(~params4,[w^4,0,0,1,0,0]);
end if;

Append(~params4,[0,0,0,1,1,0]);
Append(~params4,[0,0,0,1,w,0]);
if p mod 4 eq 1 then
  Append(~params4,[0,0,0,1,w^2,0]);
  Append(~params4,[0,0,0,1,w^3,0]);
end if;
for x in [0..p-1] do
  Append(~params4,[0,0,0,1,x,1]);
  Append(~params4,[0,0,0,1,x,w]);
end for;
for z in [1..(p-1) div 2] do
  Append(~params4,[z,0,0,1,0,1]);
  Append(~params4,[z,0,0,1,0,w]);
end for;

print #params4,6*p+8+2*Gcd(p-1,3)+Gcd(p-1,4)+Gcd(p-1,5);
gtotal:=gtotal+#params4;

/*
Descendants of 5.14.

Case 5: cb=bac=cac=0, caa=bab
*/

Append(~params5,[0,0,0,0,0,0]);
Append(~params5,[1,0,0,0,0,0]);
Append(~params5,[0,1,0,0,0,0]);
Append(~params5,[0,w,0,0,0,0]);

Append(~params5,[0,0,1,0,0,0]);
Append(~params5,[0,0,w,0,0,0]);
Append(~params5,[0,1,1,0,0,0]);
Append(~params5,[0,1,w,0,0,0]);
Append(~params5,[0,w,1,0,0,0]);
Append(~params5,[0,w,w,0,0,0]);

for x in [0..p-1] do
  Append(~params5,[x,0,0,1,0,0]);
end for;

Append(~params5,[0,0,0,0,1,0]);
Append(~params5,[0,1,0,0,1,0]);
Append(~params5,[0,w,0,0,1,0]);
if p mod 3 eq 1 then
  Append(~params5,[0,w^2,0,0,1,0]);
  Append(~params5,[0,w^3,0,0,1,0]);
  Append(~params5,[0,w^4,0,0,1,0]);
  Append(~params5,[0,w^5,0,0,1,0]);
end if;

if p mod 4 eq 1 then
for x in [0..(p-1) div 2] do
  Append(~params5,[0,x,0,1,1,0]);
  Append(~params5,[0,x,0,w,1,0]);
  Append(~params5,[0,x,0,w^2,1,0]);
  Append(~params5,[0,x,0,w^3,1,0]);
end for;
end if;

if p mod 4 eq 3 then
for x in [0..p-1] do
  Append(~params5,[0,x,0,1,1,0]);
  Append(~params5,[0,x,0,w,1,0]);
end for;
end if;

for x in [0..p-1] do
Append(~params5,[0,0,x,0,0,1]);
Append(~params5,[0,0,x,0,0,w]);
end for;
Append(~params5,[1,0,0,0,0,1]);
Append(~params5,[1,0,0,0,0,w]);
Append(~params5,[1,0,1,0,0,1]);
Append(~params5,[1,0,w,0,0,w]);

print #params5,5*p+13+2*Gcd(p-1,3)+Gcd(p-1,4);
gtotal:=gtotal+#params5;

/*
Descendants of 5.14.

Case 6: bac=cac=0, caa=bab, cb=baa
*/

Append(~params6,[0,0,0,0,0,0]);
Append(~params6,[1,0,0,0,0,0]);
if p mod 5 eq 1 then
  Append(~params6,[w,0,0,0,0,0]);
  Append(~params6,[w^2,0,0,0,0,0]);
  Append(~params6,[w^3,0,0,0,0,0]);
  Append(~params6,[w^4,0,0,0,0,0]);
end if;
Append(~params6,[0,1,0,0,0,0]);
Append(~params6,[0,w,0,0,0,0]);
if p mod 3 eq 1 then
  Append(~params6,[0,w^2,0,0,0,0]);
  Append(~params6,[0,w^3,0,0,0,0]);
  Append(~params6,[0,w^4,0,0,0,0]);
  Append(~params6,[0,w^5,0,0,0,0]);
end if;
if p mod 4 eq 3 then
  for x in [0..p-1] do
    Append(~params6,[0,x,1,0,0,0]);
    Append(~params6,[0,x,w,0,0,0]);
  end for;
end if;
if p mod 4 eq 1 then
  for x in [0..(p-1) div 2] do
    Append(~params6,[0,x,1,0,0,0]);
    Append(~params6,[0,x,w,0,0,0]);
    Append(~params6,[0,x,w^2,0,0,0]);
    Append(~params6,[0,x,w^3,0,0,0]);
  end for;
end if;
for x in [0..p-1] do
  Append(~params6,[x,0,0,1,0,0]);
  if p mod 5 eq 1 then
    Append(~params6,[x,0,0,w,0,0]);
    Append(~params6,[x,0,0,w^2,0,0]);
    Append(~params6,[x,0,0,w^3,0,0]);
    Append(~params6,[x,0,0,w^4,0,0]);
  end if;
end for;
if p mod 3 eq 2 then
  for x in [0..p-1] do
  for y in [0..p-1] do
    Append(~params6,[0,y,0,x,1,0]);
  end for;
  end for;
end if;
if p mod 3 eq 1 then
  xrange:=transversal3;
  Append(~xrange,0);
  for x in xrange do
  for y in [0..p-1] do
    Append(~params6,[0,y,0,x,1,0]);
    Append(~params6,[0,y,0,x,w,0]);
    Append(~params6,[0,y,0,x,w^2,0]);
  end for;
  end for;
end if;
if p mod 4 eq 3 then
  for y in [0..(p-1) div 2] do
    Append(~params6,[y,0,0,0,0,1]);
    Append(~params6,[y,0,0,0,0,w]);
  end for;
  for x in [1..p-1] do
    Append(~params6,[0,0,x,0,0,1]);
    Append(~params6,[0,0,x,0,0,w]);
  end for;
  for y in [1..(p-1) div 2] do
    Append(~params6,[y,0,1,0,0,1]);
    Append(~params6,[y,0,w,0,0,w]);
  end for;
end if;
if p mod 4 eq 1 then
  yrange:=transversal4;
  Append(~yrange,0);
  for y in yrange do
    Append(~params6,[y,0,0,0,0,1]);
    Append(~params6,[y,0,0,0,0,w]);
    Append(~params6,[y,0,0,0,0,w^2]);
    Append(~params6,[y,0,0,0,0,w^3]);
  end for;
  for x in [1..p-1] do
    Append(~params6,[0,0,x,0,0,1]);
    Append(~params6,[0,0,x,0,0,w]);
    Append(~params6,[0,0,x,0,0,w^2]);
    Append(~params6,[0,0,x,0,0,w^3]);
  end for;
  for y in transversal4 do
    Append(~params6,[y,0,1,0,0,1]);
    Append(~params6,[y,0,w,0,0,w]);
    Append(~params6,[y,0,w^2,0,0,w^2]);
    Append(~params6,[y,0,w^3,0,0,w^3]);
  end for;
end if;


print #params6,
p^2+3*p-3+(p+2)*Gcd(p-1,3)+(p+1)*Gcd(p-1,4)+(p+1)*Gcd(p-1,5);
gtotal:=gtotal+#params6;


/*
Descendants of 5.14.

Case 7: cb=baa=bac=cac=0
*/

for u in [0,1,w] do
for t in [0,1] do
for x in [0,1] do
  Append(~params7,[u,0,t,x,0,0]);
end for;
end for;
end for;

for u in [0,1,w] do
for x in [0,1] do
  Append(~params7,[u,1,0,x,0,0]);
end for;
end for;

for x in [0,1,w] do
  Append(~params7,[0,1,1,x,0,0]);
end for;

for u in [1,w] do
for x in [0..p-1] do
  Append(~params7,[u,1,1,x,0,0]);
end for;
end for;

for u in [0,1,w] do
for x in [0,1] do
for z in [1,w] do
  Append(~params7,[u,0,0,x,0,z]);
end for;
end for;
end for;

for u in [0..p-1] do
for x in [0,1] do
for z in [1,w] do
  Append(~params7,[u,0,1,x,0,z]);
end for;
end for;
end for;

for v in [0,1,w] do
for z in [0,1,w] do
  Append(~params7,[0,v,0,0,1,z]);
end for;
end for;

vrange:=[0,1,w];
if p mod 4 eq 1 then vrange:=[0,1,w,w^2,w^3]; end if;
for v in vrange do
  Append(~params7,[0,v,0,1,1,0]);
end for;

for v in [0..p-1] do
for z in [1,w] do
  Append(~params7,[0,v,0,1,1,z]);
end for;
end for;

for v in [0..p-1] do
for z in [0,1,w] do
  Append(~params7,[0,v,1,0,1,z]);
end for;
end for;

for v in [0..p-1] do
for x in [1,w] do
for z in [0..p-1] do
  Append(~params7,[0,v,1,x,1,z]);
end for;
end for;
end for;

print #params7,2*p^2+11*p+43+Gcd(p-1,4);
gtotal:=gtotal+#params7;

/*
Descendants of 5.14.

Case 8: cb=caa, baa=bac=cac=0
*/


urange:=[0,1,w];
if p mod 4 eq 1 then urange:=[0,1,w,w^2,w^3]; end if;
for u in urange do
  Append(~params8,[u,0,0,0,0,0]);
  Append(~params8,[u,0,0,1,0,0]);
  Append(~params8,[u,1,0,0,0,0]);
end for;

if p mod 3 eq 1 then
  urange:=transversal3;
  Append(~urange,0);
else;
  urange:=[0..p-1];
end if;
for u in urange do
  Append(~params8,[u,0,1,0,0,0]);
  Append(~params8,[u,0,1,1,0,0]);
  Append(~params8,[u,1,1,0,0,0]);
  if p mod 3 eq 1 then
    Append(~params8,[u,0,w,0,0,0]);
    Append(~params8,[u,0,w,1,0,0]);
    Append(~params8,[u,1,w,0,0,0]);
    Append(~params8,[u,0,w^2,0,0,0]);
    Append(~params8,[u,0,w^2,1,0,0]);
    Append(~params8,[u,1,w^2,0,0,0]);
  end if;
end for;

for u in [0..p-1] do
for t in [0..p-1] do
  Append(~params8,[u,1,t,1,0,0]);
end for;
end for;

for u in [0..p-1] do
for t in [0..(p-1) div 2] do
for x in [0,1] do
for z in [1,w] do
  Append(~params8,[u,0,t,x,0,z]);
end for;
end for;
end for;
end for;

vrange:=[0,1,w];
if p mod 3 eq 1 then vrange:=[0,1,w,w^2,w^3,w^4,w^5]; end if;
for v in vrange do
  Append(~params8,[0,v,0,0,1,0]);
end for;

if p mod 5 ne 1 then
  for v in [0..p-1] do
    Append(~params8,[0,v,0,1,1,0]);
  end for;
end if;

if p mod 5 eq 1 then
  vrange:=transversal5;
  Append(~vrange,0);
  for v in vrange do
  for x in [1,w,w^2,w^3,w^4] do
    Append(~params8,[0,v,0,x,1,0]);
  end for;
  end for;
end if;

if p mod 3 eq 2 then
  for v in [0..p-1] do
  for x in [0..p-1] do
    Append(~params8,[0,v,1,x,1,0]);
  end for;
  end for;
end if;

if p mod 3 eq 1 then
  xrange:=transversal3;
  Append(~xrange,0);
  for v in [0..p-1] do
  for x in xrange do
    Append(~params8,[0,v,1,x,1,0]);
    Append(~params8,[0,v,w,x,1,0]);
    Append(~params8,[0,v,w^2,x,1,0]);
  end for;
  end for;
end if;

for v in [0..p-1] do
for x in [0..(p-1) div 2] do
for z in [1,w] do
  Append(~params8,[0,v,0,x,1,z]);
end for;
end for;
end for;

for v in [0..p-1] do
for t in [1..(p-1) div 2] do
for x in [0..p-1] do
for z in [1,w] do
  Append(~params8,[0,v,t,x,1,z]);
end for;
end for;
end for;
end for;

print #params8,
p^3+4*p^2+6*p+(p+5)*Gcd(p-1,3)+3*Gcd(p-1,4)+Gcd(p-1,5);
gtotal:=gtotal+#params8;

/*
Descendants of 5.14.

Cases 9 & 10: cb=bac=caa=0, cac=kbab (k=1,w)
*/

Append(~params9,[0,0,0,0,0,0]);
Append(~params9,[1,0,0,0,0,0]);
Append(~params9,[0,1,0,0,0,0]);
Append(~params9,[0,w,0,0,0,0]);

for y in [0,1,w] do
  Append(~params9,[0,y,0,0,1,0]);
  if p mod 4 eq 1 then
    Append(~params9,[0,y,0,0,w,0]);
   end if;
end for;
yrange:=[0..p-1];
if p mod 4 eq 1 then yrange:=[0..(p-1) div 2]; end if;
for y in yrange do
  Append(~params9,[1,y,0,0,1,0]);
  if p mod 4 eq 1 then
    Append(~params9,[1,y,0,0,w,0]);
  end if;
end for;

for x in [0..(p-1) div 2] do
  Append(~params9,[x,0,0,0,0,1]);
end for;
Append(~params9,[0,1,0,0,0,1]);
Append(~params9,[0,w,0,0,0,1]);

for y in [0,1,w] do
for z in [0..(p-1) div 2] do
  Append(~params9,[0,y,1,0,z,0]);
  Append(~params9,[0,y,w,0,z,0]);
end for;
end for;

for y in [0..p-1] do
for z in [0..(p-1) div 2] do
  Append(~params9,[0,y,1,0,z,1]);
  Append(~params9,[0,y,w,0,z,1]);
end for;
end for;

for y in [0..p-1] do
for t in [0..(p-1) div 2] do
  Append(~params9,[1,y,1,0,0,t]);
  Append(~params9,[1,y,w,0,0,t]);
end for;
end for;

for y in [0..p-1] do
for z in [1..(p-1) div 2] do
for t in [0..p-1] do
  Append(~params9,[1,y,1,0,z,t]);
  Append(~params9,[1,y,w,0,z,t]);
end for;
end for;
end for;

for x in [0..(p-1) div 2] do
for y in [0..p-1] do
  Append(~params9,[y,0,0,1,0,x]);
end for;
end for;
for x in [0..(p-1) div 2] do
  Append(~params9,[1,1,0,1,0,x]);
  Append(~params9,[1,w,0,1,0,x]);
end for;

if p mod 4 eq 1 then
  for x in [0..p-1] do
  for y in [0..(p-1) div 2] do
    Append(~params9,[x,y,0,1,1,0]);
    Append(~params9,[x,y,0,1,w,0]);
  end for;
  end for;
end if;

if p mod 4 eq 3 then
  for x in [0..p-1] do
  for y in [0..p-1] do
    Append(~params9,[x,y,0,1,1,0]);
  end for;
  end for;
end if;

print #params9,
p^3+5*p^2/2+7*p+19/2+(p+4)*Gcd(p-1,4)/2;
gtotal:=gtotal+#params9;

params10:=params9;
gtotal:=gtotal+#params9;

/*
Descendants of 5.14.

Cases 11 & 12: bac=caa=0, cb=baa, cac=kbab (k=1,w)
*/

Append(~params11,[0,0,0,0,0,0]);
Append(~params11,[1,0,0,0,0,0]);
if p mod 3 eq 1 then
  Append(~params11,[w,0,0,0,0,0]);
  Append(~params11,[w^2,0,0,0,0,0]);
end if;
Append(~params11,[0,1,0,0,0,0]);
Append(~params11,[0,w,0,0,0,0]);
if p mod 4 eq 1 then
  Append(~params11,[0,w^2,0,0,0,0]);
  Append(~params11,[0,w^3,0,0,0,0]);
end if;

if p mod 4 eq 1 then
  xrange:=transversal4;
  Append(~xrange,0);
  for x in xrange do
  for y in [0..p-1] do
    Append(~params11,[x,y,0,0,1,0]);
    Append(~params11,[x,y,0,0,w,0]);
  end for;
  end for;
end if;
if p mod 4 eq 3 then
  for x in [0..(p-1) div 2] do
  for y in [0..p-1] do
    Append(~params11,[x,y,0,0,1,0]);
  end for;
  end for;
end if;

if p mod 3 eq 1 then
  for x in [0..(p-1) div 2] do
    Append(~params11,[x,0,0,0,0,1]);
    Append(~params11,[x,0,0,0,0,w]);
    Append(~params11,[x,0,0,0,0,w^2]);
  end for;
  for x in transversal3 do
    Append(~params11,[0,x,0,0,0,1]);
    Append(~params11,[0,x,0,0,0,w]);
    Append(~params11,[0,x,0,0,0,w^2]);
  end for;
end if;
if p mod 3 eq 2 then
  for x in [0..(p-1) div 2] do
    Append(~params11,[x,0,0,0,0,1]);
  end for;
  for x in [1..p-1] do
    Append(~params11,[0,x,0,0,0,1]);
  end for;
end if;

for y in [0..p-1] do
for z in [0..(p-1) div 2] do
for t in [0..(p-1) div 2] do
  Append(~params11,[0,y,1,0,z,t]);
  Append(~params11,[0,y,w,0,z,t]);
end for;
end for;
end for;
for y in [0..p-1] do
for x in [1..(p-1) div 2] do
for z in [0..(p-1) div 2] do
  Append(~params11,[x,y,1,0,z,0]);
  Append(~params11,[x,y,w,0,z,0]);
end for;
end for;
end for;
for y in [0..p-1] do
for x in [1..(p-1) div 2] do
for t in [1..(p-1) div 2] do
for z in [0..p-1] do
  Append(~params11,[x,y,1,0,z,t]);
  Append(~params11,[x,y,w,0,z,t]);
end for;
end for;
end for;
end for;


if p mod 3 eq 1 then
  for x in [0..p-1] do
    Append(~params11,[x,0,0,1,0,0]);
    Append(~params11,[x,0,0,w,0,0]);
    Append(~params11,[x,0,0,w^2,0,0]);
  end for;
  for y in transversal3 do
    Append(~params11,[1,y,0,1,0,0]);
    Append(~params11,[w,y,0,w,0,0]);
    Append(~params11,[w^2,y,0,w^2,0,0]);
  end for;
  for x in transversal6 do
  for y in [0..p-1] do
  for z in [0..p-1] do
    Append(~params11,[y,z,0,1,x,0]);
    Append(~params11,[y,z,0,w,x,0]);
    Append(~params11,[y,z,0,w^2,x,0]);
  end for;
  end for;
  end for;
  for x in [0..p-1] do
  for z in [1..(p-1) div 2] do
    Append(~params11,[x,0,0,1,0,z]);
    Append(~params11,[x,0,0,w,0,z]);
    Append(~params11,[x,0,0,w^2,0,z]);
  end for;
  end for;
  for y in transversal3 do
  for z in [1..(p-1) div 2] do
    Append(~params11,[1,y,0,1,0,z]);
    Append(~params11,[w,y,0,w,0,z]);
    Append(~params11,[w^2,y,0,w^2,0,z]);
  end for;
  end for;
end if;

if p mod 3 eq 2 then
  for x in [0..p-1] do
    Append(~params11,[x,0,0,1,0,0]);
  end for;
  for y in [1..p-1] do
    Append(~params11,[1,y,0,1,0,0]);
  end for;
  for x in [1..(p-1) div 2] do
  for y in [0..p-1] do
  for z in [0..p-1] do
    Append(~params11,[y,z,0,1,x,0]);
  end for;
  end for;
  end for;
  for x in [0..p-1] do
  for z in [1..(p-1) div 2] do
    Append(~params11,[x,0,0,1,0,z]);
  end for;
  end for;
  for y in [1..p-1] do
  for z in [1..(p-1) div 2] do
    Append(~params11,[1,y,0,1,0,z]);
  end for;
  end for;
end if;


print #params11,
(p^4+p^3+4*p^2+p-1+(p^2+2*p+3)*Gcd(p-1,3)+(p+2)*Gcd(p-1,4))/2;
gtotal:=gtotal+#params11;


params12:=params11;
gtotal:=gtotal+#params11;

/*
Descendants of 5.14.

Case 13, cb=bac=0, caa=baa, cac=-bab
*/

//pb=pc=0
Append(~params13,[0,0,0,0,0,0]);
Append(~params13,[1,0,0,0,0,0]);
Append(~params13,[0,1,0,0,0,0]);

//pb=0, pc=baa or wbaa
Append(~params13,[0,0,0,0,1,0]);
Append(~params13,[0,0,0,0,w,0]);
Append(~params13,[0,1,0,0,1,0]);
Append(~params13,[0,1,0,0,w,0]);
Append(~params13,[0,w,0,0,1,0]);
Append(~params13,[0,w,0,0,w,0]);
for y in [0..p-1] do
  Append(~params13,[1,y,0,0,1,0]);
  Append(~params13,[1,y,0,0,w,0]);
end for;

//pb=0, pc=bab
for x in [0..p-1] do
  Append(~params13,[x,0,0,0,0,1]);
end for;
Append(~params13,[-1,1,0,0,0,1]);
Append(~params13,[-1,w,0,0,0,1]);

//pb=pc=bab
Append(~params13,[0,0,0,1,0,1]);
Append(~params13,[1,0,0,1,0,1]);
Append(~params13,[0,1,0,1,0,1]);

//pb=bab, pc=-bab
for x in [0..p-1] do
  Append(~params13,[x,0,0,1,0,-1]);
end for;
Append(~params13,[2,1,0,1,0,-1]);

//pb=bab, pc=baa or wbaa
for x in [0..p-1] do
for y in [0..p-1] do
  Append(~params13,[x,y,0,1,1,0]);
  Append(~params13,[x,y,0,1,w,0]);
end for;
end for;

//pb=pc=baa
Append(~params13,[0,0,1,0,1,0]);
Append(~params13,[1,0,1,0,1,0]);
Append(~params13,[0,1,1,0,1,0]);
Append(~params13,[1,1,1,0,1,0]);

//pb=baa, pc=baa+bab
for x in [0..p-2] do //Note that -1 is excluded
  Append(~params13,[x,0,1,0,1,1]);
  y:=F!((x+1)/2); y:=Z!y;
  Append(~params13,[x,y,1,0,1,1]);
end for;
Append(~params13,[-1,0,1,0,1,1]);
Append(~params13,[-1,1,1,0,1,1]);

//pb=baa, pc=-baa
Append(~params13,[0,0,1,0,-1,0]);
Append(~params13,[0,1,1,0,-1,0]);
Append(~params13,[0,w,1,0,-1,0]);
for y in [0..p-1] do
  Append(~params13,[1,y,1,0,-1,0]);
end for;

//pb=baa, pc=-baa+bab
half:=(p+1) div 2;
Append(~params13,[1,-half,1,0,-1,1]);
Append(~params13,[1,half,1,0,-1,1]);
Append(~params13,[1,w-half,1,0,-1,1]);
if p mod 4 eq 1 then
  Append(~params13,[1,w^2-half,1,0,-1,1]);
  Append(~params13,[1,w^3-half,1,0,-1,1]);
end if;
for y in [0..p-1] do
  Append(~params13,[0,y,1,0,-1,1]);
  Append(~params13,[1-w,y,1,0,-1,1]);
end for;

//pb=pc=wbaaa
Append(~params13,[0,0,w,0,w,0]);
Append(~params13,[1,0,w,0,w,0]);
Append(~params13,[0,1,w,0,w,0]);
Append(~params13,[1,1,w,0,w,0]);

//pb=wbaa, pc=wbaa+bab
for x in [0..p-2] do //Note that -1 is excluded
  Append(~params13,[x,0,w,0,w,1]);
  y:=F!((x+1)/(2*w)); y:=Z!y;
  Append(~params13,[x,y,w,0,w,1]);
end for;
Append(~params13,[-1,0,w,0,w,1]);
Append(~params13,[-1,1,w,0,w,1]);

print #params13,2*p^2+11*p+27+Gcd(p-1,4);
gtotal:=gtotal+#params13;

/*
Descendants of 5.14.

Case 14, bac=0, cb=caa=baa, cac=-bab
*/

//1. v ne +/- y 
for y in [1,w] do
for u in [0..p-1] do
for x in [0..(p-1) div 2] do
  zrange:=[0..p-1];
  if x eq 0 then zrange:=[0..(p-1) div 2]; end if;
  for z in zrange do
    Append(~params14,[0,u,0,x,y,z]);
  end for;
end for;
end for;
end for;

//1.5 t=0, v=1, y=-1, u ne 0
for u in [1,w] do
for x in [0..(p-1) div 2] do
  zrange:=[0..p-1];
  if x eq 0 then zrange:=[0..(p-1) div 2]; end if;
  for z in zrange do
    Append(~params14,[0,u,1,x,-1,z]);
  end for;
end for;
end for;

//1.5 t=0, v=1, y=-1, u=0
for x in [0..p-1] do
for z in [0..p-1] do
  new:=true;
  for a in [1..p-1] do
    x1:=F!((x+z+x*a^2-z*a^2)/(2*a^3)); x1:=Z!x1;
    z1:=F!((x+z-x*a^2+z*a^2)/(2*a^3)); z1:=Z!z1;
    if [x1,z1] lt [x,z] then new:=false; break; end if;
  end for;
  if new then Append(~params14,[0,0,1,x,-1,z]); end if;
end for;
end for;



//2. t ne 0, v=y, x ne z+1, v ne 0
for v in [1,w] do
for x in [0..p-1] do
for z in [0..p-1] do
  if x eq (z+1) mod p then continue; end if;
  x1:=(-x) mod p; z1:=(-z) mod p;
  if [z1,x1] lt [x,z] then continue; end if;
  Append(~params14,[1,0,v,x,v,z]);
end for;
end for;
end for;

//2. t ne 0, v=y=0, x ne z+1
for x in [0..p-1] do
for z in [0..p-1] do
  if x eq (z+1) mod p then continue; end if;
  new:=true;
  for a in [1..p-1] do
    x1:=F!((x+z+x*a^3-z*a^3)/(2*a^3)); x1:=Z!x1;
    z1:=F!((x+z-x*a^3+z*a^3)/(2*a^3)); z1:=Z!z1;
    if [x1,z1] lt [x,z] then new:=false; break; end if;
  end for;
  if new then Append(~params14,[1,0,0,x,0,z]); end if;
end for;
end for;

//3. t ne 0, v=y, x eq z+1, v ne 0
for v in [1,w] do
for u in [0..(p-1) div 2] do
  Append(~params14,[1,u,v,0,v,-1]);
end for;
end for;

//3. t ne 0, v=y, x eq z+1, v eq 0, u ne 0
for z in [0..p-1] do
  Append(~params14,[1,1,0,z+1,0,z]);
end for;

//3. t ne 0, v=y, x eq z+1, v eq 0, u eq 0
for z in [0..p-1] do
  new:=true;
  for a in [1..p-1] do
    z1:=F!((-a^3+2*z+1)/(2*a^3)); z1:=Z!z1;
    if z1 lt z then new:=false; break; end if;
  end for;
  if new then Append(~params14,[1,0,0,z+1,0,z]); end if;
end for;

//4. t=0, v=y, x ne +/- z, v ne 0
for v in [1,w] do
for z in [1..(p-1) div 2] do
  Append(~params14,[0,0,v,0,v,z]);
end for;
end for;

//4. t=0, v=y=0, x ne +/- z
Append(~params14,[0,0,0,0,0,1]);
if p mod 3 eq 1 then
  Append(~params14,[0,0,0,0,0,w]);
  Append(~params14,[0,0,0,0,0,w^2]);
end if;

//4.5 t=0, v=y, x=-z ne 0
Append(~params14,[0,0,0,1,0,-1]);
Append(~params14,[0,0,w,1,w,-1]);
Append(~params14,[0,0,w^2,1,w^2,-1]);

//5. t=0, v=y ne 0, x=z
for v in [1,w] do
for u in [0,1] do
  Append(~params14,[0,u,v,0,v,0]);
end for;
end for;


//6. t=0, v=y=0, x=z
xrange:=[0,1];
if p mod 3 eq 1 then xrange:=[0,1,w,w^2]; end if;
for u in [0,1] do
for x in xrange do
  Append(~params14,[0,u,0,x,0,x]);
end for;
end for;

print #params14,p^3+2*p^2+6*p+10+(p+4)*Gcd(p-1,3);
gtotal:=gtotal+#params14;

/*
Descendants of 5.14.

Case 15, cb=baa=bac=caa=0
*/

//x=y=z=v=0
Append(~params15,[0,0,0,0,0,0]);
Append(~params15,[0,0,0,0,1,0]);
Append(~params15,[0,0,0,1,1,0]);
if p mod 3 eq 1 then
  Append(~params15,[0,0,0,1,w,0]);
end if;

//x=z=v=0, y=1,w
Append(~params15,[0,1,0,0,0,0]);
Append(~params15,[0,w,0,0,0,0]);
Append(~params15,[0,1,0,0,1,0]);
Append(~params15,[0,w,0,0,1,0]);
Append(~params15,[0,1,0,1,0,0]);
Append(~params15,[0,w,0,1,0,0]);
Append(~params15,[0,1,0,1,1,0]);
Append(~params15,[0,w,0,1,1,0]);
if  p mod 3 eq 1 then
  Append(~params15,[0,1,0,1,w,0]);
  Append(~params15,[0,w,0,1,w,0]);
  Append(~params15,[0,1,0,1,w^2,0]);
  Append(~params15,[0,w,0,1,w^2,0]);
end if;

//z=v=0, x=y=1,w
Append(~params15,[1,1,0,0,0,0]);
Append(~params15,[w,w,0,0,0,0]);
Append(~params15,[1,1,0,0,1,0]);
Append(~params15,[w,w,0,0,1,0]);
for u in [1..(p-1) div 2] do
  u1:=F!u^-1; u1:=Z!u1;
  if u1 lt u or (p-u1) lt u then continue; end if;
  Append(~params15,[1,1,0,1,u,0]);
  Append(~params15,[w,w,0,1,u,0]);
end for;

//z=v=0, x=1, y=w
Append(~params15,[1,w,0,0,0,0]);
Append(~params15,[1,w,0,0,1,0]);
Append(~params15,[1,w,0,1,0,0]);
for u in [1..(p-1) div 2] do
  Append(~params15,[1,w,0,1,u,0]);
end for;

//z=0, v=1, t=0
Append(~params15,[0,0,0,0,0,1]);
Append(~params15,[0,0,0,0,1,1]);
Append(~params15,[0,0,0,0,w,1]);
Append(~params15,[0,1,0,0,0,1]);
Append(~params15,[0,1,0,0,1,1]);
Append(~params15,[0,1,0,0,w,1]);
Append(~params15,[0,w,0,0,0,1]);
Append(~params15,[0,w,0,0,1,1]);
Append(~params15,[0,w,0,0,w,1]);
Append(~params15,[1,0,0,0,0,1]);
Append(~params15,[1,0,0,0,1,1]);
Append(~params15,[1,0,0,0,w,1]);
Append(~params15,[w,0,0,0,0,1]);
Append(~params15,[w,0,0,0,1,1]);
Append(~params15,[w,0,0,0,w,1]);
for x in [1,w] do
for y in [1,w] do
for u in [0..p-1] do
  Append(~params15,[x,y,0,0,u,1]);
end for;
end for;
end for;

//z=0, v=1, t=1
for x in [1,w] do
for y in [0..p-1] do
for u in [0..p-1] do
  Append(~params15,[x,y,0,1,u,1]);
end for;
end for;
end for;
for y in [0,1,w] do
for u in [0..p-1] do
  Append(~params15,[0,y,0,1,u,1]);
end for;
end for;

//z=v=1
for t in [0..p-2] do
for u in [t+1..p-1] do
  Append(~params15,[0,0,1,t,u,1]);
  Append(~params15,[0,1,1,t,u,1]);
  Append(~params15,[0,w,1,t,u,1]);
  for y in [0..p-1] do
    Append(~params15,[1,y,1,t,u,1]);
    Append(~params15,[w,y,1,t,u,1]);
  end for;
end for;
end for;

//Get set of non-zero squares with k^2 ~ k^-2
range:={};
for k in [1..(p-1) div 2] do
  k2:=k^2 mod p;
  km2:=F!(k2^-1); km2:=Z!km2;
  if k2 le km2 then Include(~range,k2); end if;
end for;
for t in [0..p-1] do
   Append(~params15,[0,0,1,t,t,1]);
   Append(~params15,[0,1,1,t,t,1]);
   Append(~params15,[0,w,1,t,t,1]);
   for k2 in range do
      Append(~params15,[1,k2,1,t,t,1]);
      Append(~params15,[w,w*k2,1,t,t,1]);
   end for;
   for k in [1..(p-1) div 2] do
      Append(~params15,[1,w*k^2,1,t,t,1]);
   end for;   
end for;

print #params15,
p^3+(7*p^2+17*p+59+5*Gcd(p-1,3)+(p+1)*Gcd(p-1,4))/2;
gtotal:=gtotal+#params15;

/*
Descendants of 5.14.

Case 16, cb=bac=caa=0, baa=cac
*/

//u=v=z=0
for y in [0,1,w] do
  Append(~params16,[0,y,0,0,0,0]);
  Append(~params16,[0,y,0,1,0,0]);
  Append(~params16,[0,y,0,w,0,0]);
  for t in [0..p-1] do
    Append(~params16,[1,y,0,t,0,0]);
    Append(~params16,[w,y,0,t,0,0]);
  end for;
end for;

//u=v=0, z=1,w
Append(~params16,[0,0,1,0,0,0]);
Append(~params16,[0,0,w,0,0,0]);
Append(~params16,[0,1,1,0,0,0]);
Append(~params16,[0,1,w,0,0,0]);
Append(~params16,[0,w,1,0,0,0]);
Append(~params16,[0,w,w,0,0,0]);
if p mod 4 eq 1 then
  Append(~params16,[0,w^2,1,0,0,0]);
  Append(~params16,[0,w^2,w,0,0,0]);
  Append(~params16,[0,w^3,1,0,0,0]);
  Append(~params16,[0,w^3,w,0,0,0]);
end if;
for y in [0..p-1] do
  Append(~params16,[0,y,1,1,0,0]);
  Append(~params16,[0,y,w,1,0,0]);
  Append(~params16,[0,y,1,w,0,0]);
  Append(~params16,[0,y,w,w,0,0]);
end for;
for y in [0..p-1] do
for t in [0..p-1] do
  Append(~params16,[1,y,1,t,0,0]);
  Append(~params16,[1,y,w,t,0,0]);
  Append(~params16,[w,y,1,t,0,0]);
  Append(~params16,[w,y,w,t,0,0]);
end for;
end for;

//u=0, v=1
Append(~params16,[0,0,0,0,0,1]);
Append(~params16,[0,0,1,0,0,1]);
Append(~params16,[0,0,w,0,0,1]);
if p mod 3 eq 1 then
  Append(~params16,[0,0,w^2,0,0,1]);
  Append(~params16,[0,0,w^3,0,0,1]);
  Append(~params16,[0,0,w^4,0,0,1]);
  Append(~params16,[0,0,w^5,0,0,1]);
end if;
for z in [0..p-1] do
  Append(~params16,[0,0,z,1,0,1]);
  Append(~params16,[0,0,z,w,0,1]);
end for;
for z in [0..p-1] do
for t in [0..p-1] do
  Append(~params16,[0,1,z,t,0,1]);
  Append(~params16,[0,w,z,t,0,1]);
end for;
end for;
for y in [0..p-1] do
for z in [0..p-1] do
for t in [0..p-1] do
  Append(~params16,[1,y,z,t,0,1]);
  Append(~params16,[w,y,z,t,0,1]);
end for;
end for;
end for;

//u=1
Append(~params16,[0,0,0,0,1,0]);
Append(~params16,[0,1,0,0,1,0]);
Append(~params16,[0,w,0,0,1,0]);
if p mod 3 eq 1 then
  Append(~params16,[0,w^2,0,0,1,0]);
  Append(~params16,[0,w^3,0,0,1,0]);
  Append(~params16,[0,w^4,0,0,1,0]);
  Append(~params16,[0,w^5,0,0,1,0]);
end if;
if p mod 4 eq 1 then
  for v in [1,w,w^2,w^3] do
  for y in [0..(p-1) div 2] do  
    Append(~params16,[0,y,0,0,1,v]);
  end for;
  end for;
end if;
if p mod 4 eq 3 then
  for v in [1,w] do
  for y in [0..p-1] do  
    Append(~params16,[0,y,0,0,1,v]);
  end for;
  end for;
end if;
for v in [0..p-1] do
for y in [0..p-1] do
  Append(~params16,[0,y,0,1,1,v]);
  Append(~params16,[0,y,0,w,1,v]);
end for;
end for;
for t in [0..p-1] do
for v in [0..p-1] do
for y in [0..p-1] do
  Append(~params16,[0,y,1,t,1,v]);
  Append(~params16,[0,y,w,t,1,v]);
end for;
end for;
end for;
for z in [0..p-1] do
for t in [0..p-1] do
for v in [0..p-1] do
for y in [0..p-1] do
  Append(~params16,[1,y,z,t,1,v]);
  Append(~params16,[w,y,z,t,1,v]);
end for;
end for;
end for;
end for;

print #params16,
2*p^4+4*p^3+8*p^2+14*p+11+4*Gcd(p-1,3)+3*Gcd(p-1,4);
gtotal:=gtotal+#params16;

/*
Descendants of 5.14.

Case 17, cb=bac=0, cac=baa, caa=bab
*/

mats:=[];
//get pc,pb

for y5 in [0,1,lns] do
for y6 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y3,y4,y5,y6];
if A eq 0 then
  Append(~mats,A);
  continue;
end if;

new:=1;
index:=p^3*y5+p^2*y6+p*y3+y4;

for a in [1..p-1] do
for x in CU do
c:=a*x;

B:=H22![a^-1*c^2,0,0,c];
C:=H22![a*c^2,0,0,a^2*c];

D:=B*A*C^-1;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);


ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=0; end if;

B:=H22![0,a^-1*c^2,c,0];
C:=H22![0,a*c^2,a^2*c,0];

D:=B*A*C^-1;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);


ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=0; end if;


if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~mats,A);
  //print y3,y4,y5,y6;
end if;

end for;
end for;
end for;
end for;


for AS in mats do
//Get stabilizer of AS
stab1:=[];
stab2:=[];
for a in [1..p-1] do
for x in CU do
c:=a*x;

B:=H22![a^-1*c^2,0,0,c];
C:=H22![a*c^2,0,0,a^2*c];

D:=B*AS*C^-1;

if D eq AS then
  Append(~stab1,H33![a,0,0,0,a^-1*c^2,0,0,0,c]);
  Append(~stab2,C);
end if;

B:=H22![0,a^-1*c^2,c,0];
C:=H22![0,a*c^2,a^2*c,0];

D:=B*AS*C^-1;

if D eq AS then
  Append(~stab1,H33![a,0,0,0,0,a^-1*c^2,0,c,0]);
  Append(~stab2,C);
end if;

end for;
end for;

y3:=Z!(AS[1][1]); y4:=Z!(AS[1][2]); y5:=Z!(AS[2][1]); y6:=Z!(AS[2][2]);

for y1 in [0..p-1] do
for y2 in [0..p-1] do

A:=H32![y1,y2,y3,y4,y5,y6];
if A eq 0 then
  Append(~params17,[0,0,0,0,0,0]);
  continue;
end if;

new:=1;
index:=p*y1+y2;

for ii in [1..#stab1] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);

ind1:=p*z1+z2;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;


if new eq 1 then
  Append(~params17,[y1,y2,y3,y4,y5,y6]);
end if;

end for;
end for;

end for;

print #params17,
(p^4+2*p^3+3*p^2+4*p+2)*(p-1)/Gcd(p-1,3)+3*p+4+(p^2+p+1)*Gcd(p-1,4)/2;
gtotal:=gtotal+#params17;

/*
Descendants of 5.14.

Case 18, cb=bac=0, cac=baa, caa=wbab (p=1mod3 only)
*/
if p mod 3 eq 1 then

mats:=[];
//get pc,pb

for y5 in [0,1,lns] do
for y6 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y3,y4,y5,y6];
if A eq 0 then
  Append(~mats,A);
  continue;
end if;

new:=1;
index:=p^3*y5+p^2*y6+p*y3+y4;

for a in [1..p-1] do
for x in CU do
c:=a*x;

B:=H22![a^-1*c^2,0,0,c];
C:=H22![a*c^2,0,0,a^2*c];

D:=B*A*C^-1;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);


ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=0; end if;


if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~mats,A);
  //print y3,y4,y5,y6;
end if;

end for;
end for;
end for;
end for;

for AS in mats do
//Get stabilizer of AS
stab1:=[];
stab2:=[];
for a in [1..p-1] do
for x in CU do
c:=a*x;

B:=H22![a^-1*c^2,0,0,c];
C:=H22![a*c^2,0,0,a^2*c];

D:=B*AS*C^-1;

if D eq AS then
  Append(~stab1,H33![a,0,0,0,a^-1*c^2,0,0,0,c]);
  Append(~stab2,C);
end if;

end for;
end for;

y3:=Z!(AS[1][1]); y4:=Z!(AS[1][2]); y5:=Z!(AS[2][1]); y6:=Z!(AS[2][2]);

for y1 in [0..p-1] do
for y2 in [0..p-1] do

A:=H32![y1,y2,y3,y4,y5,y6];
if A eq 0 then
  Append(~params18,[0,0,0,0,0,0]);
  continue;
end if;

new:=1;
index:=p*y1+y2;

for ii in [1..#stab1] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);

ind1:=p*z1+z2;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;


if new eq 1 then
  Append(~params18,[y1,y2,y3,y4,y5,y6]);
end if;

end for;
end for;

end for;

print #params18,
(2*p^5+2*p^4+2*p^3+2*p^2+14*p+17)/3;
gtotal:=gtotal+#params18;

end if;

/*
Descendants of 5.14.

Case 19, cb=baa=caa=cac=0
*/

Append(~params19,[0,0,0,0,0,0]);
Append(~params19,[1,0,0,0,0,0]);
Append(~params19,[w,0,0,0,0,0]);
Append(~params19,[0,1,0,0,0,0]);

for x in [0,1,w] do
for y in [0,1] do
Append(~params19,[x,y,1,0,0,0]);
end for;
end for;

Append(~params19,[0,0,0,1,0,0]);
Append(~params19,[1,0,0,1,0,0]);
Append(~params19,[w,0,0,1,0,0]);
Append(~params19,[0,1,0,1,0,0]);

for t in [0..p-1] do
  Append(~params19,[0,0,1,t,1,0]);
  Append(~params19,[0,1,1,t,1,0]);
  Append(~params19,[0,w,1,t,1,0]);
  for y in [0..p-1] do
    Append(~params19,[1,y,1,t,1,0]);
    Append(~params19,[w,y,1,t,1,0]);
  end for;
end for;

for x in [0,1,w] do
for t in [0,1,w] do
  Append(~params19,[x,0,0,t,1,0]);
end for;
end for;

for x in [1,w] do
for t in [0..p-1] do
  Append(~params19,[x,1,0,t,1,0]);
end for;
end for;

Append(~params19,[0,1,0,0,1,0]);
Append(~params19,[0,1,0,1,1,0]);
Append(~params19,[0,1,0,w,1,0]);
if p mod 4 eq 1 then
  Append(~params19,[0,1,0,w^2,1,0]);
  Append(~params19,[0,1,0,w^3,1,0]);
end if;

half:=(p+1) div 2;
for x in [0..p-1] do
  if x eq half then continue; end if;
  for z in [0,1,w] do
  for t in [0,1] do
    Append(~params19,[z,t,x,0,0,1]);
  end for;
  end for;
end for;

Append(~params19,[1,0,half,0,0,1]);
Append(~params19,[w,0,half,0,0,1]);
Append(~params19,[1,0,half,1,0,1]);
Append(~params19,[w,0,half,1,0,1]);

Append(~params19,[0,0,half,0,0,1]);
Append(~params19,[0,0,half,1,0,1]);
Append(~params19,[0,1,half,0,0,1]);
Append(~params19,[0,1,half,1,0,1]);
Append(~params19,[0,w,half,1,0,1]);


print #params19,2*p^2+11*p+27+Gcd(p-1,4);
gtotal:=gtotal+#params19;

/*
Descendants of 5.14.

Case 20, cb=baa=cac, caa=bab
*/

//z=v=y=0
for u in [0,1,w] do
  Append(~params20,[0,0,0,0,u,0]);
  Append(~params20,[0,0,0,1,u,0]);
  Append(~params20,[0,0,0,w,u,0]);
  for t in [0..p-1] do
    Append(~params20,[1,0,0,t,u,0]);
    Append(~params20,[w,0,0,t,u,0]);
  end for;
end for;

//z=v=0, y=1, x=t=0
Append(~params20,[0,1,0,0,0,0]);
Append(~params20,[0,1,0,0,1,0]);
Append(~params20,[0,1,0,0,w,0]);
if p mod 3 eq 1 then
  Append(~params20,[0,1,0,0,w^2,0]);
  Append(~params20,[0,1,0,0,w^3,0]);
  Append(~params20,[0,1,0,0,w^4,0]);
  Append(~params20,[0,1,0,0,w^5,0]);
end if;

//z=v=0, y=1, x=0, t=1,w 
for u in [0..p-1] do
  Append(~params20,[0,1,0,1,u,0]);
  Append(~params20,[0,1,0,w,u,0]);
end for;

//z=v=0, y=1, x=1,w 
for t in [0..p-1] do
for u in [0..p-1] do
  Append(~params20,[1,1,0,t,u,0]);
  Append(~params20,[w,1,0,t,u,0]);
end for;
end for;

//z=0, v=1, x=t=u=0
Append(~params20,[0,0,0,0,0,1]);
Append(~params20,[0,1,0,0,0,1]);
Append(~params20,[0,w,0,0,0,1]);
if p mod 4 eq 1 then
  Append(~params20,[0,w^2,0,0,0,1]);
  Append(~params20,[0,w^3,0,0,0,1]);
end if;

//z=0, v=1, x=t=0, u=1,w
for y in [0..p-1] do
  Append(~params20,[0,y,0,0,1,1]);
  Append(~params20,[0,y,0,0,w,1]);
end for;

//z=0, v=1, x=0, t=1,w
for y in [0..p-1] do
for u in [0..p-1] do
  Append(~params20,[0,y,0,1,u,1]);
  Append(~params20,[0,y,0,w,u,1]);
end for;
end for;

//z=0, v=1, x=1,w
for y in [0..p-1] do
for t in [0..p-1] do
for u in [0..p-1] do
  Append(~params20,[1,y,0,t,u,1]);
  Append(~params20,[w,y,0,t,u,1]);
end for;
end for;
end for;

//z=1
for v in [0..p-1] do

Append(~params20,[0,0,1,0,0,v]);
Append(~params20,[0,1,1,0,0,v]);
Append(~params20,[0,w,1,0,0,v]);
if p mod 4 eq 1 then
  Append(~params20,[0,w^2,1,0,0,v]);
  Append(~params20,[0,w^3,1,0,0,v]);
end if;

for y in [0..p-1] do
  Append(~params20,[0,y,1,0,1,v]);
  Append(~params20,[0,y,1,0,w,v]);
end for;

for y in [0..p-1] do
for u in [0..p-1] do
  Append(~params20,[0,y,1,1,u,v]);
  Append(~params20,[0,y,1,w,u,v]);
end for;
end for;

for y in [0..p-1] do
for t in [0..p-1] do
for u in [0..p-1] do
  Append(~params20,[1,y,1,t,u,v]);
  Append(~params20,[w,y,1,t,u,v]);
end for;
end for;
end for;

end for;


print #params20,
2*p^4+4*p^3+6*p^2+11*p+11+2*Gcd(p-1,3)+(p+1)*Gcd(p-1,4);
gtotal:=gtotal+#params20;

/*
Descendants of 5.14.

Case 21, cb=caa=cac=0, bab=baa
*/

Append(~params21,[0,0,0,0,0,0]);
Append(~params21,[1,0,0,0,0,0]);
Append(~params21,[w,0,0,0,0,0]);
Append(~params21,[0,1,0,0,0,0]);

for x in [0..p-1] do
  Append(~params21,[x,0,1,0,0,0]);
  Append(~params21,[x,0,w,0,0,0]);
  Append(~params21,[x,1,1,0,0,0]);
  Append(~params21,[x,1,w,0,0,0]);
end for;

Append(~params21,[1,0,0,1,0,0]);
Append(~params21,[w,0,0,1,0,0]);
for x in [0..p-1] do
  Append(~params21,[0,x,0,1,0,0]);
end for;

//pc=baa
for y in [0..p-1] do
for z in [0..p-1] do
for t in [0..p-1] do
  Append(~params21,[1,y,z,t,1,0]);
  Append(~params21,[w,y,z,t,1,0]);
end for;
end for;
end for;

for y in [0..p-1] do
for t in [0..p-1] do
  Append(~params21,[0,y,1,t,1,0]);
  Append(~params21,[0,y,w,t,1,0]);
end for;
end for;

for t in [0..p-1] do
  Append(~params21,[0,1,0,t,1,0]);
  Append(~params21,[0,w,0,t,1,0]);
  if p mod 4 eq 1 then
    Append(~params21,[0,w^2,0,t,1,0]);
    Append(~params21,[0,w^3,0,t,1,0]);
  end if;
end for;

Append(~params21,[0,0,0,0,1,0]);
Append(~params21,[0,0,0,1,1,0]);
Append(~params21,[0,0,0,w,1,0]);
if p mod 4 eq 1 then
  Append(~params21,[0,0,0,w^2,1,0]);
  Append(~params21,[0,0,0,w^3,1,0]);
end if;

for x in [0..p-1] do
if x eq 1 then continue; end if;
for z in [0..p-1] do
for t in [0,1] do
Append(~params21,[x,0,z,t,0,1]);
end for;
end for;
end for;

for x in [0..p-1] do
if x eq w then continue; end if;
for z in [0..p-1] do
for t in [0,1] do
Append(~params21,[x,0,z,t,0,w]);
end for;
end for;
end for;

half:=(p+1) div 2;
wover2:=F!(w/2); wover2:=Z!wover2;
for z in [0..p-1] do
  if z eq half then continue; end if;
  Append(~params21,[1,0,z,0,0,1]);
  Append(~params21,[1,1,z,0,0,1]);
end for;
for z in [0..p-1] do
  if z eq wover2 then continue; end if;
  Append(~params21,[w,0,z,0,0,w]);
  Append(~params21,[w,1,z,0,0,w]);
end for;

Append(~params21,[1,0,half,0,0,1]);
Append(~params21,[w,0,wover2,0,0,w]);
Append(~params21,[1,0,half,1,0,1]);
Append(~params21,[w,0,wover2,1,0,w]);
for t in [0..p-1] do
  Append(~params21,[1,1,half,t,0,1]);
  Append(~params21,[w,1,wover2,t,0,w]);
end for;

print #params21,
2*p^3+6*p^2+7*p+7+(p+1)*Gcd(p-1,4);
gtotal:=gtotal+#params21;

/*
Descendants of 5.14.

Case 22, cb=baa=caa=0, cac=wbab
*/

mats:=[];
//get pb,pc

bcrange:=[[0,1]];
for i in [0..p-1] do
  Append(~bcrange,[1,i]);
end for;

for y3 in [0,1] do
for y4 in [0..p-1] do
for y5 in [0..p-1] do
for y6 in [0..p-1] do

A:=H22![y3,y4,y5,y6];
if A eq 0 then
  Append(~mats,A);
  continue;
end if;

new:=1;
index:=p^3*y3+p^2*y4+p*y5+y6;

for bc in bcrange do
b:=bc[1]; c:=bc[2];
for m in [-1,1] do

B:=H22![w*b,m*c,w*c,w*m*b];
C:=H22![w*(w*b^2+c^2),m*2*w*b*c,2*w^2*b*c,m*w*(w*b^2+c^2)];

D:=B*A*C^-1;
if D[1][1] ne 0 then D:=D[1][1]^-1*D; end if;
if D[1][1] eq 0 and D[1][2] ne 0 then D:=D[1][2]^-1*D; end if;
if D[1][1] eq 0 and D[1][2] eq 0 and D[2][1] ne 0 then D:=D[2][1]^-1*D; end if;
if D[1][1] eq 0 and D[1][2] eq 0 and D[2][1] eq 0 and D[2][2] ne 0 then D:=D[2][2]^-1*D; end if;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);


ind1:=p^3*z3+p^2*z4+p*z5+z6;

if ind1 lt index then new:=0; end if;


if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~mats,A);
  //print y3,y4,y5,y6;
end if;

end for;
end for;
end for;
end for;

for AS in mats do
//Get stabilizer of AS
stab1:=[];
stab2:=[];
for b in [0..p-1] do
for c in [0..p-1] do
for m in [-1,1] do
if b+c eq 0 then continue; end if;

B:=H22![w*b,m*c,w*c,w*m*b];
C:=H22![w*(w*b^2+c^2),m*2*w*b*c,2*w^2*b*c,m*w*(w*b^2+c^2)];

D:=B*AS*C^-1;
if D[1][1] ne 0 then D:=D[1][1]^-1*D; end if;
if D[1][1] eq 0 and D[1][2] ne 0 then D:=D[1][2]^-1*D; end if;
if D[1][1] eq 0 and D[1][2] eq 0 and D[2][1] ne 0 then D:=D[2][1]^-1*D; end if;
if D[1][1] eq 0 and D[1][2] eq 0 and D[2][1] eq 0 and D[2][2] ne 0 then D:=D[2][2]^-1*D; end if;

if D eq AS then
  Append(~stab1,H33![1,0,0,0,w*b,m*c,0,w*c,w*m*b]);
  Append(~stab2,C);
end if;

end for;
end for;
end for;

y3:=Z!(AS[1][1]); y4:=Z!(AS[1][2]); y5:=Z!(AS[2][1]); y6:=Z!(AS[2][2]);

for y1 in [0..p-1] do
for y2 in [0..p-1] do

A:=H32![y1,y2,y3,y4,y5,y6];
if A eq 0 then
  Append(~params22,[0,0,0,0,0,0]);
  continue;
end if;

new:=1;
index:=p*y1+y2;

for ii in [1..#stab1] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);

ind1:=p*z1+z2;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;


if new eq 1 then
  Append(~params22,[y1,y2,y3,y4,y5,y6]);
end if;

end for;
end for;

end for;

print #params22,
(2*p^3+3*p^2+3*p+13-Gcd(p-1,3)+(p+1)*Gcd(p-1,4))/2;
gtotal:=gtotal+#params22;

/*
Descendants of 5.14.

Case 23, cb=baa=0, caa=bac, cac=wbab
*/
sol:=0;
if p mod 3 eq 2 then
  //look for solution of 12wx^2=-1
  for x in [1..p-1] do
    if F!(12*w*x^2+1) eq 0 then sol:=x; break; end if;
  end for;
end if;

mats:=[];
//get pc,pb

for y5 in [0,1,lns] do
for y6 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y3,y4,y5,y6];
if A eq 0 then
  Append(~mats,A);
  continue;
end if;

new:=1;
index:=p^3*y5+p^2*y6+p*y3+y4;

for a in [1..p-1] do
for m in [-1,1] do

B:=H22![a,0,0,m*a];
C:=H22![a^3,0,0,m*a^3];

D:=B*A*C^-1;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);


ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=0; end if;

if p mod 3 eq 2 and new eq 1 then
for b in [sol,-sol] do

B:=H22![-2*w*a*b,a,m*w*a,-2*m*w*a*b];
C:=H22![8/3*w^2*a^3*b,4/3*w*a^3,4/3*m*w^2*a^3,8/3*m*w^2*a^3*b];

D:=B*A*C^-1;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);

ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;

end if;

if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~mats,A);
  //print y3,y4,y5,y6;
end if;

end for;
end for;
end for;
end for;

for AS in mats do
//Get stabilizer of AS
stab1:=[];
stab2:=[];
for a in [1..p-1] do
for m in [-1,1] do

B:=H22![a,0,0,m*a];
C:=H22![a^3,0,0,m*a^3];

D:=B*AS*C^-1;

if D eq AS then
  Append(~stab1,H33![a,0,0,0,a,0,0,0,m*a]);
  Append(~stab2,C);
end if;

if p mod 3 eq 2 then
for b in [sol,-sol] do

B:=H22![-2*w*a*b,a,m*w*a,-2*m*w*a*b];
C:=H22![8/3*w^2*a^3*b,4/3*w*a^3,4/3*m*w^2*a^3,8/3*m*w^2*a^3*b];

D:=B*AS*C^-1;

if D eq AS then
  Append(~stab1,H33![4*w*a*b,-3*w*a*b,a/2,0,-2*w*a*b,a,0,m*w*a,-2*m*w*a*b]);
  Append(~stab2,C);
end if;

end for;
end if;

end for;
end for;

y3:=Z!(AS[1][1]); y4:=Z!(AS[1][2]); y5:=Z!(AS[2][1]); y6:=Z!(AS[2][2]);

for y1 in [0..p-1] do
for y2 in [0..p-1] do

A:=H32![y1,y2,y3,y4,y5,y6];
if A eq 0 then
  Append(~params23,[0,0,0,0,0,0]);
  continue;
end if;

new:=1;
index:=p*y1+y2;

for ii in [1..#stab1] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);

ind1:=p*z1+z2;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;


if new eq 1 then
  Append(~params23,[y1,y2,y3,y4,y5,y6]);
end if;

end for;
end for;

end for;

if p mod 3 eq 1 then
  expect:=p^5+p^4+p^3+p^2+p+2+(p^2+p+1)*Gcd(p-1,4)/2;
else;
  expect:=p^5/3+p^4/3+p^3/3+p^2/3+p+2+(p^2+p+1)*Gcd(p-1,4)/2;
end if;
print #params23,expect;
gtotal:=gtotal+#params23;

/*
Descendants of 5.14.

Case 24, cb=baa=0, caa=xbab+bac, cac=wbab
where x is not a value of y(y^2+3w)/(3y^2+w)

p=2 mod 3 only
*/

if p mod 3 eq 2 then

val:={};
for x in [0..p-1] do
  if F!(3*x^2+w) ne 0 then
    a:=F!(x*(x^2+3*w))*(F!(3*x^2+w))^-1;
    Include(~val,Z!a);
  end if;
end for;
for a in [1..p-1] do
  if a notin val then s:=a; break; end if;
end for;
//This is the value of s we need for the presentation

//look for solution of wb^2=-3
for y in [1..p-1] do
  if F!(w*y^2+3) eq 0 then b:=y; break; end if;
end for;

B1:=H22![2,2*b,2*w*b,2];
C1:=H22![32,-32*b,-32*w*b,32];
B2:=H22![2,-2*b,-2*w*b,2];
C2:=H22![32,32*b,32*w*b,32];
BB1:=H33![-4,s*b+3,3*s*w^-1+b,0,2,2*b,0,2*w*b,2];
BB2:=H33![-4,-s*b+3,3*s*w^-1-b,0,2,-2*b,0,-2*w*b,2];

mats:=[];
//get pc,pb

for y5 in [0,1,lns] do
for y6 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

A:=H22![y3,y4,y5,y6];
if A eq 0 then
  Append(~mats,A);
  continue;
end if;

new:=true;
index:=p^3*y5+p^2*y6+p*y3+y4;

A1:=B1*A*C1^-1;
A2:=B2*A*C2^-1;

for a in [1..(p-1) div 2] do

D:=(F!a)^-2*A;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);


ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=false; break; end if;

D:=(F!a)^-2*A1;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);

ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=false; break; end if;

D:=(F!a)^-2*A2;

z3:=Z!(D[1][1]);
z4:=Z!(D[1][2]);
z5:=Z!(D[2][1]);
z6:=Z!(D[2][2]);

ind1:=p^3*z5+p^2*z6+p*z3+z4;

if ind1 lt index then new:=false; break; end if;

end for;

if new then
  Append(~mats,A);
  //print y3,y4,y5,y6;
end if;

end for;
end for;
end for;
end for;


for AS in mats do
//Get stabilizer of AS
stab1:=[];
stab2:=[];

AS1:=B1*AS*C1^-1;
AS2:=B2*AS*C2^-1;

for a in [1..(p-1) div 2] do

D:=(F!a)^-2*AS;

if D eq AS then
  Append(~stab1,H33![1,0,0,0,1,0,0,0,1]);
  Append(~stab2,H22![a^2,0,0,a^2]);
end if;

D:=(F!a)^-2*AS1;

if D eq AS then
  Append(~stab1,BB1);
  Append(~stab2,a^2*C1);
end if;

D:=(F!a)^-2*AS2;

if D eq AS then
  Append(~stab1,H33!BB2);
  Append(~stab2,a^2*C2);
end if;

end for;

y3:=Z!(AS[1][1]); y4:=Z!(AS[1][2]); y5:=Z!(AS[2][1]); y6:=Z!(AS[2][2]);

for y1 in [0..p-1] do
for y2 in [0..p-1] do

A:=H32![y1,y2,y3,y4,y5,y6];
if A eq 0 then
  Append(~params24,[0,0,0,0,0,0]);
  continue;
end if;

new:=1;
index:=p*y1+y2;

for ii in [1..#stab1] do

B:=stab1[ii];
C:=stab2[ii];

D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);

ind1:=p*z1+z2;

if ind1 lt index then new:=0; end if;

if new eq 0 then break; end if;
end for;


if new eq 1 then
  Append(~params24,[y1,y2,y3,y4,y5,y6]);
end if;

end for;
end for;

end for;

print #params24,2*(p^5+p^4+p^3+p^2)/3+2*p+3;
gtotal:=gtotal+#params24;

end if;

print "Algebra 5.14 has",gtotal,"descendants of order p^7 and class 3";
print "2p^5+7p^4+19p^3+49p^2+128p+256+(p^2+7p+29)gcd(p-1,3)+
(p^2+7p+24)gcd(p-1,4)+(p+3)gcd(p-1,5) =",
2*p^5+7*p^4+19*p^3+49*p^2+128*p+256+(p^2+7*p+29)*Gcd(p-1,3)+
(p^2+7*p+24)*Gcd(p-1,4)+(p+3)*Gcd(p-1,5);


print Cputime(tt);

