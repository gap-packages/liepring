//Towards the end of Case 6 in the descendants of 5.1 with pL of order p^2,
//as described in note2dec5.1.pdf, we have a nasty action of GL(2,p) on
//3x2 matrices.  It turns out that every orbit has a matrix with first
//row (0,1). We find a set of representatives for the orbits, each
//representative having first row (0,1), and store the entries for
//the other two rows in a sequence mats2, with the entries in mats2
//taking the form [t,x,y,z] where row 2 is [t,x] and row 3 is [y,z]

tt:=Cputime();

readi p,"Input the prime p";
F:=GF(p);

SQ:={};
for x in [0..p-1] do
  y:=F!x;
  Include(~SQ,y^2);
end for;

for x in [2..p-1] do
  y:=F!x;
  if y notin SQ then
    ns:=x;
    break;
  end if;
end for;
//ns is least non-square

Z:=Integers();count :=0;
expect:=(p^2-1) div 2;
if (p mod 3) eq 2 then
  expect:=expect+1;
end if;

xend:=(p-1) div 2;
mats2:=[];

for t in [1,ns] do
//The experimental evidence is that t=1 is enough,
//(i.e. t=ns never arises), but I don't have a proof
//of this.  However, no time is wasted as the loop
//aborts when the expected number of orbits is found.
t1:=F!t;
for x in [0..xend] do
x1:=F!x;
for y in [1..p-1] do
y1:=F!y;
if x eq 0 then zend:=(p-1) div 2; else; zend:=p-1; end if;
for z in [0..zend] do
z1:=F!z;

//We only consider matrices where delta is not a square mod p
delta:=(t1*z1-y1*x1)^2-t1*y1;
if delta in SQ then continue; end if;

new:=1;

//We are transforming t,x,y,z by a non-singular
//2x2 matrix with entries a,b,c,d.
//We need only consider non-zero values of c, and
//we take c=1 and then there is at most one k
//such that ak,bk,k,dk is a possible matrix.

//We want a*(2*z-2*d*y-d)+b*(2*t*d^2-1-2*x*d)=0, so we choose a,b,d carefully
for d in F do
test1:=2*z1-2*d*y1-d;
test2:=2*t1*d^2-1-2*x1*d;
brange:=F;
if test1 eq 0 and test2 ne 0 then brange:=[0]; end if;
for b in brange do
arange:=F;
if test1 ne 0 then arange:=[-b*test2/test1]; end if;
for a in arange do

//check that matrix is non-singular
e:=a*d-b;
if e eq 0 then continue; end if;


k:=2*(b*d*t1-a*y1)*e^-1;
f:=e^-2*k^-1;
//multiply all entries in the matrix by f to make (1,2) entry 1

t2:=Z!(f*(d^3*t1-d*y1-d-d^2*x1+z1));
if t2 lt t then
  new:=0;
  break;
end if;
if t2 gt t then
  continue;
end if;

x2:=Z!(f*(-b*d^2*t1+b*y1+a*d+a*d^2*x1-a*z1));
if x2 lt x then
  new:=0;
  break;
end if;
if x2 gt x then
  continue;
end if;

y2:=Z!(f*(-b^2*d*t1+d*a^2*y1+a*b+b^2*x1-a^2*z1));
if y2 lt y then
  new:=0;
  break;
end if;
if y2 gt y then
  continue;
end if;

z2:=Z!(f*(b^3*t1-b*a^2*y1-a^2*b-a*b^2*x1+a^3*z1));
if z2 lt z then
  new:=0;
  break;
end if;

end for;
if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Append(~mats2,[t,x,y,z]);
  count:=count+1;
end if;

if count eq expect then break; end if;
end for;
if count eq expect then break; end if;
end for;
if count eq expect then break; end if;
end for;
if count eq expect then break; end if;
end for;

print #mats2,expect;
print Cputime(tt);
