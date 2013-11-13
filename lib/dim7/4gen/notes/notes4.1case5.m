//Revised version 17th August 2013. It seems to work
/*
This program computes the orbits of the 4x3 matrices
needed for Case 5 in the descendants of 4.1 with
L^2 of order p^3.  Read notes4.1.pdf to find what is going on!

We have 4x3 matrices A over GF(p) with rows representing pa,pb,pc,pd
so that

         [pa] = A[ba]
         [pb]    [ca]
         [pc]    [dc]
         [pd]

We get the isomorphism classes of algebras satisfying these
commutator relations by considering the orbits of matrices A
under the action A -> BAC^-1 where

  B = [al,bl,bm,-am]
      [cl,dl,dm,-cm]
      [cn,dn,dx,-cx]
      [-an,-bn,-bx,ax]

  C =  (ad-bc)[l^2,2lm,m^2]
              [ln,lx+mn,mx]
              [n^2,2nx,x^2]

with ad-bc ne 0, lx-mn ne 0.

If we consider the orbits of matrices A under the subgroup H consisting of
matrices B,C above with m=0 then we can "write down" a set orbit representatives.
(There are p^6+... of them!)  It turns out that there are 11 2x3 matrices
A1,A2,...,A11 such that every matrix A lies in the same H orbit as a matrix
with first two rows equal to one of A1,A2,...,A11. Furthermore a matrix with
first two rows Ai can only be in the same H-orbit as a matrix with first
two rows Aj if i=j.

For each Ai (i=1,2,...,11) we "write down" a set of representatives for the
H-orbits containing a matrix A with first two rows Ai. These are stored
in reps1,reps2,...,reps11. A1=0, and reps1 is just the 11 matrices with
first two rows equal to zero, and rows 3 and 4 equal to Ai (i=1,2,...,11).

Actually, we don't store reps10 and reps11, as they have size p^3*(p^3+1)/2.
For p=19 it takes 4.5 gigabytes to store all the representatives.

We have a function getreps1slow which, given a matrix A, returns a matrices
D,B,C (with (B,C) in H) such that D:=B*A*C^-1, and such that D has first
two rows equal to Ai for some i.

We use this function to build a table of triples D,B,C for all A such that
the first two rows of A are in reduced row echelon form. (So the table is
not too big.) We then use this table to define a function getreps1 which
(effectively) returns D,B,C for ALL A.  (getreps1 is quite a bit faster
than getreps1slow!)

The 11 matrices in reps1 are all in separate orbits under the action of the
full group G.

We then take reps2, reps3, ... in turn, to see how these H-orbits
fuse together under the action of G.

We attach (distinct) indices to all the representatives from
reps2, reps3, ... , reps11 so that the index attached to a representative 
from repsi is less than the index attached to repsj whenever i<j.

We then take reps2, reps3, ... in turn, to see how these H-orbits
fuse together under the action of G.

Suppose we are dealing with repsi.

For each A in repsi, we hit it (in turn) with elements for a transveral
for H in G, giving a matrix D.  (Since there are p^6 + ... H-orbit representatives, 
this means we have to deal with p^7 + ... matrices D in all.)

The essential idea is to find the index of the H-orbit representative
of the H-orbit containing D. If this index is less than the index
of A we throw A away, and go on to consider the next representative
from repsi. If we haven't gone on to the next representative from repsi
after hitting A with all the elements from the transversal, then
we add A to our list of G-orbit representatives.

The aim in the paragraph above is achieved as follows. First we use
getreps1 to find E,B,C with E=BDC^-1, with the first two rows
of E equal to Aj for some j. If the j<i we can abandon A, since
whatever the index of the H-orbit representive of D, it is bound
to be less than the index of A. If j>i then we go on to the next
transversal element, since whatever the index of the H-orbit
representative of D, it is bound to be bigger than the index
of A. If i=j we use the function getindexi to find the index
of the H-orbit representative of D.

Better look at the code now, or refer back to notes4.1.pdf.
*/


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

Z:=Integers();
V2:=VectorSpace(F,2);
V3:=VectorSpace(F,3);
V4:=VectorSpace(F,4);
H22:=Hom(V2,V2);
H23:=Hom(V2,V3); // 2x3 
H33:=Hom(V3,V3);
H43:=Hom(V4,V3);
H44:=Hom(V4,V4);


//There is a slight gain in storing the table below, rather
//than getting Magma to cube roots as you need them.

//Get table of cube roots
curoots:=[0:i in [1..p-1]];
for i in [1..p-1] do
  j:=i^3 mod p;
  curoots[j]:=i;
end for;

//We need to sort some of reps2, reps3 etc into lexicographic order,
//but Magma doesn't know how to sort matrices over GF(p).
sortmats:=function(s)
  t:=[];
  for i in [1..#s] do
    A:=s[i];
    t[i]:=[Z!A[1][1],Z!A[1][2],Z!A[1][3],Z!A[2][1],Z!A[2][2],Z!A[2][3],
           Z!A[3][1],Z!A[3][2],Z!A[3][3],Z!A[4][1],Z!A[4][2],Z!A[4][3]];
  end for;
  Sort(~t);
  u:=[];
  for i in [1..#t] do
    Append(~u,H43!t[i]);
  end for;
  return u;
end function;

//These are the matrices A1,A2,...,A11
reps1:=[H23![0,0,0,0,0,0],
        H23![0,0,0,1,0,0],
        H23![0,0,0,0,1,0],
        H23![1,0,0,0,1,0],
        H23![1,0,0,0,0,1],
        H23![0,0,0,0,0,1],
        H23![0,0,0,1,0,1],
        H23![0,0,0,w,0,1],
        H23![0,1,0,0,0,1],
        H23![0,1,0,1,0,1],
        H23![0,1,0,w,0,1]];


//Function to get the H-orbit represntative of a 2x3 matrix A
//representing pa,pb, together with H-matrices which transform
//A to the representative.
getreps1slow:=function(A1);
B:=H22![1,0,0,1];
C:=H33![1,0,0,0,1,0,0,0,1];
if A1 eq 0 then return 1,B,C; end if;
//Fiddle matrix about so as to be able to use built in echelon form
R:=H33![0,0,1,0,1,0,1,0,0];
A:=A1*R;
A,B:=EchelonForm(A);
B:=H22!B;
A:=A*R;
D:=H22![0,1,1,0];
B:=D*B;
A:=D*A;
C:=H33![1,0,0,0,1,0,0,0,1];
if A[2][3] ne 0 then
  if A[1][2] ne 0 then
    n:=A[1][1];
    if n ne 0 then
      D:=H22![1,0,2*n,1];
      B:=D*B;
      A:=D*A;
      C:=H33![1,0,0,n,1,0,n^2,2*n,1];
      A:=A*C^-1;
    end if;
    a:=A[2][1];
    if a eq 0 then return 9,B,C; end if;
    if not IsSquare(a) then a:=a*w^-1; end if;
    x:=Sqrt(a^-1);
    D:=H33![1,0,0,0,x,0,0,0,x^2];
    C:=D*C;
    E:=H22![x,0,0,x^2];
    B:=E*B;
    A:=E*A*D^-1;
    if A[2][1] eq 1 then return 10,B,C; end if;
    return 11,B,C;
  end if;
  if A[1][2] eq 0 then
    n:=A[2][2]/2;
    c:=A[2][2]*n-n^2;
    D:=H22![1,0,c,1];
    E:=H33![1,0,0,n,1,0,n^2,2*n,1];
    A:=D*A*E^-1;
    B:=D*B;
    C:=E*C;
  end if;
  if A[2][1] ne 0 then
    a:=A[2][1];
    if not IsSquare(a) then a:=a*w^-1; end if;
    x:=Sqrt(a^-1);
    D:=H33![1,0,0,0,x,0,0,0,x^2];
    C:=D*C;
    E:=H22![1,0,0,x^2];
    B:=E*B;
    A:=E*A*D^-1;
    if A[2][1] eq 1 then return 7,B,C; end if;
    return 8,B,C;
  end if;
  if A[1][1] eq 0 then return 6,B,C; end if;
  return 5,B,C;
end if;
if A[2][2] ne 0 and A[2][1] ne 0 then
  n:=A[2][1];
  E:=H33![1,0,0,n,1,0,n^2,2*n,1];
  A:=A*E^-1;
  C:=E*C;
end if;
if A[2][2] ne 0 and A[1][1] ne 0 then return 4,B,C; end if;
if A[2][2] ne 0 then return 3,B,C; end if;
return 2,B,C;
end function;

/*
//test getreps1slow to see if it is correct.
for A in H23 do
n,B,C:=getreps1slow(A);
if B*A*C^-1 ne reps1[n] then print A,B*A*C^-1,reps1[n]; readi z; end if;
end for;
*/

//Use getreps1slow to
//make a look up table for finding the H-orbit representative
//for the first two rows of a matrix A.  Also store the H-matrices
//which transform A to a matrix where the first two rows equal
//this represenative.  Save space (and time at this step)
//by only creating a table for matrices where the first two
//rows are in echelon form.
repstable:=[];
Btable:=[];
Ctable:=[];
for v in [0,1] do
for y in [0,1] do
for z in [0..p-1] do
for t in [0..p-1] do
for u in [0..p-1] do
  A:=H23![v,z,t,0,y,u];
  if EchelonForm(A) ne A then continue; end if;
  index:=1+p^3*(2*v+y)+p^2*z+p*t+u;
  n,B,C:=getreps1slow(A);
  repstable[index]:=n;
  a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
  n:=C[2][1]; x:=C[2][2];
  B:=H44![a,b,0,0,c,d,0,0,c*n,d*n,d*x,-c*x,-a*n,-b*n,-b*x,a*x];
  Btable[index]:=B;
  Ctable[index]:=C;
end for;
end for;
end for;
end for;
end for;

//Now use this table to define a function which given A, finds
//a matrix D in the same H-orbit as A, with first two rows
//equal to Ai for some i.

//The function returns i,index,B.

//Here i = repstable[index]

//B is a 2x2 matrix which transforms the first two rows of A into echelon
//form. If needed you can use B to construct a 4x4 matrix E with (E,I) in H,
//such that the first two rows of EA are in echelon form.

//index tells you where to look in Btable,Ctable to find B,C
//such that the first two rows of BEAC^-1 equal Ai.

getreps1:=function(A)
  E:=H23!0;
  E[1]:=A[1]; E[2]:=A[2];
  F,B:=EchelonForm(E);
  index:=1+p^3*(2*Z!F[1][1]+Z!F[2][2])+p^2*Z!F[1][2]+p*Z!F[1][3]+Z!F[2][3];
  return repstable[index],index,B;
end function;

/*
//Test getreps1;
for A in H23 do
  B:=Random(H23);
  A1:=H43!0;
  A1[1]:=A[1];
  A1[2]:=A[2];
  A1[3]:=B[1];
  A1[4]:=B[2];
  n,index,B:=getreps1(A);
  a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
  B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
  A1:=B*A1;
  B:=Btable[index];
  C:=Ctable[index];
  A1:=B*A1*C^-1;
  B:=H23!0;
  B[1]:=A1[1];
  B[2]:=A1[2];
  if B ne reps1[n] then print "Arghhhh"; end if;
end for;
*/

//Generate reps2
reps2:=[];
A:=H43!0;
A[2][1]:=1;
B:=A;
for q in [0..p-1] do
  B[3][2]:=q;
  for x in [0,1] do
    B[4][1]:=x;
    Append(~reps2,B);
  end for;
end for;
B:=A;
B[3][1]:=1;
B[3][2]:=1;
Append(~reps2,B);
B:=A;
B[4][2]:=1;
Append(~reps2,B);
B[4][1]:=1;
Append(~reps2,B);
B:=A;
B[3][3]:=1;
B[4][1]:=1;
Append(~reps2,B);
B:=A;
B[3][3]:=1;
for x in [0..p-1] do
  B[3][1]:=x;
  Append(~reps2,B);
end for;
B:=A;
B[3][1]:=1;
B[4][3]:=1;
for q in [0..p-1] do
  B[3][2]:=q;
  for x in [0..p-1] do
    B[4][1]:=x;
    Append(~reps2,B);
  end for;
end for;
B[3][1]:=0;
for q in [0..p-1] do
  B[3][2]:=q;
  for x in [0,1,w] do
    B[4][1]:=x;
    Append(~reps2,B);
  end for;
end for;
B:=A;
B[3][3]:=1;
B[4][2]:=1;
for x in [0..p-1] do
  B[3][1]:=x;
  Append(~reps2,B);
end for;

reps2:=sortmats(reps2);
//We need these representatives sorted into lexicographic order
//print #reps2,p^2+7*p+4;


//The function getindex2 finds the index of the H-orbit representative
//of a matrix with first two rows equal to A2
getindex2:=function(A)
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v ne 0 then
  x:=x+(u^2*z+u*v-2*t*v*z-u*v*y)/(2*v^2);
  y:=y-u*z/v;
  t:=t/v-u^2/(4*v^2);
  if x ne 0 then
    t:=t/x^2;
    return p^5+p^4*Z!y+p^2*Z!t+1;
  end if;
  if t eq 0 then return p^4*Z!y+1; end if;
  if IsSquare(t) then return p^4*Z!y+p^2+1; end if;
  return p^4*Z!y+p^2*w+1;
end if;

if z ne 0 and u ne 0 then
  x:=z*x+(t^2*z+t*u-t*u*y)*z/u^2;
  return p^5*Z!x+p^3+p;
end if;

if z ne 0 and u eq 0 then
  if t ne 0 then return p^3+p^2; end if;
  x:=x*z+y/2-y^2/4;
  return p^5*Z!x+p^3;
end if;

if u ne 0 then
  t:=t-t*y+x*u;
  if t eq 0 then return p; end if;
  return p^2+p;
end if;

if t ne 0 then return p^4*Z!y+p^2; end if;

if y ne 1 or x eq 0 then return p^4*Z!y; end if;

return p^5+p^4;

end function; 


//Generate reps3
reps3:=[];
A:=H43!0;
A[2][2]:=1;
B:=A;
B[3][2]:=1;
for x in [0..p-1] do
for q in [0..p-1] do
  B[3][1]:=x;
  B[3][3]:=q;
  Append(~reps3,B);
end for;
end for;
B[3][2]:=0;
for x in [0,1,w] do
for q in [0..p-1] do
  B[3][1]:=x;
  B[3][3]:=q;
  Append(~reps3,B);
end for;
end for;
B[3][1]:=0;
B[4][1]:=1;
for x in [0,1] do
for q in [0..p-1] do
  B[3][2]:=x;
  B[3][3]:=q;
  Append(~reps3,B);
end for;
end for;
B[3][2]:=0;
B[4][2]:=1;
for x in [0..p-1] do
for q in [0..p-1] do
  B[3][1]:=x;
  B[3][3]:=q;
  Append(~reps3,B);
end for;
end for;
B[4][1]:=0;
for x in [0,1,w] do
for q in [0..p-1] do
  B[3][1]:=x;
  B[3][3]:=q;
  Append(~reps3,B);
end for;
end for;
B:=A;
B[4][3]:=1;
for x in [0,1,w] do
  B[3][1]:=x;
  Append(~reps3,B);
end for;
for x in [0..p-1] do
for q in [1,w] do
  B[3][1]:=x;
  B[4][1]:=q;
  Append(~reps3,B);
end for;
end for;
B[3][2]:=1;
for x in [0..p-1] do
for q in [0..p-1] do
  B[3][1]:=x;
  B[4][1]:=q;
  Append(~reps3,B);
end for;
end for;
B[4][2]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
for z in [0..p-1] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[4][1]:=z;
  Append(~reps3,B);
end for;
end for;
end for;

reps3:=sortmats(reps3);
//print #reps3,p^3+3*p^2+10*p+3;

//The function getindex3 finds the index of the H-orbit representative
//of a matrix with first two rows equal to A3
getindex3:=function(A)
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v ne 0 then
  c:=z/v;
  x:=x-c*t;
  y:=y-c*u;
  t:=t/v;
  u:=u/v;
  if u ne 0 then
    x:=x/u^2;
    y:=y/u;
    t:=t/u^2;
    return p^5*Z!x+p^4*Z!y+p^2*Z!t+p+1;
  end if;
  if y ne 0 then
    x:=x/y^2;
    t:=t/y^2;
    return p^5*Z!x+p^4+p^2*Z!t+1;
  end if;
  if t ne 0 then
    e:=t^-1;
    if not IsSquare(e) then e:=e*w; end if;
    x:=x*e;
    t:=t*e;
    return p^5*Z!x+p^2*Z!t+1;
  end if;
  if x eq 0 then return 1; end if;
  if IsSquare(x) then return p^5+1; end if;
  return p^5*w+1;
end if;

if u ne 0 then
  x:=x-t*y/u;
  t:=t/u;
  if t ne 0 then
    x:=x/t^2;
    return p^5*Z!x+p^3*Z!z+p^2+p;
  end if;
  if x eq 0 then return p^3*Z!z+p; end if;
  if IsSquare(x) then return p^5+p^3*Z!z+p; end if;
  return p^5*w+p^3*Z!z+p;
end if;

if t ne 0 then
  if y eq 0 then return p^3*Z!z+p^2; end if;
  return p^4+p^3*Z!z+p^2;
end if;

if y ne 0 then
  x:=x/y^2;
  return p^5*Z!x+p^4+p^3*Z!z;
end if;

if x eq 0 then return p^3*Z!z; end if;
if IsSquare(x) then return p^5+p^3*Z!z; end if;
return p^5*w+p^3*Z!z;

end function;


//Generate reps4
reps4:=[];
A:=H43!0;
A[1][1]:=1;
A[2][2]:=1;
B:=A;
B[4][3]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
for z in [0..p-1] do
for t in [0..p-1] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[4][1]:=z;
  B[4][2]:=t;
  Append(~reps4,B);
end for;
end for;
end for;
end for;
B:=A;
B[3][3]:=1;
B[4][1]:=1;
B[4][2]:=-1;
Append(~reps4,B);
B[4][1]:=0;
for x in [0,1,w] do
  B[3][1]:=x;
  Append(~reps4,B);
end for;
B[3][1]:=0;
B[3][2]:=1;
for x in [0..p-1] do
  B[4][1]:=x;
  Append(~reps4,B);
end for;
B[4][1]:=-1;
for x in [1..p-1] do
  B[3][1]:=x;
  Append(~reps4,B);
end for;
B:=A;
for x in [0,1,w] do
  B[3][1]:=x;
  B[3][3]:=0;
  B[4][2]:=1;
  Append(~reps4,B);
  for y in [2..p-1] do
    B[3][3]:=y;
    B[4][2]:=1-2*y;
    Append(~reps4,B);
  end for;
end for;
B:=A;
B[3][2]:=1;
for x in [0..p-1] do
  B[3][1]:=x;
  B[3][3]:=0;
  B[4][2]:=1;
  Append(~reps4,B);
  for y in [2..p-1] do
    B[3][3]:=y;
    B[4][2]:=1-2*y;
    Append(~reps4,B);
  end for;
end for;
B:=A;
for x in [0,1,w] do
  B[3][1]:=x;
  for y in [0..p-1] do
    B[3][3]:=y;
    for z in [0..p-1] do
      if z eq (1-2*y) mod p then continue; end if;
      B[4][2]:=z;
      Append(~reps4,B);
    end for;
  end for;
end for;
B:=A;
B[4][1]:=1;
for x in [0..p-1] do
  B[3][1]:=x;
  for y in [0..p-1] do
    B[3][3]:=y;
    for z in [0..p-1] do
      if z eq (1-2*y) mod p then continue; end if;
      B[4][2]:=z;
      Append(~reps4,B);
    end for;
  end for;
end for;

reps4:=sortmats(reps4);
//print #reps4,p^4+p^3+3*p^2+p;

//Function to find the index of the H-orbit representative
//of a matrix with first two rows equal to A4
getindex4:=function(A)
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v ne 0 then
  x:=v^2*x+u*z^2-t*v*z-v*y*z;
  y:=z-u*z+v*y;
  t:=t*v-z-u*z+z^2;
  u:=u-2*z;
  return p^5*Z!x+p^4*Z!y+p^2*Z!t+p*Z!u+1;
end if;

if z eq 1 and u eq -1 then
  if y ne 0 then x:=x/y^2; t:=t/y; y:=1; end if;
  if y eq 0 and t ne 0 then x:=x/t^2; t:=1; end if;
  if y+t ne 0 or x eq 0 then return p^4*Z!y+p^3+p^2*Z!t+p*Z!u; end if;
  if y eq 1 then return p^5*x+p^4+p^3+p^2*Z!t+p*Z!u; end if;
  if IsSquare(x) then return p^5+p^3+p*(p-1); end if;
  return p^5*w+p^3+p*(p-1);
end if;

if z ne 1 and u eq 1-2*z then
  x:=(t^2+2*y*t-4*x+4*x*z)/(4*z-4);
  if y ne 0 then x:=x/y^2; return p^5*Z!x+p^4+p^3*Z!z+p*Z!u; end if;
  if x eq 0 then return p^3*Z!z+p*Z!u; end if;
  if IsSquare(x) then return p^5+p^3*Z!z+p*Z!u; end if;
  return p^5*w+p^3*Z!z+p*Z!u;
end if;

x:=(x*u^2-t*u*y+4*x*u*z-2*x*u-y^2*z+y^2-2*t*y*z+t*y+4*x*z^2-4*x*z+x)*(u+2*z-1)^-2;
t:=-(t+y-t*u-2*t*z+u*y)*(u+2*z-1)^-1;
if t ne 0 then x:=x/t^2; return p^5*Z!x+p^3*Z!z+p^2+p*Z!u; end if;
if x eq 0 then return p^3*Z!z+p*Z!u; end if;
if IsSquare(x) then return p^5+p^3*Z!z+p*Z!u; end if;
return p^5*w+p^3*Z!z+p*Z!u;

end function;


//Generate reps5
//Just to be perverse, we want these to be listed in a rather
//twisted lexicographic order.  The representatives are
//integer sequences [0,1,0,0,0,1,x,y,z,t,u,v] and we want
//to order them lexicographically according to the sequence
//[v,z,t,y,x,u].  The reason is that I can see a way to generate
//them in that order, whereas if I wanted the list in lexicographic
//order I would have to generate the whole list and then sort it.
//At p=19 this would take 2 minutes on my laptop
reps5:=[];
xrange:=[0,1];
if p mod 3 eq 1 then xrange:=[0,1,w,w^2 mod p]; end if;
Sort(~xrange);
//Sometimes, I suppose, (w^2 mod p) can be less than w, and I want these
//in ascending order.
for x in xrange do
for u in [0..p-1] do
  Append(~reps5,H43![1,0,0,0,0,1,x,0,0,0,u,0]);
end for;
end for;

for y in [1,w] do
for x in [0..(p-1) div 2] do
for u in [0..p-1] do
  Append(~reps5,H43![1,0,0,0,0,1,x,y,0,0,u,0]);
end for;
end for;
end for;

for y in [0..p-1] do
for x in [0..p-1] do
for u in [0..p-1] do
  Append(~reps5,H43![1,0,0,0,0,1,x,y,0,1,u,0]);
end for;
end for;
end for;

for t in [0..p-1] do
for y in [0..p-1] do
for x in [0..p-1] do
for u in [0..p-1] do
  Append(~reps5,H43![1,0,0,0,0,1,x,y,1,t,u,0]);
end for;
end for;
end for;
end for;

for z in [0..p-1] do
for t in [0..p-1] do
for y in [0..p-1] do
for x in [0..p-1] do
for u in [0..p-1] do
  Append(~reps5,H43![1,0,0,0,0,1,x,y,z,t,u,1]);
end for;
end for;
end for;
end for;
end for;

//These are sorted in a weird way
//print #reps5,p^5+p^4+p^3+p^2+2*p+p*GCD(p-1,3);


//Function to find the index of the H-orbit representative
//of a matrix with first two rows equal to A5
getindex5:=function(A)
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v ne 0 then
  x:=x*v^3;
  y:=y*v^2;
  z:=z*v;
  t:=t*v;
  return p^5+p^4*Z!z+p^3*Z!t+p^2*Z!y+p*Z!x+Z!u;
end if;
if z ne 0 then
  x:=x/z^3;
  y:=y/z^2;
  t:=t/z;
  return p^4+p^3*Z!t+p^2*Z!y+p*Z!x+Z!u;
end if;
if t ne 0 then
  x:=x/t^3;
  y:=y/t^2;
  return p^3+p^2*Z!y+p*Z!x+Z!u;
end if;
if y ne 0 then
  e:=y^-1;
  if not IsSquare(e) then e:=e*w; end if;
  a:=Sqrt(e);
  x:=x*a^3;
  y:=y*e;
  if Z!x gt (p-1) div 2 then x:=-x; end if;
  return p^2*Z!y+p*Z!x+Z!u;
end if;
if x eq 0 then return Z!u; end if;
e:=x;
x:=1;
//if p=1 mod 3 we need both these if statements.
if curoots[Z!e] eq 0 then e:=e/w; x:=F!w; end if;
if curoots[Z!e] eq 0 then x:=x^2; end if;
return p*Z!x+Z!u;
end function;


//Generate reps6
reps6:=[];
A:=H43!0;
A[2][3]:=1;
B:=A;
B[4][2]:=1;
B[4][3]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
for z in [0..p-1] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[4][1]:=z;
  Append(~reps6,B);
end for;
end for;
end for;

B[4][2]:=0;
for x in [0..(p-1) div 2] do
for y in [1,w] do
for z in [0..p-1] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[4][1]:=z;
  Append(~reps6,B);
end for;
end for;
end for;

B[3][2]:=0;
for x in [0..(p-1) div 2] do
for y in [1,w] do
  B[3][1]:=x;
  B[4][1]:=y;
  Append(~reps6,B);
end for;
end for;

B:=A;
B[4][3]:=1;
Append(~reps6,B);
B[3][1]:=1;
Append(~reps6,B);
if p mod 3 eq 1 then
  B[3][1]:=w;
  Append(~reps6,B);
  B[3][1]:=w^2;
  Append(~reps6,B);
end if;

B:=A;
B[3][3]:=1;
B[4][2]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
  B[3][1]:=x;
  B[4][1]:=y;
  Append(~reps6,B);
end for;
end for;

B:=A;
B[4][1]:=1;
B[4][2]:=1;
for x in [0..p-1] do
  B[3][1]:=x;
  Append(~reps6,B);
end for;


B[4][1]:=0;
B[3][1]:=0;
Append(~reps6,B);
B[3][1]:=1;
Append(~reps6,B);
if p mod 3 eq 1 then
  B[3][1]:=w;
  Append(~reps6,B);
  B[3][1]:=w^2;
  Append(~reps6,B);
end if;

B:=A;
B[3][3]:=1;
B[4][1]:=1;
for x in [0..p-1] do
  B[3][2]:=x;
  Append(~reps6,B);
end for;

B[3][3]:=0;
B[3][2]:=0;
Append(~reps6,B);
B[3][2]:=1;
Append(~reps6,B);
B[3][2]:=w;
Append(~reps6,B);

B:=A;
B[3][3]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
  B[3][1]:=x;
  B[3][2]:=y;
  Append(~reps6,B);
end for;
end for;

B[3][3]:=0;
for x in [0..(p-1) div 2] do
for y in [1,w] do
  B[3][1]:=x;
  B[3][2]:=y;
  Append(~reps6,B);
end for;
end for;

B[3][2]:=0;
B[3][1]:=0;
Append(~reps6,B);
B[3][1]:=1;
Append(~reps6,B);
if p mod 3 eq 1 then
  B[3][1]:=w;
  Append(~reps6,B);
  B[3][1]:=w^2;
  Append(~reps6,B);
end if;

reps6:=sortmats(reps6);
//print #reps6,p^3+3*p^2+5*p+8+3*GCD(p-1,3);


//Function to find the index of the H-orbit representative
//of a matrix with first two rows equal to A6
getindex6:=function(A)
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v ne 0 then 
  c:=z/v;
  x:=x-c*t;
  y:=y-c*u;
  z:=0;
  t:=t/v;
  u:=u/v;
  v:=1;
  if u ne 0 then
    x:=x/u^3;
    y:=y/u^2;
    t:=t/u^2;
    u:=1;
    return p^5*Z!x+p^4*Z!y+p^2*Z!t+p+1;
  end if;
  if y ne 0 then
    e:=y^-1;
    if not IsSquare(e) then e:=e*w; end if;
    a:=Sqrt(e);
    x:=x*a^3;
    y:=y*e;
    t:=t*e;
    if Z!x gt (p-1) div 2 then x:=-x; end if;
    return p^5*Z!x+p^4*Z!y+p^2*Z!t+1;
  end if;
  if t ne 0 then
    e:=t^-1;
    if not IsSquare(e) then e:=e*w; end if;
    a:=Sqrt(e);
    x:=x*a^3;
    t:=t*e;
    if Z!x gt (p-1) div 2 then x:=-x; end if;
    return p^5*Z!x+p^2*Z!t+1;
  end if;
  if x eq 0 then return 1; end if;
  e:=x;
  x:=1;
  //if p=1 mod 3 we need both these if statements.
  if curoots[Z!e] eq 0 then e:=e/w; x:=F!w; end if;
  if curoots[Z!e] eq 0 then x:=x^2; end if;
  return p^5*Z!x+1;
end if;

if u ne 0 then
  c:=y/u;
  x:=x-c*t;
  y:=0;
  t:=t/u;
  u:=1;
  if z ne 0 then
    x:=x/z^3;
    t:=t/z;
    return p^5*Z!x+p^3+p^2*Z!t+p;
  end if;
  if t ne 0 then
    x:=x/t^3;
    return p^5*Z!x+p^2+p;
  end if;
  if x eq 0 then return p; end if;
  e:=x;
  x:=1;
  //if p=1 mod 3 we need both these if statements.
  if curoots[Z!e] eq 0 then e:=e/w; x:=F!w; end if;
  if curoots[Z!e] eq 0 then x:=x^2; end if;
  return p^5*Z!x+p;
end if;

if t ne 0 then
  //x=0,t=1,u=0,v=0
  if z ne 0 then
    y:=y/z^2;
    return p^4*Z!y+p^3+p^2;
  end if;
  if y eq 0 then return p^2; end if;
  if IsSquare(y) then return p^4+p^2; end if;
  return p^4*w+p^2;
end if;

if z ne 0 then
  x:=x/z^3;
  y:=y/z^2;
  return p^5*Z!x+p^4*Z!y+p^3;
end if;

if y ne 0 then
  e:=y^-1;
  if not IsSquare(e) then e:=e*w; end if;
  a:=Sqrt(e);
  x:=x*a^3;
  y:=y*e;
  if Z!x gt (p-1) div 2 then x:=-x; end if;
  return p^5*Z!x+p^4*Z!y;
end if;

if x eq 0 then return 0; end if;
e:=x;
x:=1;
//if p=1 mod 3 we need both these if statements.
if curoots[Z!e] eq 0 then e:=e/w; x:=F!w; end if;
if curoots[Z!e] eq 0 then x:=x^2; end if;
return p^5*Z!x;

end function;

//Generate reps7
reps7:=[];
A:=H43!0;
A[2][1]:=1;
A[2][3]:=1;
B:=A;
B[4][3]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
for z in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[4][1]:=z;
  B[4][2]:=t;
  Append(~reps7,B);
end for;
end for;
end for;
end for;

B[4][2]:=0;
for y in [0..p-1] do
for z in [0..p-1] do
for t in [0..(p-1) div 2] do
  B[3][1]:=t;
  B[3][2]:=y;
  B[4][1]:=z;
  Append(~reps7,B);
end for;
end for;
end for;

B:=A;
B[4][2]:=1;
for y in [0..p-1] do
for z in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=y;
  B[3][3]:=t;
  B[4][1]:=z;
  Append(~reps7,B);
end for;
end for;
end for;

B:=A;
B[4][2]:=1;
for y in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=y;
  B[4][1]:=t;
  Append(~reps7,B);
end for;
end for;

B[4][1]:=0;
for t in [0..(p-1) div 2] do
  B[3][1]:=t;
  Append(~reps7,B);
end for;

B:=A;
B[4][1]:=1;
for y in [0..p-1] do
for t in [0..(p-1) div 2] do
  B[3][2]:=y;
  B[3][3]:=t;
  Append(~reps7,B);
end for;
end for;

B:=A;
for x in [0..p-1] do
for y in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[3][3]:=t;
  Append(~reps7,B);
end for;
end for;
end for;

B[3][3]:=0;
for y in [0..p-1] do
for t in [0..(p-1) div 2] do
  B[3][1]:=t;
  B[3][2]:=y;
  Append(~reps7,B);
end for;
end for;

reps7:=sortmats(reps7);
//print #reps7,(p^2+1)*(p+1)^2/2;

//Function to find the index of the H-orbit representative
//of a matrix with first two rows equal to A7
//The same function works as getindex8
getindex7:=function(A)
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if v ne 0 then 
  c:=z/v;
  x:=x-c*t;
  y:=y-c*u;
  z:=0;
  t:=t/v;
  u:=u/v;
  v:=1;
  u1:=-u;
  x1:=-x;
  if [Z!u1,Z!x1] lt [Z!u,Z!x] then u:=u1; x:=x1; end if;
  return p^5*Z!x+p^4*Z!y+p^2*Z!t+p*Z!u+1;
end if;

if u ne 0 then
  c:=y/u;
  x:=x-c*t;
  y:=0;
  t:=t/u;
  u:=1;
  z1:=-z;
  t1:=-t;
  x1:=-x;
  if [Z!z1,Z!t1,Z!x1] lt [Z!z,Z!t,Z!x] then z:=z1; t:=t1; x:=x1; end if;
  return p^5*Z!x+p^3*Z!z+p^2*Z!t+p;
end if;

if t ne 0 then
  if Z!z gt (p-1) div 2 then z:=-z; end if;
  return p^4*Z!y+p^3*Z!z+p^2;
end if;

x1:=-x;
z1:=-z;
if [Z!z1,Z!x1] lt [Z!z,Z!x] then x:=x1; z:=z1; end if;
return p^5*Z!x+p^4*Z!y+p^3*Z!z;

end function;


//Generate reps8
reps8:=[];
A:=H43!0;
A[2][1]:=w;
A[2][3]:=1;
B:=A;
B[4][3]:=1;
for x in [0..p-1] do
for y in [0..p-1] do
for z in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[4][1]:=z;
  B[4][2]:=t;
  Append(~reps8,B);
end for;
end for;
end for;
end for;

B[4][2]:=0;
for y in [0..p-1] do
for z in [0..p-1] do
for t in [0..(p-1) div 2] do
  B[3][1]:=t;
  B[3][2]:=y;
  B[4][1]:=z;
  Append(~reps8,B);
end for;
end for;
end for;

B:=A;
B[4][2]:=1;
for y in [0..p-1] do
for z in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=y;
  B[3][3]:=t;
  B[4][1]:=z;
  Append(~reps8,B);
end for;
end for;
end for;

B:=A;
B[4][2]:=1;
for y in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=y;
  B[4][1]:=t;
  Append(~reps8,B);
end for;
end for;

B[4][1]:=0;
for t in [0..(p-1) div 2] do
  B[3][1]:=t;
  Append(~reps8,B);
end for;

B:=A;
B[4][1]:=1;
for y in [0..p-1] do
for t in [0..(p-1) div 2] do
  B[3][2]:=y;
  B[3][3]:=t;
  Append(~reps8,B);
end for;
end for;

B:=A;
for x in [0..p-1] do
for y in [0..p-1] do
for t in [1..(p-1) div 2] do
  B[3][1]:=x;
  B[3][2]:=y;
  B[3][3]:=t;
  Append(~reps8,B);
end for;
end for;
end for;

B[3][3]:=0;
for y in [0..p-1] do
for t in [0..(p-1) div 2] do
  B[3][1]:=t;
  B[3][2]:=y;
  Append(~reps8,B);
end for;
end for;

reps8:=sortmats(reps8);
//print #reps8,(p^2+1)*(p+1)^2/2;



//Generate reps9
//Just to be perverse, we want these to be listed in a rather
//twisted lexicographic order.  The representatives are
//integer sequences [0,1,0,0,0,1,x,y,z,t,u,v] and we want
//to order them lexicographically according to the sequence
//[z,u,y,t,x,v].  The reason is that I can see a way to generate
//them in that order, whereas if I wanted the list in lexicographic
//order I would have to generate the whole list and then sort it.
//At p=19 sorting a list of this size takes two minutes on my laptop
reps9:=[];

xrange:=[0,1];
if p mod 3 eq 1 then xrange:=[0,1,w,w^2 mod p]; end if;
Sort(~xrange);
//Sometimes, I suppose, (w^2 mod p) can be less than w, and I want these
//in ascending order.
for x in xrange do
for v in [0..p-1] do
  Append(~reps9,H43![0,1,0,0,0,1,x,0,0,0,0,v]);
end for;
end for;

for t in [1,w] do
for x in [0..(p-1) div 2] do
for v in [0..p-1] do
  Append(~reps9,H43![0,1,0,0,0,1,x,0,0,t,0,v]);
end for;
end for;
end for;

for y in [1,w] do
for t in [0..p-1] do
for x in [0..(p-1) div 2] do
for v in [0..p-1] do
  Append(~reps9,H43![0,1,0,0,0,1,x,y,0,t,0,v]);
end for;
end for;
end for;
end for;

for y in [0..p-1] do
for t in [0..p-1] do
for x in [0..p-1] do
for v in [0..p-1] do
  Append(~reps9,H43![0,1,0,0,0,1,x,y,0,t,1,v]);
end for;
end for;
end for;
end for;

for u in [0..p-1] do
for y in [0..p-1] do
for t in [0..p-1] do
for x in [0..p-1] do
for v in [0..p-1] do
  Append(~reps9,H43![0,1,0,0,0,1,x,y,1,t,u,v]);
end for;
end for;
end for;
end for;
end for;

//These reps are sorted in a wierd way
//print #reps9,p^5+p^4+p^3+2*p^2+2*p+p*GCD(p-1,3);

//Function to find the index of the H-orbit representative
//of a matrix with first two rows equal to A9
getindex9:=function(A);
x:=A[3][1];
y:=A[3][2];
z:=A[3][3];
t:=A[4][1];
u:=A[4][2];
v:=A[4][3];

if z ne 0 then
  x:=x/z^3;
  y:=y/z^2;
  t:=t/z^2;
  u:=u/z;
  return p^5+p^4*Z!u+p^3*Z!y+p^2*Z!t+p*Z!x+Z!v;
end if;
if u ne 0 then
  x:=x/u^3;
  y:=y/u^2;
  t:=t/u^2;
  return p^4+p^3*Z!y+p^2*Z!t+p*Z!x+Z!v;
end if;
if y ne 0 then
  e:=y^-1;
  if not IsSquare(e) then e:=e*w; end if;
  a:=Sqrt(e);
  x:=x*a^3;
  y:=y*e;
  t:=t*e;
  if Z!x gt (p-1) div 2 then x:=-x; end if;
  return p^3*Z!y+p^2*Z!t+p*Z!x+Z!v;
end if;
if t ne 0 then
  e:=t^-1;
  if not IsSquare(e) then e:=e*w; end if;
  a:=Sqrt(e);
  x:=x*a^3;
  t:=t*e;
  if Z!x gt (p-1) div 2 then x:=-x; end if;
  return p^2*Z!t+p*Z!x+Z!v;
end if;
if x eq 0 then return Z!v; end if;
e:=x;
x:=1;
//if p=1 mod 3 we need both these if statements.
if curoots[Z!e] eq 0 then e:=e/w; x:=F!w; end if;
if curoots[Z!e] eq 0 then x:=x^2; end if;
return p*Z!x+Z!v;
end function;

/*
//Save space by not storing reps10
//Generate reps10
reps10:=[];
for x in [0..(p-1) div 2] do
zrange:=[0..p-1];
if x eq 0 then zrange:=[0..(p-1) div 2]; end if;
for y in [0..p-1] do
for z in zrange do
urange:=[0..p-1];
if x+z eq 0 then urange:=[0..(p-1) div 2]; end if;
for t in [0..p-1] do
for u in urange do
for v in [0..p-1] do
Append(~reps10,H43![0,1,0,1,0,1,x,y,z,t,u,v]);
end for;
end for;
end for;
end for;
end for;
end for;

//These reps are sorted
//print #reps10,p^3*(p^3+1)/2;
*/


//Function to find the index of the H-orbit representative
//of a matrix with first two rows equal to A10
//This function also works as getindex11
getindex10:=function(A);
x:=Z!A[3][1];
y:=Z!A[3][2];
z:=Z!A[3][3];
t:=Z!A[4][1];
u:=Z!A[4][2];
v:=Z!A[4][3];
x1:=-x mod p;
z1:=-z mod p;
u1:=-u mod p;
if [x1,z1,u1] lt [x,z,u] then
  x:=x1; z:=z1; u:=u1;
end if;
return p^5*x+p^4*y+p^3*z+p^2*t+p*u+v;
end function;

/*
//Save space by generating the elemements of reps11 as we need them
//Generate reps11
reps11:=[];
for x in [0..(p-1) div 2] do
zrange:=[0..p-1];
if x eq 0 then zrange:=[0..(p-1) div 2]; end if;
for y in [0..p-1] do
for z in zrange do
urange:=[0..p-1];
if x+z eq 0 then urange:=[0..(p-1) div 2]; end if;
for t in [0..p-1] do
for u in urange do
for v in [0..p-1] do
Append(~reps11,H43![0,1,0,w,0,1,x,y,z,t,u,v]);
end for;
end for;
end for;
end for;
end for;
end for;

//These reps are sorted
//print #reps11,p^3*(p^3+1)/2;
*/


paramscase5:=[];

//Right transversal for H in G. (Simultaneously a left transversal.)
translft:=[];
transrt:=[];
translft[1]:=H44![0,0,0,-1,0,0,1,0,0,1,0,0,-1,0,0,0];
transrt[1]:=H33![0,0,1,0,1,0,1,0,0];
for i in [1..p-1] do
  Append(~translft,H44![1,0,0,-i,0,1,i,0,0,0,1,0,0,0,0,1]);
  Append(~transrt,H33![1,2*i,i^2,0,1,i,0,0,1]);
end for;


//reps1
A:=H43!0;
for i in [1..11] do
  B:=reps1[i];
  A[3]:=B[1];
  A[4]:=B[2];
  Append(~paramscase5,A);
end for;

print #paramscase5;

//reps2
for i in [1..#reps2] do
  A:=reps2[i];
  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 2 then new:=false; break; end if;
    if n gt 2 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex2(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;


//reps3
for i in [1..#reps3] do
  A:=reps3[i];
  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 3 then new:=false; break; end if;
    if n gt 3 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex3(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps4
for i in [1..#reps4] do
  A:=reps4[i];
  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 4 then new:=false; break; end if;
    if n gt 4 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex4(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps5
for i in [1..#reps5] do
  A:=reps5[i];
  new:=true;
  //This index is weird --- see comments on code that generates reps5
  index:=p^5*Z!A[4][3]+p^4*Z!A[3][3]+p^3*Z!A[4][1]+p^2*Z!A[3][2]+p*Z!A[3][1]+Z!A[4][2];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 5 then new:=false; break; end if;
    if n gt 5 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex5(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps6
for i in [1..#reps6] do
  A:=reps6[i];
  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 6 then new:=false; break; end if;
    if n gt 6 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex6(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps7
for i in [1..#reps7] do
  A:=reps7[i];
  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 7 then new:=false; break; end if;
    if n gt 7 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex7(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps8
for i in [1..#reps8] do
  A:=reps8[i];
  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 8 then new:=false; break; end if;
    if n gt 8 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex7(D); //Not a mistake!
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps9
for i in [1..#reps9] do
  A:=reps9[i];
  new:=true;
  //This index is weird --- see comments on code which generate reps9
  index:=p^5*Z!A[3][3]+p^4*Z!A[4][2]+p^3*Z!A[3][2]+p^2*Z!A[4][1]+p*Z!A[3][1]+Z!A[4][3];
  //Hit with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 9 then new:=false; break; end if;
    if n gt 9 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex9(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;
end for;

print #paramscase5;

//reps10
//Generate the elements of reps10 as we use them
for x in [0..(p-1) div 2] do
zrange:=[0..p-1];
if x eq 0 then zrange:=[0..(p-1) div 2]; end if;
for y in [0..p-1] do
for z in zrange do
urange:=[0..p-1];
if x+z eq 0 then urange:=[0..(p-1) div 2]; end if;
for t in [0..p-1] do
for u in urange do
for v in [0..p-1] do
A:=H43![0,1,0,1,0,1,x,y,z,t,u,v];

  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit A with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 10 then new:=false; break; end if;
    if n gt 10 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex10(D);
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;

end for;
end for;
end for;
end for;
end for;
end for;

print #paramscase5;

//reps11
//Generate elements of reps11 as we need them
for x in [0..(p-1) div 2] do
zrange:=[0..p-1];
if x eq 0 then zrange:=[0..(p-1) div 2]; end if;
for y in [0..p-1] do
for z in zrange do
urange:=[0..p-1];
if x+z eq 0 then urange:=[0..(p-1) div 2]; end if;
for t in [0..p-1] do
for u in urange do
for v in [0..p-1] do
A:=H43![0,1,0,w,0,1,x,y,z,t,u,v];

  new:=true;
  index:=p^5*Z!A[3][1]+p^4*Z!A[3][2]+p^3*Z!A[3][3]+p^2*Z!A[4][1]+p*Z!A[4][2]+Z!A[4][3];
  //Hit A with transversal for H
  for j in [1..p] do
    B:=translft[j];
    C:=transrt[j];
    D:=B*A*C^-1;
    n,ind,B:=getreps1(D);
    if n lt 11 then new:=false; break; end if;
    if n gt 11 then continue; end if;
    a:=B[1][1]; b:=B[1][2]; c:=B[2][1]; d:=B[2][2];
    B:=H44![a,b,0,0,c,d,0,0,0,0,d,-c,0,0,-b,a];
    D:=B*D;
    B:=Btable[ind];
    C:=Ctable[ind];
    D:=B*D*C^-1;
    ind1:=getindex10(D);  //Not a mistake!
    if ind1 lt index then new:=false; break; end if;
  end for;
  if new then Append(~paramscase5,A); end if;

end for;
end for;
end for;
end for;
end for;
end for;

print "Parameter list length",#paramscase5;

if p eq 3 then 
   print "Total number of algebras is ", 550; 
elif p mod 3 eq 1 then 
   print "Total number of algebras is p^5+p^4+4p^3+6p^2+18p+19 = ",  
       p^5+p^4+4*p^3+6*p^2+18*p+19; 
elif p mod 3 eq 2 then 
   print "Total number of algebras is p^5+p^4+4p^3+6p^2+16p+17 = ",  
       p^5+p^4+4*p^3+6*p^2+16*p+17; 
end if;

print Cputime(tt);

