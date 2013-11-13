//Descendants of 4.1 with L^2 of order p^2: Case 4
//[d,a]=[c,b]=1, [d,b]=[c,a], [d,c]=[b,a]^w
/*
This program computes the orbits of the non-singular 2x2 matrices
needed for algebra 7.311

The orbit reps are stored in mats.

We let S be the set of non-singular 2x2 matrices and we let G be
the group of non-singular matrices of the form

                 [a      b   ]
                 [+/-wb  +/-a]

(Note that this matrix is non-singular unless a = b = 0, so that
G has order 2(p^2 - 1).)

We let G L (Z_p\{0}) act on S as follows: if P in G and c in Z_p\{0},
and A in S then (P,c) sends A to

            cPAP^-1.

We attach a numerical index to each matrix, and we
print out the matrix in each orbit of lowest index.

*/

readi p,"Input prime p";
//Get primitive element w
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

range:={[0,1]};
for i in [0..p-1] do
  Include(~range,[1,i]);
end for;

Z:=Integers();
M:=MatrixAlgebra(F,2);
A:=M!0;
B:=M!0;
mats:=[];
for s in range do
x:=s[1]; y:=s[2];
A[1][1]:=F!x;
A[1][2]:=F!y;
for z in [0..p-1] do
A[2][1]:=F!z;
for t in [0..p-1] do
A[2][2]:=F!t;
if Rank(A) eq 2 then
index:=p^3*x+p^2*y+p*z+t;
new:=1;

for r in range do
a:=r[1]; b:=r[2];

B[1][1]:=F!a;
B[2][2]:=F!a;
B[1][2]:=F!b;
B[2][1]:=F!(w*b);
C:=B*A*B^-1;

u:=C[1][1];
if u eq 0 then u:=C[1][2]; end if;
C:=u^-1*C;
  
x1:=Z!(C[1][1]);
y1:=Z!(C[1][2]);
z1:=Z!(C[2][1]);
t1:=Z!(C[2][2]);
ind1:=p^3*x1+p^2*y1+p*z1+t1;

if ind1 lt index then new:=0; break; end if;

B[2][1]:=-B[2][1];
B[2][2]:=-B[2][2];
C:=B*A*B^-1;

u:=C[1][1];
if u eq 0 then u:=C[1][2]; end if;
C:=u^-1*C;
  
x1:=Z!(C[1][1]);
y1:=Z!(C[1][2]);
z1:=Z!(C[2][1]);
t1:=Z!(C[2][2]);
ind1:=p^3*x1+p^2*y1+p*z1+t1;

if ind1 lt index then new:=0; break; end if;

end for;

if new eq 1 then
  n:=#mats;
  mats[n+1]:=[x,y,z,t];
end if;

end if;
end for;
end for;
end for;

print #mats;
print "(p+1)^2/2 =",(p+1)^2/2;
