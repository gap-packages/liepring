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
V2:=VectorSpace(F,2);
H21:=Hom(V2,V1);
H22:=Hom(V2,V2);

mats:=[];

//x here corresponds to k in the .pdf file, and to x in the GAP file
for x in [0..p-2] do
if x eq 3 then continue; end if;

C:=H22![x-1,1,-1,0];
D:=H22![x^2-2*x,x-1,1-x,-1];

//z corresponds to z in both the .pdf file and in the GAP file
for z in [0..p-1] do
A:=H21![1,z];
new:=1;

G:=C*A;
y:=G[1][1];
if y ne 0 then
  G:=y^-1*G;
  z1:=Z!(G[2][1]);
  if z1 lt z then new:=0; end if;
end if;

G:=D*A;
y:=G[1][1];
if y ne 0 then
  G:=y^-1*G;
  z1:=Z!(G[2][1]);
  if z1 lt z then new:=0; end if;
end if;

for c in [0..p-2] do
if F!(c*x+c^2-c+1) eq 0 then continue; end if;
E:=H22![(1+c*x)*(c*x-2*c+1),c*(c*x+2-c),-c*(c*x+2-c),-(-1+c)*(c+1)];
G:=E*A;
y:=G[1][1];
if y ne 0 then
  G:=y^-1*G;
  z1:=Z!(G[2][1]);
  if z1 lt z then new:=0; break; end if;
end if;
end for;

if new eq 1 then
  Append(~mats,[x,z]);
end if;

end for;

end for;

