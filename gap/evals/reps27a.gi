BindGlobal( "GetReps1Slow", function(F, W, A1)
    local A, B, C, D, E, R, n, a, c, x;

    B := [[1,0],[0,1]] * One(F);
    C := [[1,0,0],[0,1,0],[0,0,1]] * One(F);
    if A1 = 0*A1 then return [1,B,C]; fi;

    R := [[0,0,1],[0,1,0],[1,0,0]] * One(F);
    A := A1*R;
    A:=EchelonForm(A); B := A[2]; A := A[1];
    A := A*R;
    D := [[0,1],[1,0]] * One(F);
    B := D*B;
    A := D*A;
    C := [[1,0,0],[0,1,0],[0,0,1]] * One(F);
    if A[2][3] <> Zero(F) then
        if A[1][2] <> Zero(F) then
            n:=A[1][1];
            if n <> Zero(F) then
                D:=[[1,0],[2*n,1]] * One(F);
                B:=D*B;
                A:=D*A;
                C:=[[1,0,0],[n,1,0],[n^2,2*n,1]] * One(F);
                A:=A*C^-1;
            fi;
            a:=A[2][1];
            if a = Zero(F) then return [9,B,C]; fi;
            if IsBool(IsSquareGF(F, a)) then a:=a*W^-1; fi;
            x:=IsSquareGF(F, a^-1);
            D:=[[1,0,0],[0,x,0],[0,0,x^2]] * One(F);
            C:=D*C;
            E:=[[x,0],[0,x^2]] * One(F);
            B:=E*B;
            A:=E*A*D^-1;
            if A[2][1] = One(F) then return [10,B,C]; fi;
            return [11,B,C];
       else
            n:=A[2][2]/2;
            c:=A[2][2]*n-n^2;
            D:=[[1,0],[c,1]] * One(F);
            E:=[[1,0,0],[n,1,0],[n^2,2*n,1]] * One(F);
            A:=D*A*E^-1;
            B:=D*B;
            C:=E*C;
       fi;
       if A[2][1] <> Zero(F) then
           a:=A[2][1];
           if IsBool(IsSquareGF(F, a)) then a:=a*W^-1; fi;
           x:=IsSquareGF(F, a^-1);
           D:=[[1,0,0],[0,x,0],[0,0,x^2]] * One(F);
           C:=D*C;
           E:=[[1,0],[0,x^2]] * One(F);
           B:=E*B;
           A:=E*A*D^-1;
           if A[2][1] = One(F) then return [7,B,C]; fi;
           return [8,B,C];
        fi;
        if A[1][1] = Zero(F) then return [6,B,C]; fi;
        return [5,B,C];
    fi;
    if A[2][2] <> Zero(F) and A[2][1] <> Zero(F) then
        n:=A[2][1];
        E:=[[1,0,0],[n,1,0],[n^2,2*n,1]] * One(F);
        A:=A*E^-1;
        C:=E*C;
    fi;
    if A[2][2] <> Zero(F) and A[1][1] <> Zero(F) then return [4,B,C]; fi;
    if A[2][2] <> Zero(F) then return [3,B,C]; fi;
    return [2,B,C];
end );

BindGlobal( "GetReps1", function(P, A, repstable)
   local E, F, B, index;
   E := NullMat(2, 3, GF(P));
   E[1]:=A[1]; E[2]:=A[2];
   F := EchelonForm(E);
   B := F[2]; F := List(F[1], IntVecFFE);
   index:=1+P^3*(2*F[1][1]+F[2][2])+P^2*F[1][2]+P*F[1][3]+F[2][3];
   return [repstable[index],index,B];
end );

#############################################################################

BindGlobal( "GetReps2", function(P)
local reps2, W, F, A, B, q, x, l;
reps2:=[];
F:=GF(P);
W:=PrimitiveRootMod(P);
A:=NullMat(4,3,F);
A[2][1]:=1*One(F);
B:=List(A,ShallowCopy);
for q in [0..P-1] do
  B[3][2]:=q*One(F);
  for x in [0,1] do
    B[4][1]:=x*One(F);
    Add(reps2,B);
    B:=List(B,ShallowCopy);
  od;
od;
B:=List(A,ShallowCopy);
B[3][1]:=1*One(F);
B[3][2]:=1*One(F);
Add(reps2,B);
B:=List(A,ShallowCopy);
B[4][2]:=1*One(F);
Add(reps2,B);
B:=List(B,ShallowCopy);
B[4][1]:=1*One(F);
Add(reps2,B);
B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
B[4][1]:=1*One(F);
Add(reps2,B);
B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
for x in [0..P-1] do
  B[3][1]:=x*One(F);
  Add(reps2,B);
  B:=List(B,ShallowCopy);
od;
B:=List(A,ShallowCopy);
B[3][1]:=1*One(F);
B[4][3]:=1*One(F);
for q in [0..P-1] do
  B[3][2]:=q*One(F);
  for x in [0..P-1] do
    B[4][1]:=x*One(F);
    Add(reps2,B);
    B:=List(B,ShallowCopy);
  od;
od;
B[3][1]:=0*One(F);
for q in [0..P-1] do
  B[3][2]:=q*One(F);
  for x in [0,1,W] do
    B[4][1]:=x*One(F);
    Add(reps2,B);
    B:=List(B,ShallowCopy);
  od;
od;
B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
B[4][2]:=1*One(F);
for x in [0..P-1] do
  B[3][1]:=x*One(F);
  Add(reps2,B);
  B:=List(B,ShallowCopy);
od;

Sort(reps2);
l := P^2+7*P+4;
if Length(reps2) <> l then Error("reps2 wrong"); fi;
#if Length(reps2) <> Length(Set(reps2)) then Error("reps2 has duplicates"); fi;
return reps2;
end );

#############################################################################

BindGlobal( "GetReps3", function(P)
local reps3, F, A, B, x, q, y, z, l, W; 
reps3:=[];
W := PrimitiveRootMod(P);
F := GF(P);
A:=NullMat(4,3,F);
A[2][2]:=1*One(F);
B:=List(A,ShallowCopy);
B[3][2]:=1*One(F);
for x in [0..P-1] do
for q in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][3]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B[3][2]:=0*One(F);
for x in [0,1,W] do
for q in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][3]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B[3][1]:=0*One(F);
B[4][1]:=1*One(F);
for x in [0,1] do
for q in [0..P-1] do
  B[3][2]:=x*One(F);
  B[3][3]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B[3][2]:=0*One(F);
B[4][2]:=1*One(F);
for x in [0..P-1] do
for q in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][3]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B[4][1]:=0*One(F);
for x in [0,1,W] do
for q in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][3]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B:=List(A,ShallowCopy);
B[4][3]:=1*One(F);
for x in [0,1,W] do
  B[3][1]:=x*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
for x in [0..P-1] do
for q in [1,W] do
  B[3][1]:=x*One(F);
  B[4][1]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B[3][2]:=1*One(F);
for x in [0..P-1] do
for q in [0..P-1] do
  B[3][1]:=x*One(F);
  B[4][1]:=q*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
B[4][2]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
for z in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  Add(reps3,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

Sort(reps3);
l := P^3+3*P^2+10*P+3;
if Length(reps3) <> l then Error("reps3 wrong"); fi;
#if Length(reps3) <> Length(Set(reps3)) then Error("reps3 has duplicates"); fi;
return reps3;
end );

#############################################################################

BindGlobal( "GetReps4", function(P)
local reps4, F, W, A, B, x, y, z, t, l;
reps4:=[];
W:=PrimitiveRootMod(P);
F:=GF(P);
A:=NullMat(4,3,F);
A[1][1]:=1*One(F);
A[2][2]:=1*One(F);
B:=List(A,ShallowCopy);
B[4][3]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
for z in [0..P-1] do
for t in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  B[4][2]:=t*One(F);
  Add(reps4,B);
  B:=List(B,ShallowCopy);
od;
od;
od;
od;
B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
B[4][1]:=1*One(F);
B[4][2]:=-1*One(F);
Add(reps4,B);
B:=List(B,ShallowCopy);
B[4][1]:=0*One(F);
for x in [0,1,W] do
  B[3][1]:=x*One(F);
  Add(reps4,B);
  B:=List(B,ShallowCopy);
od;
B[3][1]:=0*One(F);
B[3][2]:=1*One(F);
for x in [0..P-1] do
  B[4][1]:=x*One(F);
  Add(reps4,B);
  B:=List(B,ShallowCopy);
od;
B[4][1]:=-1*One(F);
for x in [1..P-1] do
  B[3][1]:=x*One(F);
  Add(reps4,B);
  B:=List(B,ShallowCopy);
od;
B:=List(A,ShallowCopy);
for x in [0,1,W] do
  B[3][1]:=x*One(F);
  B[3][3]:=0*One(F);
  B[4][2]:=1*One(F);
  Add(reps4,B);
  B:=List(B,ShallowCopy);
  for y in [2..P-1] do
    B[3][3]:=y*One(F);
    B[4][2]:=(1-2*y)*One(F);
    Add(reps4,B);
    B:=List(B,ShallowCopy);
  od;
od;
B:=List(A,ShallowCopy);
B[3][2]:=1*One(F);
for x in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][3]:=0*One(F);
  B[4][2]:=1*One(F);
  Add(reps4,B);
  B:=List(B,ShallowCopy);
  for y in [2..P-1] do
    B[3][3]:=y*One(F);
    B[4][2]:=(1-2*y)*One(F);
    Add(reps4,B);
    B:=List(B,ShallowCopy);
  od;
od;
B:=List(A,ShallowCopy);
for x in [0,1,W] do
  B[3][1]:=x*One(F);
  for y in [0..P-1] do
    B[3][3]:=y*One(F);
    for z in [0..P-1] do
      if z <> (1-2*y) mod P then 
        B[4][2]:=z*One(F);
        Add(reps4,B);
        B:=List(B,ShallowCopy);
      fi;
    od;
  od;
od;
B:=List(A,ShallowCopy);
B[4][1]:=1*One(F);
for x in [0..P-1] do
  B[3][1]:=x*One(F);
  for y in [0..P-1] do
    B[3][3]:=y*One(F);
    for z in [0..P-1] do
      if z <> (1-2*y) mod P then 
        B[4][2]:=z*One(F);
        Add(reps4,B);
        B:=List(B,ShallowCopy);
      fi;
    od;
  od;
od;

Sort(reps4);
l := P^4+P^3+3*P^2+P;
if Length(reps4) <> l then Error("reps4 wrong"); fi;
#if Length(reps4) <> Length(Set(reps4)) then Error("reps4 has duplicates"); fi;
return reps4;
end );

#############################################################################

#//Just to be perverse, we want these to be listed in a rather
#//twisted lexicographic order.  The representatives are
#//integer sequences [0,1,0,0,0,1,x,y,z,t,u,v] and We Want
#//to order them lexicographically according to the sequence
#//[v,z,t,y,x,u].  The reason is that I can see a Way to generate
#//them in that order, Whereas if I Wanted the list in lexicographic
#//order I Would have to generate the Whole list and then sort it.
#//At p=19 this Would take 2 minutes on my laptop

BindGlobal( "GetReps5", function(P)
local reps5, W, x, y, u, t, z, l, xrange;
reps5:=[];
W := PrimitiveRootMod(P);
xrange:=[0,1]; if P mod 3 = 1 then xrange:=[0,1,W,W^2] mod P; fi; Sort(xrange);
for x in xrange do
for u in [0..P-1] do
  Add(reps5,MyCutVector([1,0,0,0,0,1,x,0,0,0,u,0], 4)*One(GF(P)));
od;
od;

for y in [1,W] do
for x in [0..(P-1)/2] do
for u in [0..P-1] do
  Add(reps5,MyCutVector([1,0,0,0,0,1,x,y,0,0,u,0], 4)*One(GF(P)));
od;
od;
od;

for y in [0..P-1] do
for x in [0..P-1] do
for u in [0..P-1] do
  Add(reps5,MyCutVector([1,0,0,0,0,1,x,y,0,1,u,0], 4)*One(GF(P)));
od;
od;
od;

for t in [0..P-1] do
for y in [0..P-1] do
for x in [0..P-1] do
for u in [0..P-1] do
  Add(reps5,MyCutVector([1,0,0,0,0,1,x,y,1,t,u,0], 4)*One(GF(P)));
od;
od;
od;
od;

for z in [0..P-1] do
for t in [0..P-1] do
for y in [0..P-1] do
for x in [0..P-1] do
for u in [0..P-1] do
  Add(reps5,MyCutVector([1,0,0,0,0,1,x,y,z,t,u,1], 4)*One(GF(P)));
od;
od;
od;
od;
od;

l := P^5+P^4+P^3+P^2+2*P+P*Gcd(P-1,3);
if Length(reps5) <> l then Error("reps5 wrong"); fi;
#if Length(reps5) <> Length(Set(reps5)) then Error("reps5 has duplicates"); fi;
return reps5;
end );

#############################################################################

BindGlobal( "GetReps6", function(P)
local reps6, F, W, A, B, x, y, z, l;
reps6:=[];
W:=PrimitiveRootMod(P);
F:=GF(P);
A:=NullMat(4,3,F);
A[2][3]:=1*One(F);
B:=List(A,ShallowCopy);
B[4][2]:=1*One(F);
B[4][3]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
for z in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B[4][2]:=0*One(F);
for x in [0..(P-1)/2] do
for y in [1,W] do
for z in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B[3][2]:=0*One(F);
for x in [0..(P-1)/2] do
for y in [1,W] do
  B[3][1]:=x*One(F);
  B[4][1]:=y*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;
od;
B:=List(A,ShallowCopy);
B[4][3]:=1*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
B[3][1]:=1*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
if P mod 3 = 1 then
  B[3][1]:=W*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
  B[3][1]:=W^2*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
fi;

B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
B[4][2]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
  B[3][1]:=x*One(F);
  B[4][1]:=y*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;
od;

B:=List(A,ShallowCopy);
B[4][1]:=1*One(F);
B[4][2]:=1*One(F);
for x in [0..P-1] do
  B[3][1]:=x*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;

B[4][1]:=0*One(F);
B[3][1]:=0*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
B[3][1]:=1*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
if P mod 3 = 1 then
  B[3][1]:=W*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
  B[3][1]:=W^2*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
fi;
B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
B[4][1]:=1*One(F);
for x in [0..P-1] do
  B[3][2]:=x*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;

B[3][3]:=0*One(F);
B[3][2]:=0*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
B[3][2]:=1*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
B[3][2]:=W*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);

B:=List(A,ShallowCopy);
B[3][3]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;
od;

B[3][3]:=0*One(F);
for x in [0..(P-1)/2] do
for y in [1,W] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
od;
od;

B[3][2]:=0*One(F);
B[3][1]:=0*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
B[3][1]:=1*One(F);
Add(reps6,B);
B:=List(B,ShallowCopy);
if P mod 3 = 1 then
  B[3][1]:=W*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
  B[3][1]:=W^2*One(F);
  Add(reps6,B);
  B:=List(B,ShallowCopy);
fi;

Sort(reps6);
l := P^3+3*P^2+5*P+8+3*Gcd(P-1,3);
if Length(reps6) <> l then Error("reps6 wrong"); fi;
#if Length(reps6) <> Length(Set(reps6)) then Error("reps6 has duplicates"); fi;
return reps6;
end );

#############################################################################

BindGlobal( "GetReps7", function(P)
local reps7, F, W, A, B, x, y, z, t, l;
reps7:=[];
W:=PrimitiveRootMod(P);
F:=GF(P);
A:=NullMat(4,3,F);
A[2][1]:=1*One(F);
A[2][3]:=1*One(F);
B:=List(A,ShallowCopy);
B[4][3]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
for z in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  B[4][2]:=t*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;
od;
od;

B[4][2]:=0*One(F);
for y in [0..P-1] do
for z in [0..P-1] do
for t in [0..(P-1)/2] do
  B[3][1]:=t*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B:=List(A,ShallowCopy);
B[4][2]:=1*One(F);
for y in [0..P-1] do
for z in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=y*One(F);
  B[3][3]:=t*One(F);
  B[4][1]:=z*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;
od;
B:=List(A,ShallowCopy);
B[4][2]:=1*One(F);
for y in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=y*One(F);
  B[4][1]:=t*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;

B[4][1]:=0*One(F);
for t in [0..(P-1)/2] do
  B[3][1]:=t*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;

B:=List(A,ShallowCopy);
B[4][1]:=1*One(F);
for y in [0..P-1] do
for t in [0..(P-1)/2] do
  B[3][2]:=y*One(F);
  B[3][3]:=t*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;

B:=List(A,ShallowCopy);
for x in [0..P-1] do
for y in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[3][3]:=t*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B[3][3]:=0*One(F);
for y in [0..P-1] do
for t in [0..(P-1)/2] do
  B[3][1]:=t*One(F);
  B[3][2]:=y*One(F);
  Add(reps7,B);
  B:=List(B,ShallowCopy);
od;
od;

Sort(reps7);
l := (P^2+1)*(P+1)^2/2;
if Length(reps7) <> l then Error("reps7 wrong"); fi;
#if Length(reps7) <> Length(Set(reps7)) then Error("reps7 has duplicates"); fi;
return reps7;
end );

#############################################################################

BindGlobal( "GetReps8", function(P)
local reps8, W, F, A, B, x, y, z, t, l;
reps8:=[];
W:=PrimitiveRootMod(P);
F:=GF(P);
A:=NullMat(4,3,F);
A[2][1]:=W*One(F);
A[2][3]:=1*One(F);
B:=List(A,ShallowCopy);
B[4][3]:=1*One(F);
for x in [0..P-1] do
for y in [0..P-1] do
for z in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  B[4][2]:=t*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;
od;
od;

B[4][2]:=0*One(F);
for y in [0..P-1] do
for z in [0..P-1] do
for t in [0..(P-1)/2] do
  B[3][1]:=t*One(F);
  B[3][2]:=y*One(F);
  B[4][1]:=z*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B:=List(A,ShallowCopy);
B[4][2]:=1*One(F);
for y in [0..P-1] do
for z in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=y*One(F);
  B[3][3]:=t*One(F);
  B[4][1]:=z*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B:=List(A,ShallowCopy);
B[4][2]:=1*One(F);
for y in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=y*One(F);
  B[4][1]:=t*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;

B[4][1]:=0*One(F);
for t in [0..(P-1)/2] do
  B[3][1]:=t*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;

B:=List(A,ShallowCopy);
B[4][1]:=1*One(F);
for y in [0..P-1] do
for t in [0..(P-1)/2] do
  B[3][2]:=y*One(F);
  B[3][3]:=t*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;

B:=List(A,ShallowCopy);
for x in [0..P-1] do
for y in [0..P-1] do
for t in [1..(P-1)/2] do
  B[3][1]:=x*One(F);
  B[3][2]:=y*One(F);
  B[3][3]:=t*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;
od;

B[3][3]:=0*One(F);
for y in [0..P-1] do
for t in [0..(P-1)/2] do
  B[3][1]:=t*One(F);
  B[3][2]:=y*One(F);
  Add(reps8,B);
  B:=List(B,ShallowCopy);
od;
od;

Sort(reps8);
l := (P^2+1)*(P+1)^2/2;
if Length(reps8) <> l then Error("reps8 wrong"); fi;
#if Length(reps8) <> Length(Set(reps8)) then Error("reps8 has duplicates"); fi;
return reps8;
end );

#############################################################################
#//Just to be perverse, we want these to be listed in a rather
#//twisted lexicographic order.  The representatives are
#//integer sequences [0,1,0,0,0,1,x,y,z,t,u,v] and we want
#//to order them lexicographically according to the sequence
#//[z,u,y,t,x,v].  The reason is that I can see a way to generate
#//them in that order, whereas if I wanted the list in lexicographic
#//order I would have to generate the whole list and then sort it.
#//At p=19 sorting a list of this size takes two minutes on my laptop

BindGlobal( "GetReps9", function(P)
local reps9, W, xrange, x, v, t, y, l, u;
reps9:=[];
W := PrimitiveRootMod(P);
xrange:=[0,1];
if P mod 3 = 1 then xrange:=[0,1,W,W^2] mod P; fi; Sort(xrange);
for x in xrange do
for v in [0..P-1] do
  Add(reps9,MyCutVector([0,1,0,0,0,1,x,0,0,0,0,v], 4)*One(GF(P)));
od;
od;

for t in [1,W] do
for x in [0..(P-1)/2] do
for v in [0..P-1] do
  Add(reps9,MyCutVector([0,1,0,0,0,1,x,0,0,t,0,v], 4)*One(GF(P)));
od;
od;
od;

for y in [1,W] do
for t in [0..P-1] do
for x in [0..(P-1)/2] do
for v in [0..P-1] do
  Add(reps9,MyCutVector([0,1,0,0,0,1,x,y,0,t,0,v], 4)*One(GF(P)));
od;
od;
od;
od;

for y in [0..P-1] do
for t in [0..P-1] do
for x in [0..P-1] do
for v in [0..P-1] do
  Add(reps9,MyCutVector([0,1,0,0,0,1,x,y,0,t,1,v], 4)*One(GF(P)));
od;
od;
od;
od;
for u in [0..P-1] do
for y in [0..P-1] do
for t in [0..P-1] do
for x in [0..P-1] do
for v in [0..P-1] do
  Add(reps9,MyCutVector([0,1,0,0,0,1,x,y,1,t,u,v], 4)*One(GF(P)));
od;
od;
od;
od;
od;

l := P^5+P^4+P^3+2*P^2+2*P+P*Gcd(P-1,3);
if Length(reps9) <> l then Error("reps9 wrong"); fi;
#if Length(reps9) <> Length(Set(reps9)) then Error("reps9 has duplicates"); fi;
return reps9;
end );

#############################################################################

BindGlobal( "GetReps10", function(P)
local reps10, x, zrange, y, z, urange, t, u, v, l; 
reps10:=[];
for x in [0..(P-1)/2] do
zrange:=[0..P-1];
if x = 0 then zrange:=[0..(P-1)/2]; fi;
for y in [0..P-1] do
for z in zrange do
urange:=[0..P-1];
if x+z = 0 then urange:=[0..(P-1)/2]; fi;
for t in [0..P-1] do
for u in urange do
for v in [0..P-1] do
Add(reps10,MyCutVector([0,1,0,1,0,1,x,y,z,t,u,v], 4)*One(GF(P)));
od;
od;
od;
od;
od;
od;

l := P^3*(P^3+1)/2;
if Length(reps10) <> l then Error("reps10 wrong"); fi;
#if Length(reps10)<>Length(Set(reps10)) then Error("reps10 has duplicates"); fi;
return reps10;
end );

#############################################################################

BindGlobal( "GetReps11", function(P)
local reps11, W, x, zrange, y, z, urange, t, u, v, l;
reps11:=[];
W := PrimitiveRootMod(P);
for x in [0..(P-1)/2] do
zrange:=[0..P-1];
if x = 0 then zrange:=[0..(P-1)/2]; fi;
for y in [0..P-1] do
for z in zrange do
urange:=[0..P-1];
if x+z = 0 then urange:=[0..(P-1)/2]; fi;
for t in [0..P-1] do
for u in urange do
for v in [0..P-1] do
Add(reps11,MyCutVector([0,1,0,W,0,1,x,y,z,t,u,v], 4)*One(GF(P)));
od;
od;
od;
od;
od;
od;

l := P^3*(P^3+1)/2;
if Length(reps11) <> l then Error("reps11 wrong"); fi;
#if Length(reps11)<>Length(Set(reps11)) then Error("reps11 has duplicates"); fi;
return reps11;
end );

