Orbits538 := function(p)

    F := GF(p);
    w := First([2..p-1], x -> Order(x*One(GF(p)) = p-1);
    V := MatrixSpace(F, 2, 2);
    Z:=Integers();

    mats1 := [];
    for y1 in [0..(p-1) div 2] do
        for y2 in [0..p-1] do
            for y3 in [0..p-1] do
                for y4 in [0..p-1] do

                    new :=1;
                    index:=p^3*y1+p^2*y2+p*y3+y4;

                    A:=[[y1,y2],[y3,y4]]*One(F);

                    for a in [0..p-1] do
                        for b in [0..p-1] do

                            if a ne b and a ne p-b then

                                B:=[[a,b],[b,a]]*One(F);
                                C:=[[a^4-b^4,2*a*b*(a^2-b^2)],
                                    [2*a*b*(a^2-b^2),a^4-b^4]]*One(F);
                                D:=B*A*C^-1;

                                z1:=IntFFE(D[1][1]);
                                z2:=IntFFE(D[1][2]);
                                z3:=IntFFE(D[2][1]);
                                z4:=IntFFE(D[2][2]);

                                ind1:=p^3*z1+p^2*z2+p*z3+z4;

                                if ind1 lt index then new:=0; end if;

                                B:=V![a,b,-b,-a];
                                C:=V![-a^4+b^4,-2*a*b*(a^2-b^2),2*a*b*(a^2-b^2),a^4-b^4];
                                D:=B*A*C^-1;

                                z1:=Z!(D[1][1]);
                                z2:=Z!(D[1][2]);
                                z3:=Z!(D[2][1]);
                                z4:=Z!(D[2][2]);

                                ind1:=p^3*z1+p^2*z2+p*z3+z4;

                                if ind1 lt index then new:=0; end if;

                            fi;
                            if new eq 0 then break; end if;
                        od;
                        if new eq 0 then break; end if;
                    od;

                    if new eq 1 then Include(~mats1,[y1,y2,y3,y4]); end if;
                od;
            od;
            if #mats1 eq (Gcd(p-1,3)*(p^2+3*p+11)+1)/2 then break; end if;
        od;
        if #mats1 eq (Gcd(p-1,3)*(p^2+3*p+11)+1)/2 then break; end if;
    od;

//Case babb=wbaaa
mats2:={};

for y1 in [0..(p-1) div 2] do
for y2 in [0..p-1] do
for y3 in [0..p-1] do
for y4 in [0..p-1] do

new :=1;
index:=p^3*y1+p^2*y2+p*y3+y4;

A:=V![y1,y2,y3,y4];

for a in [0..p-1] do
for b in [0..p-1] do

if a+b ne 0 then

B:=V![a,b,w*b,a];
C:=V![a^4-w^2*b^4,2*a*b*(a^2-w*b^2),2*w*a*b*(a^2-w*b^2),a^4-w^2*b^4];
D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;


B:=V![a,b,-w*b,-a];
C:=V![-a^4+w^2*b^4,-2*a*b*(a^2-w*b^2),2*w*a*b*(a^2-w*b^2),a^4-w^2*b^4];
D:=B*A*C^-1;

z1:=Z!(D[1][1]);
z2:=Z!(D[1][2]);
z3:=Z!(D[2][1]);
z4:=Z!(D[2][2]);

ind1:=p^3*z1+p^2*z2+p*z3+z4;

if ind1 lt index then new:=0; end if;

end if;

if new eq 0 then break; end if;
end for;
if new eq 0 then break; end if;
end for;

if new eq 1 then
  Include(~mats2,[y1,y2,y3,y4]);
end if;

end for;
end for;
if #mats2 eq (Gcd(p-1,3)*(p^2+p+1)+5)/2 then break; end if;
end for;
if #mats2 eq (Gcd(p-1,3)*(p^2+p+1)+5)/2 then break; end if;
end for;

print #mats1,(Gcd(p-1,3)*(p^2+3*p+11)+1)/2;
print #mats2,(Gcd(p-1,3)*(p^2+p+1)+5)/2;

