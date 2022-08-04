

RequirePackage("liepring");
RequirePackage("anupq");

BindGlobal( "CheckFamily", function(L, P)
    local f, a, b;
    if P = 3 and PClassOfLiePRing(L) > 2 then return true; fi;
    f := NumberOfLiePRingsInFamily(L);
    a := EvaluatePorcPoly(f,P);
    b := LiePRingsInFamily(L,P);
    if b = fail then b := 0; else b := Length(b); fi;
    return a=b;
end );

BindGlobal( "CheckPorc", function(d, P)
    local LL, cc, L1, L2, f1, f2, a, b1, b2, i, L;

    LL := LiePRingsByLibrary(d);
    cc := List(LL, x -> []);

    # two special cases
    if d = 7 then 
        if P = 3 then return true; fi;
        L1 := LL[1798];
        L2 := LL[1799];
        f1 := NumberOfLiePRingsInFamily(L1);
        f2 := NumberOfLiePRingsInFamily(L2);
        a := EvaluatePorcPoly(f1+f2,P);
        b1 := Length(LiePRingsInFamily(L1,P));
        b2 := Length(LiePRingsInFamily(L2,P));
        if a <> b1+b2 then 
            Error("dim ",d," num ",1798,"+",1799," prime ",P); 
        fi;
        Print("dim ",d," special with prime ",P," done \n");
    fi;

    # all other cases
    for i in [1..Length(LL)] do
        if d<>7 or not (i in [1798,1799]) then  
            L := LL[i];
            if CheckFamily(L, P) = false then 
                Error("dim ",d," num ",i," prime ",P); 
            fi;
            Print("dim ",d," num ",i," prime ",P," done \n");
        fi;
    od;
end );

#for d in [2..7] do
#    for P in [3,5,7,11,13,17,19] do
#        CheckPorc(d,P);
#    od;
#od;

BindGlobal( "CheckCodes", function(d, P)
    local LL, cc, dd, i, c, j, k;
    LL := LiePRingsByLibrary(d);
    cc := List(LL, x -> []);
    for i in [1..Length(LL)] do
        c := LiePRingsInFamily(LL[i], P, "group");
        if not IsBool(c) then 
            for j in [1..Length(c)] do
                c[j] := StandardPresentation(c[j]);
                c[j] := PcGroupFpGroup(c[j]);
                c[j] := CodePcGroup(c[j]);
                for k in [1..i-1] do
                    if c[j] in cc[k] then Error("isom found",i," -> ",k); fi;
                od;
                for k in [1..j-1] do
                    if c[j] = c[k] then Error("isom found",i," -> ",i); fi;
                od;
            od;
            cc[i] := c;
        fi;
        Print(i, " done \n");
    od;
    cc := Sum(List(cc, Length));
    dd := NumberSmallGroups(P^d);
    Print(" number Lie rings ",cc," -- number groups ",dd,"\n");
end );


BindGlobal( "NrZerosByPrime", function( arg )
    local pp, vars, polys, p, l, z, v, a, b, e, w, W;
    if Length(arg) = 3 then 
        pp := arg[1];
        polys := arg[2];
        p := arg[3];
    else
        polys := arg[1];
        p := arg[2];
    fi;
    vars := Set(Flat(List(polys, VarsOfPoly)));
    w := IndeterminateByName("w");
    W := PrimitiveRootMod(p);
    l := Length(pp);
    e := Length(vars);
    z := 0;
    for v in GF(p)^e do
        b := Concatenation(IntVecFFE(v),[W]);
        a := List(polys, f -> Value(f, Concatenation(vars,[w]), b) mod p);
        if a = 0*a then z := z + 1; fi;
    od;
    return p^(l-e)*z;
end );

BindGlobal( "TryInterpol", function( pp, polys, n )
    local pr, wt, i, p, r, a;

    # get values
    pr := [7];
    wt := [];
    for i in [1..n] do
        wt[i] := NrZerosByPrime(pp, polys, pr[i]);
        pr[i+1] := NextPrimeInt(pr[i]);
    od;

    # sum up
    p := IndeterminateByName("p");
    r := 0;
    for i in [1..n] do
        a := Difference([1..n], [i]);
        a := Product(List(a, k -> (p-pr[k])/(pr[i]-pr[k])));
        r := r + wt[i]*a;
    od;
    return r;
end );

BindGlobal( "Check2", function(p)
    local w, z, a, o;

    w := PrimitiveRootMod(p)*One(GF(p));
    z := Zero(GF(p));
    o := One(GF(p));

    a := Length(Filtered(GF(p)^2, x -> x[1]*(x[2]-x[2]^2)-1 = z));
    if a <> p-2 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> -x[1]*(x[2]-x[2]^2)+1=z));
    if a <> p-2 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[1]^2-x[1]*x[2]+x[2]=z));
    if a <> p-1 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[1]*(x[2]-1)=z and x[1]*(x[1]-o)=z));
    if a <> p+1 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[1]*(x[2]+o)=z and x[1]*(x[1]-o)=z));
    if a <> p+1 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[2]^2-o=z and x[1]-x[2]=z));
    if a <> 2 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[2]^2-o=z and w^2*x[1]^2-x[2]=z));
    if a <> Gcd(p-1,4) then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[1]^2+x[2]+o=z));
    if a <> p then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[1]^2+w+x[2]=z));
    if a <> p then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[2]^2-3*x[2]+3*o=z and x[1]+x[2]-4*o=z));
    if a <> Gcd(p-1,3)-1 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> x[2]^2-3*o=z and x[1]+x[2]-o=z));
    if a <> 2-(Gcd(p-1,3)-Gcd(p-1,4)+1)^2/2 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> (x[1]-x[2])*(x[2]+o)=z));
    if a <> 2*p-1 then Error("hier"); fi;

    a := Length(Filtered(GF(p)^2, x -> (w*x[1]-x[2])*(x[2]+o)=z));
    if a <> 2*p-1 then Error("hier"); fi;

end );

BindGlobal( "Check3", function(p)
    local w, v, o, a;

    w := PrimitiveRootMod(p)*One(GF(p));
    v := Zero(GF(p));
    o := One(GF(p));

    a := Length(Filtered(GF(p)^3, x -> (x[3]+o)*(x[1]*x[2]-x[3])=v));
    if a <> 2*p^2-p+1 then Error("hier 1"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[3]+o=v and x[1]*x[2]+o=v));
    if a <> p-1 then Error("hier 2"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[1]*x[2]=x[3]));
    if a <> p^2 then Error("hier 3"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[1]*x[2]=v));
    if a <> 2*p^2-p then Error("hier 4"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[3]=v and x[1]*x[2]=v));
    if a <> 2*p-1 then Error("hier 5"); fi;

    a := Length(Filtered(GF(p)^3, x -> -x[3]^3+w*x[3]=1/2*x[1]));
    if a <> p^2 then Error("hier 6"); fi;

    a := Length(Filtered(GF(p)^3, x -> w^3*x[3]-1/2*w^2*x[1]=x[3]^3));
    if a <> p^2 then Error("hier 7"); fi;

    a := Length(Filtered(GF(p)^3, x -> -x[2]*x[3]+x[1]=v and w=x[2]));
    if a <> p then Error("hier 8"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[2]*x[3]=v and x[1]=v and w=x[2]));
    if a <> 1 then Error("hier 9"); fi;

    a := Length(Filtered(GF(p)^3, x -> (x[2]-o)*((x[1]-x[2])*(x[3]-o)-o)=v));
    if a <> p^2 + (p-1)^2 then Error("hier 10"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[2]=v and x[1]*x[3]=x[1]+o));
    if a <> p-1 then Error("hier 11"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[3]=v and 
                                         x[1]*x[2]-x[2]^2-x[1]+2*x[2]=o));
    if a <> 2*p-1 then Error("hier 12"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[2]=o and x[1]*x[3]-x[1]=x[3]));
    if a <> p-1 then Error("hier 13"); fi;

    a := Length(Filtered(GF(p)^3, x -> (x[1]-x[2])*(x[3]-1)=o));
    if a <> p*(p-1) then Error("hier 14"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[2]*x[3]+x[2]=o and x[1]=v));
    if a <> p-1 then Error("hier 15"); fi;

    a := Length(Filtered(GF(p)^3, x -> (x[3]+o)*(x[1]+x[2])=o));
    if a <> p*(p-1) then Error("hier 16"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[2]=v and x[1]*(x[3]+1)=o));
    if a <> p-1 then Error("hier 17"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[3]^2+2*x[3]+2*o=v and 
                                               x[2]=v and x[1]+x[3]=-o));
    if a <> Gcd(p-1,4)-2 then Error("hier 18"); fi;

    a := Length(Filtered(GF(p)^3, x -> x[2]=o and (x[1]+1)*(x[3]+1)=o));
    if a <> p-1 then Error("hier 19"); fi;

end );


BindGlobal( "Check4", function(p)
    local w, v, o, a;

    w := PrimitiveRootMod(p)*One(GF(p));
    v := Zero(GF(p));
    o := One(GF(p));

    a := Length(Filtered(GF(p)^4, b -> b[1]*b[2]=b[3]*b[4]));
    if a <> p^3+p^2-p then Error("hier 1"); fi;

    a := Length(Filtered(GF(p)^4, x -> x[4]=v and x[2]*x[3]=v));
    if a <> p*(2*p-1) then Error("hier 2"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[4]=v and b[2]*b[3]=v and 
                                                           b[1]*b[2]=v));
    if a <> p^2+p-1 then Error("hier 3"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[4]=v and b[2]*b[3]=v and b[1]=v));
    if a <> 2*p-1 then Error("hier 4"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[2]*b[3]+b[4]^2=v and b[1]+b[4]=v));
    if a <> p^2 then Error("hier 5"); fi;

    a := Length(Filtered(GF(p)^4, b -> (b[1]+b[4])*(b[1]*b[4]-b[2]*b[3])=v));
    if a <> 2*p^3-p then Error("hier 6"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[3]=v and b[4]*b[1]*(b[1]+b[4])=v));
    if a <> p*(3*p-2) then Error("hier 7"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[3]=v and b[1]*(b[1]+b[4])=v));
    if a <> p*(2*p-1) then Error("hier 8"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[3]=v and b[2]=v and 
                                                    b[1]*(b[1]+b[4])=v));
    if a <> (2*p-1) then Error("hier 9"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[2]*b[3]*b[4]=v and b[1]=v));
    if a <> 3*p^2-3*p+1 then Error("hier 10"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[1]*b[4]=b[2]*b[3] and 
             w*b[2]*b[4]=b[1]*b[3] and b[3]*(w*b[2]^2-b[1]^2)=v));
    if a <> 2*p^2-1 then Error("hier 11"); fi;

    a := Length(Filtered(GF(p)^4, b -> (b[1]+b[4])*(b[1]*b[4]-b[2]*b[3])=v
         and w*b[2]*b[4]=b[1]*b[3] 
         and (w*b[2]^2-b[1]^2)=(b[1]*b[4]-b[2]*b[3])));
    if a <> 2*p^2-1 then Error("hier 12"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[4]*(b[3]-b[4])=v and 
         b[4]*(b[1]-b[2])=v and (b[1]-1)*b[3]-(b[2]-1)*b[4]=v and 
         (b[1]-1)*(b[1]-b[2])=v));
    if a <> 2*p^2-1 then Error("hier 13"); fi;

    a := Length(Filtered(GF(p)^4, b -> b[2]*b[3]+b[4]^2=v and 
         b[1]+b[4]=v and w*b[4]^2=b[3]^2 and w*b[2]+b[3]=v));
    if a <> 1 then Error("hier 14"); fi;



end );

