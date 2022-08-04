
##
## Valuationfunctions
##

BindGlobal( "SquaresModP", function(P)
    return Set(List([0..P-1], X -> X^2 mod P));
end );

BindGlobal( "IsSquareModP", function(P, r)
    local i;
    r := r mod P;
    if r = 0 then return 0; fi;
    for i in [1..P-1] do
        if r = (i^2 mod P) then return i; fi;
    od;
    return false;
end );

BindGlobal( "LeastNonSquareModP", function(P)
    local S, i;
    S := Set(List([1..(P-1)/2], X -> (X^2) mod P));
    for i in [2..P-1] do
        if not i in S then return i; fi;
    od;
end );

##
## "x ne 0,x~x^-1" 
##
BindGlobal( "ValsFunction1", function(P, range)
    local r, i;
    r := [];
    for i in range do
        if ForAll(r, j -> (i*j mod P) <> 1) then Add(r, i); fi;
    od;
    return r;
end );

##
## "1+4x not a square"
##
BindGlobal( "ValsFunction2", function(P)
    local s, r;
    s := SquaresModP(P);
    r := Filtered( [0..P-1], j -> not (((1+4*j) mod P) in s) );
    return r;
end );

##
## "x ne 0, x~x' if x^5=x'^5 mod p"
##
BindGlobal( "ValsFunction5", function(P, e, v)
    local r, i, k, l;
    r := [1];
    l := [1];
    for i in [2..P-1] do
        k := ((i^e) mod P);
        if not k in l then
            Add(r, i);
            Add(l, k);
        fi;
    od;
    if v = 0 then Add(r, 0); fi;
    return r;
end );

##
##  "y FIXED such that 1-wy^2 is not a square mod p"
##
BindGlobal( "ValsFunction10", function(P, W)
    local s, i, j;
    s := SquaresModP(P);
    for i in [1..P-1] do
        j := (1-W*i^2) mod P;
        if not j in s then return i; fi;
    od;
end );

BindGlobal( "ValsFunction10a", function(P, W)
    local s, i, j;
    s := SquaresModP(P);
    for i in [1..P-1] do
        j := (i^2-W) mod P;
        if not j in s then return i; fi;
    od;
end );

##
##  z FIXED so that z^2-4 not a square mod p"
##
BindGlobal( "ValsFunction11", function(P)
    local s, i, j;
    s := SquaresModP(P);
    for i in [0..P-1] do
        j := (i^2-4) mod P;
        if not j in s then return i; fi;
    od;
end );

##
##  "all x,y, [x,y]~[x',y'] if y^2-wx^2=y'^2-wx'^2 mod p"
##
BindGlobal( "ValsFunction6", function(P)
    local W, r, l, i, j, k;
    W := PrimitiveRootMod(P);
    r := [];
    l := [];
    for i in [0..P-1] do
        for j in [0..P-1] do
            k := (j^2 - W * i^2) mod P;
            if not k in l then
                Add(r, [i,j]);
                Add(l, k);
            fi;
        od;
    od;
    return r;
end );

##
##  "See note6.178" and "See Notes5.3, Case 4
##
BindGlobal( "Subgroup6178", function( P )
    local U, W, a, b, m;
    U := Subgroup(GL(2,P), []);
    W := PrimitiveRootMod(P);
    for a in [0..P-1] do
        for b in [0..P-1] do
            m := [[a,b],[W*b,a]] * One(GF(P));
            if RankMat(m) = 2 then
                U := ClosureGroup(U, [m]);
                m := [[a,b],[-W*b,-a]] * One(GF(P));
                U := ClosureGroup(U, [m]);
            fi;
            if Size(U) = 2*(P^2-1) then return U; fi;
        od;
    od;
end );

BindGlobal( "Orbits6178", function( P )
    local G, U, f, e;
    G := GL(2,P);
    U := Subgroup6178(P);
    f := function( elm, mat )
       return Determinant(mat)^-1 * (mat*elm*mat^-1);
    end;
    e := Orbits(U, Elements(G), f);
    if Length(e) <> P^2 + (P+1)/2 - Gcd(P-1,4)/2 then
        Error("wrong number of orbits for 6.178");
    fi;
    return List(e, X -> X[1]);
end );

BindGlobal( "ValsFunction8", function(P)
    local e;
    e := Orbits6178( P );
    return List( e, X -> Permuted(IntVecFFE(Flat(X)),(1,4,3,2)) );
end );

BindGlobal( "ValsFunction8a", function(P)
    local e;
    e := Orbits6178( P );
    return List( e, X -> IntVecFFE(Flat(X)));
end );



##
## See note6.62
##
BindGlobal( "Subgroup662", function( P )
    local U, W, a, b, m;
    U := Subgroup(GL(2,P), []);
    W := PrimitiveRootMod(P);
    for a in [0..P-1] do
        for b in [0..P-1] do
            m := [[a,b],[W*b,a]] * One(GF(P));
            if RankMat(m) = 2 then
                U := ClosureGroup(U, [m]);
            fi;
            if Size(U) = (P^2-1) then return U; fi;
        od;
    od;
    return U;
end );

BindGlobal( "Orbits662", function( P )
    local G, U, W, f, e;
    G := GL(2,P);
    U := Subgroup662(P);
    W := PrimitiveRootMod(P);
    f := function( elm, mat )
       local a, b, new;
       a := mat[1][1]+mat[1][2]*elm[2][1];
       b := mat[1][2]*elm[2][2];
       new := [[a, b],[W*b, a]];
       return (mat*elm*new^-1);
    end;
    e := Filtered(Elements(G),
             X -> X[1][1] = One(GF(P)) and X[1][2] = Zero(GF(P)));
    e := Orbits(U, e, f);
    if Length(e) <> P then
        Error("wrong number of orbits for 6.62");
    fi;
    return List(e, X -> X[1]);
end );

BindGlobal( "ValsFunction9", function(P)
    local e;
    e := Orbits662( P );
    return List( e, X -> IntVecFFE(X[2]) );
end );

##
## See note5.38
##

BindGlobal( "Orbits538Case1", function(P)
    local e, o, r, f, i, a, b, c, d, B, C, D, j;
    e := Elements(MatrixSpace(GF(P),2,2));
    o := List(e, X-> true);
    r := [];
    f := (Gcd(P-1,3)*(P^2+3*P+11)+1)/2;
    for i in [1..Length(e)] do
        if o[i] = true then 
            Add(r, e[i]);
            if Length(r) = f then 
                return List(r, X -> IntVecFFE(Flat(X)));
            fi;
            for a in [0..P-1] do
                for b in [0..P-1] do
                    if a <> b and a <> P-b then 
                        c := a^4-b^4;
                        d := 2*a*b*(a^2-b^2);
                        B := [[a,b],[b,a]] * One(GF(P));
                        C := [[c, d], [d, c]] *One(GF(P));
                        D := B * e[i] * C^-1;
                        j := Position(e, D);
                        o[j] := i;
                        B := [[a,b],[-b,-a]] * One(GF(P));
                        C := [[-c, -d], [d, c]] *One(GF(P));
                        D := B * e[i] * C^-1;
                        j := Position(e, D);
                        o[j] := i;
                    fi;
                od;
            od;
        fi;
    od;
end );

BindGlobal( "Orbits538Case2", function(P)
    local e, o, r, f, i, a, b, c, d, B, C, D, j, W;
    e := Elements(MatrixSpace(GF(P),2,2));
    W := PrimitiveRootMod(P);
    o := List(e, X-> true);
    r := [];
    f := ((Gcd(P-1,3) * (P^2+P+1))+5)/2;
    for i in [1..Length(e)] do
        if o[i] = true then 
            Add(r, e[i]);
            if Length(r) = f then 
                return List(r, X -> IntVecFFE(Flat(X)));
            fi;
            for a in [0..P-1] do
                for b in [0..P-1] do
                    if a+b <> 0 then 
                        c := a^4 - W^2*b^4;
                        d := 2*a*b*(a^2 - W*b^2);
                        B := [[a,b],[W*b,a]] * One(GF(P));
                        C := [[c, d], [W*d, c]] *One(GF(P));
                        D := B * e[i] * C^-1;
                        j := Position(e, D);
                        o[j] := i;
                        B := [[a,b],[-W*b,-a]] * One(GF(P));
                        C := [[-c, -d], [W*d, c]] *One(GF(P));
                        D := B * e[i] * C^-1;
                        j := Position(e, D);
                        o[j] := i;
                    fi;
                od;
            od;
        fi;
    od;
end );

BindGlobal( "ValsFunction12", function(P, case)
    if case = 1 then return Orbits538Case1(P); fi;
    if case = 2 then return Orbits538Case2(P); fi;
end );

BindGlobal( "ValsFunction13", function(P)
    local o, r, i;
    o := [];
    r := [];
    for i in [1..P-1] do
        if not i in o then 
            Add(r, i);
            Append( o, [i, -i, 1/i, -1/i] mod P);
        fi;
    od;
    return r;
end );

##
## "x ne 1-w^i, x~-x+2[1-w^i]"
##
BindGlobal( "ValsFunction14", function(P,W,i)
   local r, o, j;
   r := [];
   o := (1-W^i) mod P;
   for j in [0..P-1] do
       if j <> o and j <= ((-j+2*o) mod P) then 
           Add(r, j); 
       fi;
   od;
   return r;
end );

##
## "x ne 0, all y, 4x+y^2 not a square mod p",
##
BindGlobal( "ValsFunction15", function(P)
    local s, r, X, Y;
    s := SquaresModP(P);
    r := [];
    for X in [1..P-1] do
        for Y in [0..P-1] do
            if not ((4*X+Y^2) mod P) in s then Add(r, [X,Y]); fi; 
        od;
    od;
    return r;
end );

##
##  "all x,y,z,t, [x,y,z,t]~[t+1,z+1,y-1,x-1]"
##
BindGlobal( "ValsFunction16", function(P)
    local r, x, y, z, t, a, b;
    r := [];
    for x in [0..P-1] do
        for y in [0..P-1] do
            for z in [0..P-1] do
                for t in [0..P-1] do
                    a := [x,y,z,t];
                    b := [t+1,z+1,y-1,x-1] mod P;
                    if not b in r then Add(r,a); fi;
                od;
            od;
        od;
    od;
    return r;
end );

##
##  "all x,y,z, [x,y,z]~[-x,-y,-z]"  
##
BindGlobal( "ValsFunction17", function(P)
    local a, b, c;
    a := Cartesian([1..(P-1)/2], [0..P-1], [0..P-1]);
    b := Cartesian([0], [1..(P-1)/2], [0..P-1]);
    c := Cartesian([0],[0],[0..(P-1)/2]);
    return Concatenation(a,b,c);
end );

##
##  "t ne 0, all x,z, [x,z]~[tz,x/t]"
## neu: "y ne 0, [x,y,z]~[zy,y,x/y]" 
##
BindGlobal( "ValsFunction18", function(P)
    local r, x, y, z, t, a, b;
    r := [];
    for x in [0..P-1] do
        for y in [1..P-1] do
            for z in [0..P-1] do
                a := [x,y,z];
                b := [z*y,y,x/y] mod P;
                if not b in r then Add(r,a); fi;
            od;
        od;
    od;
    return r;
end );

##
## "t ne 0,1, all x,z such that [x+t][1+z]=1 mod p, [x,z]~[tz,x/t]"
## neu: "y ne 0,1, (x+y)(1+z)=1, [x,y,z]~[zy,y,x/y]" 
##
BindGlobal( "ValsFunction18a", function(P)
    local r, x, y, z, t, a, b;
    r := [];
    for x in [0..P-1] do
        for y in [2..P-1] do
            for z in [0..P-1] do
                if (((x+y)*(1+z)) mod P) = 1 then  
                    a := [x,y,z];
                    b := [z*y,y,x/y] mod P;
                    if not b in r then Add(r,a); fi;
                fi;
            od;
        od;
    od;
    return r;
end );

##
## See Notes5.3, Case 6 and 7
##
BindGlobal( "CanoForm53", function(A, F)
    local A2, A1, r2, a, b, B;

    # cut
    A2 := A{[3,4]};
    A1 := A{[1,2]};
    r2 := RankMat(A2);

    # the easy cases
    if r2 = 0 then return A; fi;
    if r2 = 2 then return Concatenation(0*A1, A2); fi;

    # otherwise there is some work to do
    if A2[1] <> 0*A2[1] then 
        a := NormedRowVector(A2[1]);
    else
        a := NormedRowVector(A2[2]);
    fi;
    if a[1] <> 0*a[1] then 
        b := [0,1]*One(F);
    else
        b := [1,0]*One(F);
    fi;
    B := 0*A1;
    B[1] := SolutionMat([b,a],A[1])[1] * b;
    B[2] := SolutionMat([b,a],A[2])[1] * b;
    return Concatenation(B,A2);
end );

BindGlobal( "QuotElms", function(a, F)
    if a[1] <> 0*a[1] then 
        return List(Elements(F), X -> X*[0,1]);
    else
        return List(Elements(F), X -> X*[1,0]);
    fi;
end );

BindGlobal( "CanoForms53", function(A2, F)
    local c, a, e, w, v;

    c := [];
    if A2[1] <> 0 * A2[1] then 
        a := A2[1];
    else
        a := A2[2];
    fi;
    e := QuotElms(a, F);
    for w in e do
        for v in e do
            Add(c, [w,v,A2[1],A2[2]]);
        od;
    od;

    return c;
end );

BindGlobal( "Elms53Case6", function(P)
    local F, e, a, b, c, W, V, U;

    F := GF(P);
    e := [];
    for a in [0..P-1] do
        for b in [0..P-1] do
            c := ((a^2+b^2) mod P);
            if c <> 0 then 
                W := [[a, -b],[b,a]]*One(F);
                V := [[a^2-b^2, -4*a*b],[a*b, a^2-b^2]]*One(F);
                U := (c*W)^-1;
                Add(e, [W, V, U]);

                W := [[a, -b],[-b,-a]]*One(F);
                V := [[a^2-b^2, -4*a*b],[-a*b, -(a^2-b^2)]]*One(F);
                U := (-c*W)^-1;
                Add(e, [W, V, U]);
            fi;
        od; 
    od;
    return e;
end );

BindGlobal( "Elms53Case7", function(P)
    local F, e, W, a, b, c, Z, V, U;

    F := GF(P);
    e := [];
    W := PrimitiveRootMod(P);
    for a in [0..P-1] do
        for b in [0..P-1] do
            c := ((a^2+W*b^2) mod P);
            if c <> 0 then 
                Z := [[a, b],[-W*b,a]]*One(F);
                V := [[a^2-W*b^2, 4*W*a*b],[-a*b, a^2-W*b^2]]*One(F);
                U := (c*Z)^-1;
                Add(e, [Z, V, U]);

                Z := [[a, b],[W*b,-a]]*One(F);
                V := [[a^2-W*b^2, 4*W*a*b],[a*b, -(a^2-W*b^2)]]*One(F);
                U := (-c*Z)^-1;
                Add(e, [Z, V, U]);
            fi;
        od; 
    od;
    return e;
end );

BindGlobal( "Orbits53", function(P, E)
    local F, N, e, t, i, g, h, j, k, c1, c2, c3, B, C, l, r, m, s;

    F := GF(P);
    N := NullMat(2,2,F);

    # case 1 : A2 = 0, A1 arbitrary
    e := Elements(MatrixSpace(F, 2, 2));
    t := List(e, X -> true);
    #Print("got ",Length(e)," elements in case 1 \n");
    for i in [1..Length(e)] do
        if t[i] = true then 
            for g in E do
                h := g[1]*e[i]*g[3];
                j := Position(e, h);
                if j > i then t[j] := false; fi;
            od;
        fi;
    od;
    c1 := e{Filtered([1..Length(e)], X -> t[X]=true)};
    c1 := List(c1, X -> Concatenation(X, N));
    #Print(" -- reduced to ",Length(c1)," orbits \n");

    # case 2 : A2 invertible, A1 = 0
    e := Elements(MatrixSpace(F, 2, 2));
    t := List(e, X -> true);
    #Print("got ",Length(e)," elements in case 2 \n");
    for i in [1..Length(e)] do
        if t[i] = true then 
            t[i] := [];
            for g in E do
                h := g[2]*e[i]*g[3];
                j := Position(e, h);
                if j = i then Add(t[i], g); fi;
                if j > i then t[j] := false; fi;
            od;
        fi;
    od;
    k := Filtered([1..Length(e)], X -> t[X]<>false); e := e{k}; t := t{k};
    c2 := Filtered(e, X -> RankMat(X)=2);
    c2 := List(c2, X -> Concatenation(N, X));
    #Print(" -- determined to ",Length(c2)," orbits \n");

    # case 3 : A2 has rank 1 
    k := Filtered([1..Length(e)], X -> RankMat(e[X])=1); e := e{k}; s := t{k};
    c3 := [];
    for i in [1..Length(e)] do
        m := CanoForms53(e[i],F);
        t := List(m, X -> true);
        #Print("got ",Length(m)," elements in case 3 \n");
        for j in [1..Length(m)] do
            if t[j] = true then 
                for g in s[i] do
                    B := g[1]*m[j]{[1,2]}*g[3];
                    C := Concatenation(B, m[j]{[3,4]});
                    C := CanoForm53(C, F);
                    k := Position(m,C);
                    if IsBool(k) then Error("hier"); fi;
                    if k > j then t[k] := false; fi;
                od;
            fi;
        od;
        m := m{Filtered([1..Length(m)], X -> t[X]=true)};
        #Print(" -- reduced to ",Length(m)," orbits \n");
        Append(c3, m);
    od;

    return Concatenation(c1, c2, c3);
end );

BindGlobal( "ValsFunction19", function(P)
    local E, e, l, r;

    # acting group and orbits
    E := Elms53Case6(P);
    e := Orbits53( P, E );

    # check number
    l := P mod 12;
    if l = 1 then 
        r := 3*P^2+(19/2)*P+15+P^3/2;
    elif l = 5 then 
        r := 3*P^2+(19/2)*P+12+P^3/2;
    elif l = 7 then 
        r := 2*P^2+(5/2)*P+2+P^3/2;
    elif l = 11 then 
        r := 2*P^2+(5/2)*P+3+P^3/2;
    fi;
    if r <> Length(e) then Error("wrong number of orbits"); fi;

    # permute and return
    return List( e, X -> IntVecFFE(Flat(X){[6,7,8,3,1,2,4,5]}));
end );

BindGlobal( "ValsFunction19a", function(P)
    local E, e, l, r;

    # acting group and orbits
    E := Elms53Case7(P);
    e := Orbits53( P, E );

    # check number
    l := P mod 12;
    if l = 1 then 
        r := 2*P^2+(5/2)*P+2+P^3/2;
    elif l = 5 then 
        r := 2*P^2+(5/2)*P+3+P^3/2;
    elif l = 7 then 
        r := 3*P^2+(19/2)*P+13+P^3/2;
    elif l = 11 then 
        r := 3*P^2+(19/2)*P+10+P^3/2;
    fi;
    if r <> Length(e) then Error("wrong number of orbits"); fi;

    # permute and return
    return List( e, X -> IntVecFFE(Flat(X){[6,7,8,3,1,2,4,5]}));
end );

##
## "x ne -1,3, See Notes6.114"
##
BindGlobal( "ValsFunction20", function(P)
    local r, x, m, c, z, t, A, M, B, y, w;
    r := [];
    for x in Difference([0..P-2],[3]) do

        # acting elements
        m := [];
        m[1] := [[x-1,1],[-1,0]] * One(GF(P));
        m[2] := [[x^2-2*x, x-1],[1-x,-1]] * One(GF(P));
        for c in [0..P-2] do
            if ((c*x+c^2-c+1) mod P) <> 0 then 
                Add(m, [[(1+c*x)*(c*x-2*c+1),c*(c*x+2-c)],
                        [-c*(c*x+2-c),-(-1+c)*(c+1)]] * One(GF(P)) );
            fi;
        od;

        # loop
        for z in [0..P-1] do
            t := true;
            A := [[1],[z]] * One(GF(P));

            for M in m do
                if t = true then 
                    B := M*A;
                    y := B[1][1];
                    if y <> 0*y then 
                        w := IntFFE( y^-1 * B[2][1]);
                        if w < z then t := false; fi;
                    fi;
                fi;
            od;

            if t = true then Add(r, [x,z]); fi;
        od;
    od;
    return r;
end );

##
## "See Notes 6.173"
##
BindGlobal( "Mats6173", function(P, W)
    local mats, rang, l, m, n, r, s1, s2, a, b, x, y, z, u, v, w;
    mats := [];
    rang := List([0..P-1], x -> [1,x]); Add(rang, [0,1]);
    for l in [0..P-1] do
        for m in [0..P-1] do
            if (l^2) mod P <> m then 
                n := 1;
                s1 := []; s2 := [];
                for r in rang do
                    if n = 1 then 
                        a:=r[1]; 
                        b:=r[2];
                        x:=(a^2+2*a*b*l+b^2*m) mod P;
                        if x <> 0 then 
                            u := ((W*a*b+a^2*l+W*b^2*l+a*b*m)/x) mod P;
                            v := ((W^2*b^2+2*W*a*b*l+a^2*m)/x) mod P;
                            w := (P-u) mod P;
                            if [u,v] < [l,m] or [w,v] < [l,m] then n := 0; fi;
                            if [u,v] = [l,m] then Add(s1, r); fi;
                            if [w,v] = [l,m] then Add(s2, r); fi;
                        fi;
                    fi;
                od;
                if n=1 then 
                    Add(mats, rec( elm := [l,m], stb1 := s1, stb2 := s2 )); 
                fi;
            fi;
        od;
    od;
    if Length(mats) <> P+1 then Error("wrong number of mats"); fi;
    return mats;
end );

BindGlobal( "ValsFunction21", function(P)
    local W, R, M, T, y, z, t, l, n, A, B, C, D, a, b, m, u, v, r, c;

    # set up
    W := PrimitiveRootMod(P);
    R := [[0,0,0,0,0]];
    M := Mats6173(P, W);
    T := [];
    for y in [0..P-1] do
        for z in [0..P-1] do
            if y+z > 0 then 
               l := [0];
            else
               l := [0..P-1];
            fi;
            for t in l do
                if [y,z,t] <> [0,0,0] then Add( T, [y,z,t] ); fi;
            od;
        od;
    od;
   
    # case 1 : u = v = 0
    for y in [0..P-1] do
        for z in [0..P-1] do
            for t in [0..P-1] do
                A := [y,z,t];
                if (t = 0 or z = 0) and A <> [0,0,0] then 
                    n := 1;
                    for a in [1..P-1] do
                        for b in [1,-1] do
                            if n = 1 then 
                                B := [[a,0,0],[0,a*b,0],[0,0,a^2]] *One(GF(P));
                                C := (a^4*b) *One(GF(P));
                                D := IntVecFFE(C^-1 * A * B);
                                if D < A then n := 0; fi;
                            fi;
                        od;
                    od;
                    if n = 1 then Add(R, Concatenation(A, [0,0])); fi;
                fi;
            od;
        od;
    od;

    l := (P+1+(P+3)*Gcd(P-1,3)+Gcd(P-1,4))/2;
    if Length(R) <> l then Error("wrong number in case 1"); fi;

    # case 2 : x = y = z = 0
    for m in M do
        Add(R, Concatenation([0,0,0],m.elm));
    od;

    # case 3 : other
    for m in M do
        u := m.elm[1];
        v := m.elm[2];
        for A in T do 
            n := 1;
            for r in m.stb1 do
                a := r[1]; 
                b := r[2];
                for c in [1..P-1] do
                    if n = 1 then 
                        B := [[a, W*b, 0],
                              [b, a, 0],
                              [0,0,(a^2-W*b^2)*c]] * One(GF(P));
                        C := ((a^2-W*b^2)*(a^2+2*a*b*u+b^2*v)*c^3)
                              *One(GF(P));
                        D := IntVecFFE( C^-1 * A * B );
                        if D < A then n := 0; fi;
                    fi;
                od;
            od;
                
            for r in m.stb2 do
                a := r[1]; 
                b := r[2];
                for c in [1..P-1] do
                    if n = 1 then     
                        B := [[a, -W*b, 0],
                              [b, -a, 0],
                              [0,0,(a^2-W*b^2)*c]] *One(GF(P));
                        C := (-(a^2-W*b^2)*(a^2+2*a*b*u+b^2*v)*c^3)
                             * One(GF(P));
                        D := IntVecFFE(C^-1*A*B);
                        if D < A then n := 0; fi;
                    fi;
                od;
            od;

            if n = 1 then Add(R, Concatenation(A, m.elm)); fi;
        od;
    od;
             
    l := 3*P + 3 + (P^2+2*P+3)*Gcd(P-1,3)/2;
    if l <> Length(R) then Error("wrong number in case 2"); fi;

    return R;
end );
 

##
##  "See notes5.12"
##
BindGlobal( "Subgroup512Case1", function( P )
    local U, a, b, m, l;
    U := Subgroup(GL(2,P), []);
    for a in [0..P-1] do
        for b in [0..P-1] do
            if a <> b and a <> ((P-b) mod P) then 
                m := [[a,b],[b,a]] * One(GF(P));
                l := [[a,b],[-b,-a]] *One(GF(P));
                U := ClosureGroup(U, [m, l]);
                if Size(U) = 2*(P-1)^2 then return U; fi;
            fi;
        od;
    od;
    return U;
end );

BindGlobal( "Subgroup512Case2", function( P )
    local U, W, a, b, m, l;
    U := Subgroup(GL(2,P), []);
    W := PrimitiveRootMod(P);
    for a in [0..P-1] do
        for b in [0..P-1] do
            if (a^2-W*b^2) mod P <> 0 then 
                m := [[a,W*b],[b,a]] * One(GF(P));
                l := [[a,W*b],[-b,-a]] *One(GF(P));
                U := ClosureGroup(U, [m, l]);
                #if Size(U) = 2*(P-1)*(P+1) then return U; fi;
            fi;
        od;
    od;
    return U;
end );

BindGlobal( "Orbits512", function(P, case)
    local U, A, f, o, l;
    if case = 1 then 
        U := Subgroup512Case1(P);
    else
        U := Subgroup512Case2(P);
    fi;
    A := MatrixSpace(GF(P), 2, 2);
    f := function( elm, mat )
        return (mat*elm*mat^-1)/Determinant(mat);
    end;
    o := Orbits(U, Elements(A), f);
    if case = 1 then 
        l := P^2 + (7*P+15)/2;
        if Length(o) <> l then Error("wrong number of orbits"); fi;
    else
        l := P^2 + (3*P+5)/2;
        if Length(o) <> l then Error("wrong number of orbits"); fi;
    fi;
    return List(o, X -> IntVecFFE(Flat(X[1])));
end );

BindGlobal( "ValsFunction22", function(P)
    return Orbits512(P, 1);
end );

BindGlobal( "ValsFunction22a", function(P)
    return Orbits512(P, 2);
end );

##
## "See Notes6.150
##
BindGlobal( "Less6150", function( P, A, D )
    return P^2*A[3]+P*A[1]+A[2] > P^2*D[3]+P*D[1]+D[2];
end );

# case 3
BindGlobal( "ValsFunction23", function(P)
    local m, k, x, l, y, z, A, n, a, c, B, C, D;
    m := [[0,0,0]];
    k := LeastNonSquareModP(P);
    for x in [0..P-1] do
        for y in [0..P-1] do
            l := [0,1,k]; if x > 0 then l := [0]; fi;
            for z in l do
                A := [x,y,z];
                if A <> [0,0,0] then 
                    n := 1;
                    for a in [1..P-1] do
                        for c in [0..P-1] do
                            if n = 1 then 
                                B := [[a,0,0],[0,a,0],[c,-c,a^2]] * One(GF(P));
                                C := (a^4) *One(GF(P));
                                D := IntVecFFE(C^-1*A*B); 
                                if Less6150(P, A, D) then n := 0; fi;
                                B := [[0,a,0],[a,0,0],[c,-c,a^2]] * One(GF(P));
                                C := (-(a^4)) *One(GF(P));
                                D := IntVecFFE(C^-1*A*B); 
                                if Less6150(P, A, D) then n := 0; fi;
                            fi;
                        od;
                    od;
                    if n = 1 then Add(m, A); fi;
                fi;
            od;
        od;
    od;
    l := (P+1+(P+3)*Gcd(P-1,3)+Gcd(P-1,4))/2;
    if Length(m) <> l then Error("wrong number"); fi;
    return m;
end );

# case 4
BindGlobal( "ValsFunction23a", function(P)
    local m, k, x, l, y, z, A, n, a, b, c, B, C, D;
    m := [[0,0,0]];
    k := LeastNonSquareModP(P);
    for x in [0..P-1] do
        for y in [0..P-1] do
            l := [0,1,k]; if x+y > 0 then l := [0]; fi;
            for z in l do
                A := [x,y,z];
                if A <> [0,0,0] then
                    n := 1;
                    for a in [1..P-1] do
                        for b in [1,-1] do
                            if n = 1 then
                                B := [[a*b,0,0],[0,a,0],[0,0,a^2]] * One(GF(P));
                                C := (a^4) *One(GF(P));
                                D := IntVecFFE(C^-1*A*B);
                                if D < A then n := 0; fi;
                                B := [[0,a,0],[a*b,0,0],[0,0,a^2]] * One(GF(P));
                                C := (a^4) *One(GF(P));
                                D := IntVecFFE(C^-1*A*B);
                                if D < A then n := 0; fi;
                            fi;
                        od;
                    od;
                    if n = 1 then Add(m, A); fi;
                fi;
            od;
        od;
    od;
    l := 3+Gcd(P-1,3)*(P+3+Gcd(P-1,4))/4;
    if Length(m) <> l then Error("wrong number"); fi;
    return m;
end );

# case 5
BindGlobal( "ValsFunction23b", function(P)
    local W, m, k, x, l, y, z, A, n, a, b, c, B, C, D;
    W := PrimitiveRootMod(P);
    m := [[0,0,0]];
    for x in [0..P-1] do
        for y in [0..P-1] do
            l := [0,1]; if x+y > 0 then l := [0]; fi;
            for z in l do
                A := [x,y,z];
                if A <> [0,0,0] then
                    n := 1;
                    for a in [1..P-1] do
                        for b in [1,-1] do
                            if n = 1 then
                                B := [[a*b,0,0],[0,a,0],[0,0,a^2]] * One(GF(P));
                                C := (a^4) *One(GF(P));
                                D := IntVecFFE(C^-1*A*B);
                                if D < A then n := 0; fi;
                                B := [[0,W*a,0],[a*b,0,0],[0,0,a^2]]*One(GF(P));
                                C := (W^2*a^4) *One(GF(P));
                                D := IntVecFFE(C^-1*A*B);
                                if D < A then n := 0; fi;
                            fi;
                        od;
                    od;
                    if n = 1 then Add(m, A); fi;
                fi;
            od;
        od;
    od;
    l := 2+Gcd(P-1,3)*(P+7-Gcd(P-1,4))/4;
    if Length(m) <> l then Error("wrong number"); fi;
    return m;
end );

# case 8
BindGlobal( "ValsFunction23c", function(P)
    local o, m, k, x, l, y, z, t, A, n, a, c, B, C, D;
    o := One(GF(P));
    k := LeastNonSquareModP(P);
    m := [];
    for x in [2..P-1] do
        Add(m, [x,0,0,0]);
        for y in [0..P-1] do
            for z in [0..P-1] do
                l := [0,1,k]; if y+z > 0 then l := [0]; fi;
                for t in l do
                    A := [y,z,t];
                    if A <> [0,0,0] then
                        n := 1;
                        for a in [1..P-1] do
                            if n = 1 then
                                B := [[a,0,0],[0,a,0],[0,0,a^2]]*o;
                                C := (a^4) * o;
                                D := IntVecFFE(C^-1*A*B);
                                if D < A then n := 0; fi;
                                B := [[0,x*a,0],[a,0,0],[0,0,-x*a^2]]*o;
                                C := (x^2*a^4) * o;
                                D := IntVecFFE(C^-1*A*B);
                                if D < A then n := 0; fi;
                            fi;
                        od;
                        if n = 1 then Add(m, [x,y,z,t]); fi;
                    fi;
                od;
            od;
        od;
    od;
    l := (5*P-7+(P^2-5)*Gcd(P-1,3)-Gcd(P-1,4))/2;
    if Length(m) <> l then Error("wrong number"); fi;
    return m;
end );

##
## "See Notes6.178
##
BindGlobal( "Mats6178", function(P)
    local W, mats, y1, y2, y3, y4, A, B, n, a, b, c, Q, C, D, E, l;
    W := PrimitiveRootMod(P);
    mats:=[];

    for y1 in [0,1] do
        for y2 in [0..P-1] do
            for y3 in [0..P-1] do
                y4:=(P-y1) mod P;
                A := [[y1,y2],[y3,y4]];
                B := [y1,y2,y3];
                if RankMat(A) = 2 then 
                    n := 1;
                    for a in [0..P-1] do
                        for b in [0..P-1] do
                            if a+b > 0 then 
                                for c in [1,-1] do
                                    Q := [[a,b],[W*b*c, a*c]];
                                    C := (c*(a^2-W*b^2)) * One(GF(P));
                                    D := Q * A * Q^-1 * C^-1;
                                    E := IntVecFFE([D[1][1],D[1][2],D[2][1]]);
                                    if E < B then n := 0; fi;
                                od;
                            fi;
                        od;
                    od;
                    if n = 1 then Add(mats, A); fi;
                fi;
            od;
        od;
    od;
    l := (3*P-1)/2;
    if Length(mats) <> l then Error("wrong number"); fi;
    return mats;
end );

BindGlobal( "ValsFunction24", function(P)
    local W, S, mats, A, x, y, z, r, s, u, t, u1, val, res, res1, res2, l;

    W := PrimitiveRootMod(P);
    S := Set(List([0..(P-1)/2], x -> (x^2) mod P ));
    res := [];
    res1 := [];
    res2 := [];
    mats := Mats6178(P);
    for A in mats do

        x := A[1][1];
        y := A[1][2];
        z := A[2][1];
        r := (x+y)/2 mod P;
        s := (z-x)/2 mod P;
        u := (x^2+y*z)*(2*(W*y+z)^2+W*(x^2+y*z)) mod P;

        # case 4
        if P mod 4 <> 1 and x > 0 and y > 0 and (W*y+z) mod P = 0 then 
            for t in [0..(P-1)/2] do
                Add(res, [t, y]);
            od;
        fi;


        # case 5
        if x = 0 and (z-(W*y)) mod P <> 0 and (z-(W*(P-y))) mod P <> 0 then 

            for t in [0..(P-1)/2] do
                Add(res1, [y,z,0,0,t]);
            od;

            if u = 0 then 
                if (W*y+2*z) mod P = 0 then
                    for t in [1..(P-1)/2] do
                        Add(res1, [y,z,t,0,0]);
                    od;
                fi;
                if (2*W*y+z) mod P = 0 then
                    for t in [1..(P-1)/2] do
                        Add( res1, [y,z,0,t,0]);
                    od;
                fi;
            fi;

            if u <> 0 and u in S then
                u1:= u * y^-2 mod P;
                val := First([1..(P-1)/2], x -> (x^2) mod P = u1);
                for t in [1..(P-1)/2] do
                    Add( res1, [y,z,0,t,val]);
                od;
            fi;
        fi;

        # case 6
        if x = 1 and (z+W*y) mod P <> 0 then
            for t in [0..P-1] do
                Add( res2, [y,z,0,0,t]);
            od;
            if u in S then
                val := First([1..P-1], x -> (x^2) mod P = u); 
                if (-W*y-z+val) mod P = 0 and (-W*y^2-2*y*z-1) mod P = 0 then 
                    for t in [1..(P-1)/2] do
                        Add( res2, [y,z,t,0,val]);
                    od;
                else
                    for t in [1..(P-1)/2] do
                        Add( res2, [y,z,0,t,val]);
                    od;
                fi;
                if (-W*y-z-val) mod P = 0 and (-W*y^2-2*y*z-1) mod P = 0 then 
                    for t in [1..(P-1)/2] do
                        Add( res2, ([y,z,t,0,-val] mod P));
                    od;
                else
                    for t in [1..(P-1)/2] do
                        Add( res2, ([y,z,0,t,-val] mod P));
                    od;
                fi;
            fi;
        fi;
    od;

    if P mod 4 = 3 then 
        l := (P+1)/2;
        if Length(res) <> l then Error("wrong number case 4"); fi;
    fi;

    l := (3*P^2 - 3*P -4)/2;
    if Length(res1)+Length(res2) <> l then Error("wrong number case 5/6"); fi;

    return [res, res1, res2];
end );

##
## "See note2dec5.1"
##
BindGlobal( "Range51", function(P)
    local SQ, ns, t, x, y, z, delta, zend, r;
    SQ := SquaresModP(P);
    ns := LeastNonSquareModP(P);
    r := [];
    for t in [1,ns] do
        for x in [0..(P-1)/2] do
            for y in [1..P-1] do
                if x = 0 then 
                    zend:=(P-1)/2; 
                else
                    zend:= (P-1); 
                fi;
                for z in [0..zend] do
                    delta := ((t*z-y*x)^2-t*y) mod P;
                    if not delta in SQ then 
                        Add(r, [t,x,y,z]);
                    fi;
                od;
            od;
        od;
    od;
    return r;
end );
    
BindGlobal( "ValsFunction25", function(P)
    local mats, F, l, A, t, x, y, z, t1, x1, z1, y1, new, d, test1, test2,
          brange, arange, a, b, e, k, f, t2, x2, y2, z2;

    # set up
    mats := [];
    F := GF(P);

    # expected number 
    l := (P^2-1)/2; if (P mod 3) = 2 then l := l+1; fi;

    # loop
    for A in Range51(P) do
        t := A[1]; x := A[2]; y := A[3]; z := A[4];
        t1 := t*One(F); x1 := x*One(F); y1 := y*One(F); z1 := z*One(F);
        new:=1;

        for d in Elements(F) do
            test1 := 2*z1-2*d*y1-d;
            test2 := 2*t1*d^2-One(F)-2*x1*d;
            if new = 0 then 
                brange := [];
            elif test1=Zero(F) and test2<>Zero(F) then 
                brange := [0]; 
            else 
                brange := Elements(F);
            fi;
            for b in brange do
                if new = 0 then 
                   arange := [];
                elif test1<>Zero(F) then 
                   arange := [-b*test2/test1]; 
                else
                   arange := Elements(F);
                fi;
                for a in arange do
                    e:=a*d-b;
                    if e<>Zero(F) and new = 1 then 
                        k := 2*(b*d*t1-a*y1)*e^-1;
                        f := e^-2*k^-1;
                        t2 := IntFFE(f*(d^3*t1-d*y1-d-d^2*x1+z1));
                        if t2 < t then new:=0; fi;
                        x2 := IntFFE(f*(-b*d^2*t1+b*y1+a*d+a*d^2*x1-a*z1));
                        if t2=t and x2 < x then new:=0; fi;
                        y2 := IntFFE(f*(-b^2*d*t1+d*a^2*y1+a*b+b^2*x1-a^2*z1));
                        if t2=t and x2=x and y2 < y then new := 0; fi; 
                        z2 := IntFFE(f*(b^3*t1-b*a^2*y1-a^2*b-a*b^2*x1+a^3*z1));
                        if t2=t and x2=x and y2=y and z2 < z then new := 0; fi;
                    fi;
                od;
            od;
        od;
        if new = 1 then Add(mats, [x,y,z,t]); fi;
        if Length(mats) = l then return mats; fi;
    od;
end );

##
## "See Notes6.163a/b"
##
BindGlobal( "ValsFunction26", function(P)
    local SQ, lns, mats1, mats, u, x, n, A, u1, x1, l, yrange, zrange, 
          t, y, z, new, a, c, e, y1, z1, t1;  

    SQ := Set(List([1..(P-1)/2], x -> (x^2) mod P));
    lns := LeastNonSquareModP(P);
    mats1 := [];

    mats := [];
    for u in [1,lns] do
        for x in [1..P-1] do
            n := 1;
            A := [u,x];
            if x in SQ then 
                u1:=1;
                x1:=(u*x) mod P;
            else
                u1:=lns;
                x1:=(u*x/lns) mod P;
            fi;
            if [u1,x1] >= [u,x] then Add( mats, [u,x] ); fi;
        od;
    od;
    l := 3*(P-1)/2;
    if Length(mats) <> l then Error("wrong number of mats"); fi;

    for A in mats do
        u:=A[1]; x:=A[2];

     
        if u*x mod P <> 1 then 
            yrange:=[0]; 
            zrange:=[0..(P-1)/2]; 
        else
            yrange:=[0..(P-1)/2];
            zrange:=[0..P-1];
        fi;
        for t in [0..P-1] do
            for y in yrange do
                for z in zrange do
                    new:=1;

                    for a in [1,P-1] do
                        for c in [0..P-1] do
                            if new = 1 then 
                                if (u*x) mod P = 1 then 
                                    e:=(c*u^-1) mod P;
                                else 
                                    e:=(-(c*u^-1+c*t)) mod P;
                                fi;
                                y1:=(e+c*u^-1+a*y+c*t) mod P;
                                z1:=(x*e-x*c*u^-1-2*e*u^-1+z*a+e*t) mod P;
                                if [t,y1,z1]<[t,y,z] then new := 0; fi;
                            fi;
                        od;
                    od;
                        
                    for a in [1..P-1] do
                        if u = ((a^2*x) mod P) then
                            for c in [0..P-1] do
                                if new = 1 then 
                                    if (u*x) mod P = 1 then 
                                        e:=(c*u^-1) mod P;
                                    else
                                        e:=(x*c-2*c*u^-1+z*a+c*t) mod P;
                                    fi;
                                    y1:=(-(x*c-e-2*c*u^-1+z*a+c*t)) mod P;
                                    z1:=(-(x*c*u^-1+e*u^-1+a^-1*y+e*t)) mod P;
                                    t1:=(-t+u^-1-x) mod P;
                                    if [t1,y1,z1] < [t,y,z] then new := 0; fi;
                                fi;
                            od;
                        fi;
                    od;
                    if new = 1 then Add( mats1, [x,y,z,t,u]); fi;
                od;
            od;
        od;
    od;
    l := 2*P^2-(5*P-1)/2;
    if l <> Length(mats1) then Error("wrong number"); fi;
    return mats1;
end );

BindGlobal( "ValsFunction26a", function(P)
    local SQ, lns, mats2, mats, u, x, n, A, u1, x1, l, yrange, zrange,
          t, z, y, new, a, c, e, y1, z1, t1;

    SQ := Set(List([1..(P-1)/2], x -> (x^2) mod P));
    lns := LeastNonSquareModP(P);
    mats2 := [];

    mats := [];
    for u in [1,lns] do
        for x in [1..P-1] do
            if (u*x) mod P <> 1 then 
                n := 1;
                A := [u,x];
                if x in SQ then
                    u1:=1;
                    x1:=(u^-1*x^-1) mod P;
                else
                    u1:=lns;
                    x1:=(u^-1*x^-1/lns) mod P;
                fi;
                if [u1,x1] >= [u,x] then Add( mats, [u,x] ); fi;
            fi;
        od;
    od;
    l := P-3+Gcd(P-1,4)/2;
    if Length(mats) <> l then Error("wrong number of mats"); fi;

    for A in mats do
        u:=A[1]; x:=A[2];

        for t in [0..P-1] do
            for y in [0..(P-1)/2] do
                for z in [0..P-1] do
                    new:=1;
                    for a in [1,-1] do
                        for e in [0..P-1] do
                            if new = 1 then 
                                y1:=(2*e+a*y+u*e*t) mod P;
                                z1:=(-2*x*e+a*z+e*t) mod P;
                                if [y1,z1]<[y,z] then new := 0; fi;
                            fi;
                        od;
                    od;
                    if new=1 and ((u*x+1) mod P)=0 and (P mod 4)=1 then
                        for a in [2..P-2] do
                            if (a^2+1) mod P = 0 then 
                                for e in [0..P-1] do
                                    if new = 1 then
                                        y1:=(-(2*e+z*a*u+u*e*t)) mod P;
                                        z1:=(x*(2*e+a*y+u*e*t)) mod P;
                                        if [y1,z1] < [y,z] then new := 0; fi;
                                    fi;
                                od;
                            fi;
                        od;
                    fi;

                    if new = 1 then Add(mats2, [x,y,z,t,u]); fi;
                od;
            od;
        od;
    od;
    l := (P^3-5*P+P*Gcd(P-1,4))/2;
    if Length(mats2) <> l then Error("wrong number"); fi;
    return mats2;
end );

##
## "See Notes4.1"
##
BindGlobal( "ValsFunction27", function(P)
    local F, A, B, W, range, mats, s, x, y, z, t, AA, new, r, a, b, C, u, 
          x1, y1, z1, t1, l;

    F := GF(P);
    A := NullMat(2, 2, F);
    B := NullMat(2, 2, F);
    W := PrimitiveRootMod(P);
    range := Concatenation( [[0,1]], List([0..P-1], x -> [1,x]));

    mats := [];
    for s in range do
        x:=s[1]; 
        y:=s[2];
        A[1][1]:= x * One(F);
        A[1][2]:= y * One(F);
        for z in [0..P-1] do
            A[2][1]:= z * One(F);
            for t in [0..P-1] do
                A[2][2]:= t * One(F);
                if RankMat(A) = 2 then
                    AA := [x,y,z,t];
                    new:=1;

                    for r in range do
                        if new = 1 then 
                            a:=r[1]; 
                            b:=r[2];

                            B[1][1]:= a * One(F);
                            B[2][2]:= a * One(F);
                            B[1][2]:= b * One(F);
                            B[2][1]:= (W*b) * One(F);
                            C:=B*A*B^-1;
                            u:=C[1][1];
                            if u = 0*u then u:=C[1][2]; fi;
                            C:=u^-1*C;
    
                            x1:= IntFFE(C[1][1]);
                            y1:= IntFFE(C[1][2]);
                            z1:= IntFFE(C[2][1]);
                            t1:= IntFFE(C[2][2]);
                            if [x1,y1,z1,t1] < AA then new := 0; fi;
    
                            B[2][1]:=-B[2][1];
                            B[2][2]:=-B[2][2];
                            C:=B*A*B^-1;
    
                            u:=C[1][1];
                            if u = 0*u then u:=C[1][2]; fi;
                            C:=u^-1*C;
    
                            x1:= IntFFE(C[1][1]);
                            y1:= IntFFE(C[1][2]);
                            z1:= IntFFE(C[2][1]);
                            t1:= IntFFE(C[2][2]);
                            if [x1,y1,z1,t1] < AA then new := 0; fi;
                        fi;
                    od;

                    if new = 1 then Add(mats, AA ); fi;
                fi;
            od;
        od;
    od;
    l := (P+1)^2/2;
    if Length(mats) <> l then Error("wrong number"); fi;
    return mats;
end );

