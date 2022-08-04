
##
## Specialize Lie ring with given values
##
## P : a prime or false
## para : a list of parameters (not p or w)
## vals : the values of the parameters
##
BindGlobal( "SpecialiseLiePRingNC", function( L, P, para, vals )
    local S, T, W, p, w, sp, sx, i, j, K;

    # get table
    S := SCTable(Zero(L));

    # check for cases
    sp := (IsInt(P) and not IsInt(S.prime));
    sx := (IsBound(S.param) and Length(Intersection(S.param,para))>0);
    if sp = false and sx = false then return L; fi;
    
    # do some work
    T := StructuralCopy(S);
    Unbind(T.fam);

    # specialize parameters
    if sp then 
        T.prime := P;
        W := PrimitiveRootMod(P);
        p := IndeterminateByName("p");
        w := IndeterminateByName("w");
        para := Concatenation( [p,w], para );
        vals := Concatenation( [P,W], vals );
    fi;

    if not IsBound(T.param) then 
        T.param := []; 
    else 
        T.param := Difference(T.param, para);
    fi;
    
    # evaluate
    for i in [1..Length(T.tab)] do
        for j in [1..Length(T.tab[i])] do
            if IsPolynomial(T.tab[i][j]) then 
                T.tab[i][j] := Value( T.tab[i][j], para, vals );
            fi;
            if IsRationalFunction(T.tab[i][j]) then 
                T.tab[i][j] := ValueRatFun( T.tab[i][j], para, vals );
            fi;
        od;
    od; 

    # get entries in normal form
    for i in Reversed([1..Length(T.tab)]) do
        T.tab[i] := WordByExps@(LRReduceExp(T,ExpsByWord(T.dim,T.tab[i])));
    od;

    # create new Lie ring
    K := LiePRingBySCTableNC(T);
    SetIsLiePRing(K, true);

    # transfer data
    if IsBound(L!.ShortPresentation) then 
        SetShortPresentation(K, L!.ShortPresentation); fi;
    if IsBound(L!.LibraryConditions) then 
        SetLibraryConditions(K, L!.LibraryConditions); fi;
    if IsBound(L!.LibraryName) then 
        SetLibraryName(K, L!.LibraryName); fi;
    if HasDimensionOfLiePRing(L) then 
        SetDimensionOfLiePRing(K, DimensionOfLiePRing(L) ); fi;
    if HasClassOfLiePRing(L) then 
        SetClassOfLiePRing(K, ClassOfLiePRing(L) ); fi;
    if HasPClassOfLiePRing(L) then 
        SetPClassOfLiePRing(K, PClassOfLiePRing(L) ); fi;
    if HasMinimalGeneratorNumberOfLiePRing(L) then 
        SetMinimalGeneratorNumberOfLiePRing(K, 
                    MinimalGeneratorNumberOfLiePRing(L) ); fi;
    SetLiePValues(K, [para,vals]);
    return K;
end );

##
## check if P satisfies conditions of L and if cl(L) < P
##
BindGlobal( "IsValidPrime", function(L, P)
    local c;
    if not IsPrimeInt(P) then return false; fi;
    if not HasLibraryConditions(L) then return true; fi;
    if P=2 or PClassOfLiePRing(L) > P-1 then return false; fi;
    c := L!.LibraryConditions[2];
    if c = "p=1 mod 3" then  
        if (P mod 3) <> 1 then return false; fi;
    elif c = "p=2 mod 3" then 
        if (P mod 3) <> 2 then return false; fi;
    elif c = "p=1 mod 4" then 
        if (P mod 4) <> 1 then return false; fi;
    elif c = "p=3 mod 4" then 
        if (P mod 4) <> 3 then return false; fi;
    elif c = "p=1 mod 5" then 
        if (P mod 5) <> 1 then return false; fi;
    elif c = "p ne 1 mod 5" then 
        if (P mod 5) = 1 then return false; fi;
    elif c = "p=1 mod 7" then 
        if (P mod 7) <> 1 then return false; fi;
    elif c = "p=1 mod 8" then 
        if (P mod 8) <> 1 then return false; fi;
    elif c = "p=1 mod 9" then 
        if (P mod 9) <> 1 then return false; fi;
    elif c = "p=1 mod 4, p=1 mod 3" then
        if (P mod 4) <> 1 or (P mod 3) <> 1 then return false; fi;
    elif c = "p=3 mod 4, p=1 mod 3" then
        if (P mod 4) <> 3 or (P mod 3) <> 1 then return false; fi;
    elif c = "p=±1 mod 8" then 
        if not (P mod 8) in [1,7] then return false; fi;
    elif c = "p=±3 mod 8" then 
        if not (P mod 8) in [3,5] then return false; fi;
    fi;
    return true;
end );

##
## Specialize with checking 
##
BindGlobal( "SpecialiseLiePRing", function( L, P, para, vals )
    local S, p, w;

    # get table
    S := SCTable(Zero(L));
    p := IndeterminateByName("p");
    w := IndeterminateByName("w");

    # check for consistency 
    if p in para or w in para then 
        Print("don't include p or w in para \n");
        return fail;
    elif IsInt(P) and IsInt(S.prime) and P <> S.prime then 
        Print("cannot specialize prime to ",P," \n");
        return fail;
    elif IsInt(P) and not IsValidPrime(L,P) then 
        Print("prime does not satisfy condition \n");
        return fail;
    elif Length(para) <> Length(vals) then 
        Print("parameters and values don't match \n");
        return fail;
    fi;
   
    return SpecialiseLiePRingNC( L, P, para, vals );
end );

##
## Specialize prime only
##
BindGlobal( "SpecialisePrimeOfLiePRing", function( L, P )
    return SpecialiseLiePRing( L, P, [], [] );
end );
BindGlobal( "SpecialisePrimeOfLiePRingNC", function( L, P )
    return SpecialiseLiePRingNC( L, P, [], [] );
end );

##
##
##
BindGlobal( "TranslatedLiePRings", function( L, P, para, vals, flag )
    local res, i, R;
    res := List(vals, x -> true);
    for i in [1..Length(vals)] do
        R := SpecialiseLiePRingNC(L, P, para, vals[i]);
        if flag = "group" or flag = "code" then 
            R := PGroupByLiePRing(R);
            if flag = "code" then 
                R := CodePcGroup(R);
            fi;
        fi;
        res[i] := R;
    od; 
    return res;
end );

## 
## Lie rings in a family
##

BindGlobal( "LiePRingsInFamily", function( arg )
    local L, P, W, p, w, d, c, l, para, flag, vals, fals, X, Y, Z, 
          x, y, z, t, j, k, m, n, r, s, u, v;

    # catch arguments
    L := arg[1];
    P := arg[2];
    W := PrimitiveRootMod(P);
    if IsBound(arg[3]) then flag := arg[3]; else flag := false; fi;

    # get Lie ring
    p := PrimeOfLiePRing(L);
    w := IndeterminateByName("w");
    d := DimensionOfLiePRing(L);
    c := LibraryConditions(L);
    l := LibraryName(L);
    para := ParametersOfLiePRing(L);
    vals := [];
    fals := 0;

    # check
    if IsInt(p) then return fail; fi;
    if not IsValidPrime(L,P) then return fail; fi;
    if not IsInt(P) then return fail; fi;

    # the 0-parameters case
    if Length(para) = 0 then 
        return TranslatedLiePRings( L, P, [], [[]], flag );
    fi;

    # the case of trivial conditions
    if Length(c[1]) = 0 then 
        vals := Tuples([0..P-1], Length(para));
        return TranslatedLiePRings( L, P, para, vals, flag );
    fi;

    # assign variable names
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    t := IndeterminateByName("t");
    j := IndeterminateByName("j");
    k := IndeterminateByName("k");
    m := IndeterminateByName("m");
    n := IndeterminateByName("n");
    r := IndeterminateByName("r");
    s := IndeterminateByName("s");
    u := IndeterminateByName("u");
    v := IndeterminateByName("v");

    # special cases with notes - include para to have better overview
    if l = "6.178" then # See note6.178 
        para := [x,y,z,t];
        vals := ValsFunction8(P);
    elif l = "6.62" then # See note6.62 
        para := [x,y];
        vals := ValsFunction9(P);
    elif l = "7.62" then # See note5.38 
        para := [x,y,z,t];
        vals := ValsFunction12(P, 1);
    elif l = "7.63" then # See note5.38 
        para := [x,y,z,t];
        vals := ValsFunction12(P, 2);
    elif l = "7.729" then # See Notes5.12 
        para := [x,y,z,t];
        vals := ValsFunction22(P);
    elif l = "7.730" then # See Notes5.12 
        para := [x,y,z,t];
        vals := ValsFunction22a(P);
    elif l = "7.757" then # See Notes5.14, Case 1 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 1);
    elif l = "7.758" then # See Notes5.14, Case 2 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 2);
    elif l = "7.759" then # See Notes5.14, Case 3 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 3);
    elif l = "7.760" then # See Notes5.14, Case 4 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 4);
    elif l = "7.761" then # See Notes5.14, Case 5 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 5);
    elif l = "7.762" then # See Notes5.14, Case 6 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 6);
    elif l = "7.763" then # See Notes5.14, Case 7 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 7);
    elif l = "7.764" then # See Notes5.14, Case 8 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 8);
    elif l = "7.765" then # See Notes5.14, Case 9 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 9);
    elif l = "7.766" then # See Notes5.14, Case 10 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 10);
    elif l = "7.767" then # See Notes5.14, Case 11 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 11);
    elif l = "7.768" then # See Notes5.14, Case 12 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 12);
    elif l = "7.769" then # See Notes5.14, Case 13 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 13);
    elif l = "7.770" then # See Notes5.14, Case 14 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 14);
    elif l = "7.771" then # See Notes5.14, Case 15 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 15);
    elif l = "7.772" then # See Notes5.14, Case 16 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 16);
    elif l = "7.773" then # See Notes5.14, Case 17 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 17);
    elif l = "7.774" then # See Notes5.14, Case 18 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 18);
    elif l = "7.775" then # See Notes5.14, Case 19 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 19);
    elif l = "7.776" then # See Notes5.14, Case 20 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 20);
    elif l = "7.777" then # See Notes5.14, Case 21 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 21);
    elif l = "7.778" then # See Notes5.14, Case 22 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 22);
    elif l = "7.779" then # See Notes5.14, Case 23 
        para := [x,y,z,t,u,v];
        vals := ValsFunction28(P, 23);
    elif l = "7.780" then # See Notes5.14, Case 24 
        para := [x,y,z,t,s,u,v];  
        vals := ValsFunction28(P, 24);
    elif l = "7.1691" then # See Notes6.150, Case 3 
        para := [x,y,z];
        vals := ValsFunction23(P);
    elif l = "7.1692" then # See Notes6.150, Case 4 
        para := [x,y,z];
        vals := ValsFunction23a(P);
    elif l = "7.1693" then # See Notes6.150, Case 5 
        para := [x,y,z];
        vals := ValsFunction23b(P);
    elif l = "7.1709" then # See Notes6.150, Case 8 
        para := [x,y,z,t];
        vals := ValsFunction23c(P);
    elif l = "7.1763" then # See Notes6.163a 
        para := [x,y,z,t,u];
        vals := ValsFunction26(P);
    elif l = "7.1764" then # See Notes6.163b 
        para := [x,y,z,t,u];
        vals := ValsFunction26a(P);
    elif l = "7.1777" then # See Notes6.173 
        para := [y,z,t,u,v];
        vals := ValsFunction21(P);
    elif l = "7.1797" then # See Notes6.178 
        para := [x,y];
        if fals = 0 then fals := ValsFunction24(P); fi; vals := fals[1];
    elif l = "7.1798" then # See Notes6.178a 
        para := [y,z,t,u,v];
        if fals = 0 then fals := ValsFunction24(P); fi; vals := fals[2];
    elif l = "7.1799" then # See Notes6.178b 
        para := [y,z,t,u,v];
        if fals = 0 then fals := ValsFunction24(P); fi; vals := fals[3];
    elif l = "7.3068" then # See Notes4.1 Case 4 
        para := [x,y,z,t];
        vals := ValsFunction27(P);
    elif l = "7.3285" then # See Notes4.1, Case 5 
        para := [x,y,z,t,j,k,m,n,r,s,u,v];
        vals := ValsFunction27a(P);
    elif l = "7.3286" then # See Notes4.1, Case 6 
        para := [x,y,z,t,j,k,m,n,r,s,u,v];
        vals := ValsFunction27b(P);
    elif l = "7.3387" then # See Notes5.3, Case 4 
        para := [x,y,z,t];
        vals := ValsFunction8a(P);
    elif l = "7.3437" then # See Notes5.3, Case 6 
        para := [x,y,z,t,r,s,u,v];
        vals := ValsFunction19(P);
    elif l = "7.3438" then # See Notes5.3, Case 7 
        para := [x,y,z,t,r,s,u,v];
        vals := ValsFunction19a(P);
    elif l = "7.4670" then # See note1dec5.1 
        para := [x,y];
        vals := ValsFunction9(P);
    elif l = "7.4723" then # See note2dec5.1 
        para := [x,y,z,t];
        vals := ValsFunction25(P);
    elif l = "7.1389" then # See Notes6.114
        para := [x,y];
        vals := ValsFunction20(P);
    fi;
    if vals <> [] then return TranslatedLiePRings( L, P, para, vals, flag ); fi;

    c := c[1];

    # 1 parameter
    if Length(para) = 1 then 

        if c = "1+4x not a square" then 
            vals := List( ValsFunction2(P), X -> [X]);
        elif c = "unique x so that 1-wx^2 is not a square" then 
            vals := [[ValsFunction10(P,W)]];
        elif c = "unique x so that x^2-w is not a square" then 
            vals := [[ValsFunction10a(P,W)]];
        elif c = "x ne -1" then 
            vals := List([0..P-2], X -> [X]);
        elif c = "x ne -1,-1/2" then 
            vals := List(Difference([0..P-2], [(P-1)/2]), X -> [X]);
        elif c = "x ne -2" then 
            vals := List(Difference([0..P-1], [P-2]), X -> [X]);
        elif c = "x ne 0" then 
            vals := List([1..P-1], X -> [X]);
        elif c = "x ne 0, equivalence classes {x,-x,1/x,-1/x}" then 
            vals := List(ValsFunction13(P), X -> [X]);
        elif c = "x ne 0, x~-x" then 
            vals := List( [1..(P-1)/2], X -> [X] );
        elif c = "x ne 0, x~x^-1" then 
            vals := List( ValsFunction1(P, [1..P-1]), X -> [X]);
        elif c = "x ne 0,-2, x~-x-2" then 
            vals := List(Union([1..(P-3)/2], [P-1]), X -> [X]);
        elif c = "x ne 0, x~ax if a^3=1" then 
            vals := List(ValsFunction5(P, 3, 1), X -> [X]);
        elif c = "x ne 0, x~ax if a^4=1" then 
            vals := List(ValsFunction5(P, 4, 1), X -> [X]);
        elif c = "x ne 0, x~ax if a^5=1" then 
            vals := List(ValsFunction5(P, 5, 1), X -> [X]);
        elif c = "x ne 0, x~ax if a^6=1" then 
            vals := List(ValsFunction5(P, 6, 1), X -> [X]);
        elif c = "x ne 0, x~ax if a^7=1" then 
            vals := List(ValsFunction5(P, 7, 1), X -> [X]);
        elif c = "x ne 0,-1" then 
            vals := List([1..P-2], X -> [X]);
        elif c = "x ne 0,-1, x~-1-x" then 
            vals := List( [1..(P-1)/2], X -> [X] );
        elif c = "x ne 0,-1,2,1/2" then 
            vals := List(Difference([1..P-2], [2, (P+1)/2]), X -> [X]);
        elif c = "x ne 0,-1/4" then 
            vals := List(Difference([1..P-1], [((P-1)/4) mod P]), X -> [X]);
        elif c = "x ne 0,-w, x~-w-x" then 
            vals := List(Filtered([1..P-1], Y -> ((-W-Y) mod P) >= Y), X->[X]);
        elif c = "x ne 0,-w,2w,w/2" then 
            vals := List(Difference([1..P-1], ([-W,2*W,W/2] mod P)), X -> [X]);
        elif c = "x ne 0,1" then 
            vals := List([2..P-1], X -> [X]);
        elif c = "x ne 1" then 
            vals := List(Difference( [0..P-1], [1] ), X -> [X]);
        elif c = "x~-1-x" then 
            vals := List([0..(P-1)/2], X -> [X]);
        elif c = "x~1-x" then 
            vals := List([1..(P+1)/2], X -> [X]);
        elif c = "x~-x" then 
            vals := List( [0..(P-1)/2], X -> [X] );
        elif c = "x~-x-2" then 
            vals := List(Union([0..(P-3)/2], [P-1]), X -> [X]);
        elif c = "x~w-x" then 
            vals := List(Filtered([0..P-1], Y->((W-Y) mod P)>=Y), X -> [X]);
        elif c = "x~ax if a^3=1" then 
            vals := List(ValsFunction5(P, 3, 0), X -> [X]);
        elif c = "x~ax if a^4=1" then 
            vals := List(ValsFunction5(P, 4, 0), X -> [X]);
        elif c = "x=0,1" then 
            vals := List([0,1], X -> [X]);
        elif c = "x=0,1,w" then 
            vals := List([0,1,W], X -> [X]);
        elif c = "x=0,1,w,w^2,w^3" then 
            vals := List([0,1,W,W^2,W^3], X -> [X mod P]);
        elif c = "x=1,w,...,w^(2gcd(p-1,3)-1)" then 
            vals := List([0..2*Gcd(P-1,3)-1], X -> [(W^X) mod P]);
        elif c = "x=1,w,...,w^(gcd(p-1,8)-1)" then 
            vals := List([0..Gcd(P-1,8)-1], X -> [(W^X) mod P]);
        elif c = "x=w,w^2" then 
            vals := List([W,W^2], X -> [X mod P]);
        elif c = "x=w,w^2,...,w^((p-3)/2)" then 
            vals := List([1..(P-3)/2], X -> [(W^X) mod P]);
        elif c = "x=w,w^2,...,w^6" then 
            vals := List([1..6], X -> [(W^X) mod P]);
        elif c = "x=w,w^2,w^3,w^4" then 
            vals := List([1..4], X -> [(W^X) mod P]);
        elif c = "x=w^2,w^3,w^4,w^5" then 
            vals := List([2..5], X -> [(W^X) mod P]);
        elif c = "x=w^3,w^4,...,w^8" then 
            vals := List([3..8], X -> [(W^X) mod P]);
        elif c = "x=w^4,w^5,w^6,w^7" then 
            vals := List([4..7], X -> [(W^X) mod P]);
        fi;
        if P>3 and vals = [] then Error("conditions not found - para 1"); fi;
        return TranslatedLiePRings( L, P, para, vals, flag );
    fi;
  
    # 2 parameter 
    if Length(para) = 2 then 

        if c = "[x,y]~[-x,y]" then
            vals := Cartesian([0..(P-1)/2], [0..(P-1)]);
        elif c = "[x,y]~[x,-y]" then
            vals := Cartesian([0..(P-1)], [0..(P-1)/2]);
        elif c = "[x,y]~[y,x]" then
            vals := Concatenation(List([0..P-1], X -> List([0..X],Y->[X,Y])));
        elif c = "[x,y]~[x',y'] if y^2-wx^2=y'^2-wx'^2" then
            vals := ValsFunction6(P);
        elif c = "[x,y]~[ax,y] if a^3=1" then
            vals := Cartesian( ValsFunction5(P,3,0), [0..(P-1)]);
        elif c = "[x,y]~[x,ay] if a^4=1" then
            vals := Cartesian( [0..(P-1)], ValsFunction5(P,4,0));
        elif c = "[x,y]~[±x,ay] if a^3=1" then
            vals := Cartesian( [0..(P-1)/2], ValsFunction5(P,3,0));
        elif c = "x ne -1, (1+x)y=1" then
            vals := List([0..(P-2)], X -> [X, (1+X)^-1 mod P]);
        elif c = "x ne -2" then
            vals := Cartesian( Difference([0..(P-1)], [(P-2)]), [0..(P-1)]);
        elif c = "x ne -2w" then
            vals := Cartesian( Difference([0..(P-1)], [(P-2*W) mod P]), [0..(P-1)]);
        elif c = "x ne 0" then
            vals := Cartesian( [1..(P-1)], [0..(P-1)]);
        elif c = "x ne 0, 4x+y^2 not a square" then
            vals := ValsFunction15(P);
        elif c = "x ne 0, [x,y]~[-x,-y-2]" then
            vals := Cartesian( [1..(P-1)/2], [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[-x,-y]" then
            vals := Cartesian( [1..(P-1)/2], [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[-x,y]~[x,-y]" then
            vals := Cartesian( [1..(P-1)/2], [0..(P-1)/2]);
        elif c = "x ne 0, [x,y]~[a^4x,ay] if a^5=1" then
            vals := Cartesian (ValsFunction5(P,5,1), [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[ax,a^2y] if a^3=1" then
            vals := Cartesian (ValsFunction5(P,3,1), [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[ax,a^2y] if a^6=1" then
            vals := Cartesian (ValsFunction5(P,6,1), [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[ax,a^3y] if a^5=1" then
            vals := Cartesian (ValsFunction5(P,5,1), [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[ax,a^3y] if a^6=1" then
            vals := Cartesian (ValsFunction5(P,6,1), [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[ax,ay] if a^3=1" then
            vals := Cartesian (ValsFunction5(P,3,1), [0..(P-1)]);
        elif c = "x ne 0, [x,y]~[x,-y]" then
            vals := Cartesian( [1..(P-1)], [0..(P-1)/2]);
        elif c = "x ne 0, [x,y]~[x,-y]~[-x,iy] if i^2=-1" then
            vals := Cartesian( [1..(P-1)/2], [0..(P-1)/2]);
        elif c = "x ne 0, unique y so that 1-wy^2 is not a square, [x,y]~[-x,y]" then
            vals := Cartesian([1..(P-1)/2], [ValsFunction10(P,W)]);
        elif c = "x ne 0, unique y so that wy^2=2, [x,y]~[-x,y]" then
            Y := IsSquareModP(P,((2/W) mod P)); if not IsInt(Y) then Error("no int"); fi;
            vals := Cartesian([1..(P-1)/2], [Y]);;
        elif c = "x ne 0, unique y so that y^2=2, [x,y]~[-x,y]" then
            Y := IsSquareModP(P,2); if not IsInt(Y) then Error("no int"); fi;
            vals := Cartesian([1..(P-1)/2], [Y]);;
        elif c = "x ne 0, unique y so that y^2=2w^2, [x,y]~[-x,y]" then
            Y := IsSquareModP(P,((2*W^2) mod P)); if not IsInt(Y) then Error("no int"); fi;
            vals := Cartesian([1..(P-1)/2], [Y]);;
        elif c = "x ne 0, unique y so that y^2=2w^3, [x,y]~[-x,y]" then
            Y := IsSquareModP(P,((2*W^3) mod P)); if not IsInt(Y) then Error("no int"); fi;
            if IsInt(Y) then vals := Cartesian([1..(P-1)/2], [Y]); fi;
        elif c = "x ne 0, y=1,-1" then
            vals := Cartesian( [1..(P-1)], [1,(P-1)]);
        elif c = "x ne 0, y=1,-1, [x,y]~[-x,y]" then
            vals := Cartesian( [1..(P-1)/2], [1,(P-1)]);
        elif c = "x ne 0, y=1,w, [x,y]~[-x,y]" then
            vals := Cartesian( [1..(P-1)/2], [1,W]);
        elif c = "x ne 0, y=w,w^2,...,w^((p-3)/2)" then
            vals := Cartesian( [1..(P-1)], List([1..(P-3)/2], X -> ((W^X) mod P)));
        elif c = "x ne 0, y=w,w^2,...,w^6, [x,y]~[ax,y] if a^7=1" then
            vals := Cartesian( ValsFunction5(P,7,1), List([1..6], X -> ((W^X) mod P)));
        elif c = "x ne 0, y=w,w^2,w^3,w^4, [x,y]~[ax,y] if a^5=1" then
            vals := Cartesian( ValsFunction5(P,5,1), List([1..4], X -> ((W^X) mod P)));
        elif c = "x ne 0, y=w^2,w^3,w^4,w^5, [x,y]~[ax,y] if a^3=1" then
            vals := Cartesian( ValsFunction5(P,3,1), List([2..5], X -> ((W^X) mod P)));
        elif c = "x ne 1" then
            vals := Cartesian( Difference([0..(P-1)], [1]), [0..(P-1)]);
        elif c = "x ne 1-w, [x,y]~[x,-y]~[-x+2(1-w),iy] if i^2=-1" then
            vals := Cartesian(ValsFunction14(P,W,1), [0..(P-1)/2]);
        elif c = "x ne 1-w^2, [x,y]~[x,-y]~[-x+2(1-w^2),iy] if i^2=-1" then
            vals := Cartesian(ValsFunction14(P,W,2), [0..(P-1)/2]);
        elif c = "x ne 1-w^3, [x,y]~[x,-y]~[-x+2(1-w^3),iy] if i^2=-1" then
            vals := Cartesian(ValsFunction14(P,W,3), [0..(P-1)/2]);
        elif c = "x,y ne 0, [x,y]~[-x,-y]" then
            vals := Cartesian([1..(P-1)], [1..(P-1)/2]);
        elif c = "x,y ne 0, [x,y]~[ax,a^3y] if a^4=1" then
            vals := Cartesian(ValsFunction5(P,4,1), [1..(P-1)]);
        elif c = "x,y ne 0, [x,y]~[xy^-2,y^-1]" then
            vals := Cartesian([1..(P-1)], ValsFunction1(P,[1..P-1]));
        elif c = "x,y ne 0, [x,y]~[y,x]" then
            vals := Concatenation(List([1..(P-1)], X -> List([1..X], Y -> [X,Y])));
        elif c = "x,y ne 0, x ne y, [x,y]~[y,x]" then
            vals := Concatenation(List([1..(P-1)], X -> List([1..X-1], Y -> [X,Y])));
        elif c = "x=0,-1, y=0,1" then
            vals := Cartesian( [0,P-1], [0,1]);
        elif c = "x=0,1" then
            vals := Cartesian( [0,1], [0..(P-1)]);
        elif c = "x=0,1, y=0,1" then
            vals := Cartesian( [0,1], [0,1]);
        elif c = "x=0,1,w" then
            vals := Cartesian( [0,1,W], [0..(P-1)]);
        elif c = "x=0,1,w, y=0,1" then
            vals := Cartesian( [0,1,W], [0,1]);
        elif c = "x=0,1,w, y=1,-1" then
            vals := Cartesian( [0,1,W], [1,P-1]);
        elif c = "x=0,1,w, y=w,w^2,...,w^((p-3)/2)" then
            vals := Cartesian( [0,1,W], List([1..(P-3)/2], X -> ((W^X) mod P)));
        elif c = "x=1,w" then
            vals := Cartesian( [1,W], [0..(P-1)]);
        elif c = "x=1,w,...,w^(2gcd(p-1,3)-1), y ne 0, [x,y]~[x,ay] if a^6=1" then
            vals := List([0..(2*Gcd(P-1,3)-1)], X -> ((W^X) mod P));
            vals := Cartesian( vals, ValsFunction5(P,6,1));
        elif c = "x=1,w,...,w^(gcd(p-1,8)-1), y ne 0, [x,y]~[x,ay] if a^8=1" then
            vals := List([0..(Gcd(P-1,8)-1)], X -> ((W^X) mod P));
            vals := Cartesian( vals, ValsFunction5(P,8,1));
        elif c = "xy ne 1, [x,y]~[y,x]" then
            vals := Concatenation(List([0..P-1], X -> List([0..X],Y->[X,Y])));
            vals := Filtered(vals, X -> (X[1]*X[2] mod P)<>1); 
        elif c = "y ne 0, [x,y]~[-x,y]" then
            vals := Cartesian([0..(P-1)/2], [1..(P-1)]);
        elif c = "y ne 0, [x,y]~[a^2x,ay] if a^4=1" then
            vals := Cartesian([0..(P-1)], ValsFunction5(P,4,1));
        elif c = "y ne 0, [x,y]~[x,-y]" then
            vals := Cartesian([0..(P-1)], [1..(P-1)/2]);
        elif c = "y ne 0, [x,y]~[x,ay] if a^3=1" then
            vals := Cartesian([0..(P-1)], ValsFunction5(P,3,1));
        elif c = "y ne 0, [x,y]~[x,ay] if a^4=1" then
            vals := Cartesian([0..(P-1)], ValsFunction5(P,4,1));
        elif c = "y ne 0, [x,y]~[x,ay] if a^5=1" then
            vals := Cartesian([0..(P-1)], ValsFunction5(P,5,1));
        elif c = "y ne 1/2, [x,y]~[-x,1-y]" then
            vals := Cartesian([0..(P-1)], [1..(P-1)/2]);
        elif c = "y=1,w, [x,y]~[-x,y]" then
            vals := Cartesian([0..(P-1)/2], [1,W]);
        elif c = "y=1,w, [x,y]~[ax,y] if a^8=1" then
            vals := Cartesian(ValsFunction5(P,8,0), [1,W]);
        elif c = "y=1,w,...,w^5, [x,y]~[ax,y] if a^6=1" then
            vals := List([0..5], X -> ((W^X) mod P));
            vals := Cartesian(ValsFunction5(P,6,0), vals);
        elif c = "y=1,w,w^2, [x,y]~[ax,y] if a^3=1" then
            vals := List([0..2], X -> ((W^X) mod P));
            vals := Cartesian(ValsFunction5(P,3,0), vals);
        elif c = "y=w^2,w^3, [x,y]~[ax,y] if a^8=1" then
            vals := List([2,3], X -> ((W^X) mod P));
            vals := Cartesian(ValsFunction5(P,8,0), vals);
        elif c = "y=w^2,w^3,w^4,w^5" then
            vals := Cartesian([0..(P-1)], List([2..5], X -> ((W^X) mod P)));
        elif c = "y=w^4,w^5,w^6,w^7, [x,y]~[ax,y] if a^8=1" then
            vals := List([4..7], X -> ((W^X) mod P));
            vals := Cartesian(ValsFunction5(P,8,0), vals);
        fi;
        if P>3 and vals = [] then Error("conditions not found - para 2"); fi;
        return TranslatedLiePRings( L, P, para, vals, flag );
    fi;

    # 3 parameter 
    if Length(para) = 3 then 

        if c = "[x,y,z]~[-x,-y,-z]" then
            vals := ValsFunction17(P);
        elif c = "[x,y,z]~[-x,y,z]" then
            vals := Cartesian([0..(P-1)/2],[0..P-1],[0..P-1]);
        elif c = "[x,y,z]~[x,y,-z]" then
            vals := Cartesian([0..P-1],[0..P-1],[0..(P-1)/2]);
        elif c = "[x,y,z]~[z,y,x]" then
            vals := Cartesian([0..(P-1)],[0..(P-1)],[0..(P-1)]);
            vals := Filtered(vals, X -> X[1]<=X[3]); 
        elif c = "unique z so that z^2-4 is not a square, [x,y,z]~[y,x,z]" then
            Z := ValsFunction11(P);
            vals := Concatenation(List([0..P-1], X -> List([0..X],Y->[X,Y,Z])));
        elif c = "x ne -2, [x,y,z]~[x,y,-z]" then
            vals := Cartesian(Difference([0..P-1],[P-2]),[0..P-1],[0..(P-1)/2]);
        elif c = "x ne -2w, [x,y,z]~[x,y,-z]" then
            vals := Cartesian(Difference([0..P-1],[(-2*W) mod P]),[0..P-1],[0..(P-1)/2]);
        elif c = "x ne 0" then
            vals := Cartesian([1..(P-1)], [0..(P-1)], [0..(P-1)]);
        elif c = "x ne 0, [x,y,z]~[-x,-y,z]" then
            vals := Cartesian([1..(P-1)/2], [0..(P-1)], [0..(P-1)]);
        elif c = "x ne 0, [x,y,z]~[a^3x,a^4y,az] if a^5=1" then
            vals := Cartesian(ValsFunction5(P,5,1), [0..(P-1)], [0..(P-1)]);
        elif c = "x ne 0, [x,y,z]~[ax,a^2y,az] if a^4=1" then
            vals := Cartesian(ValsFunction5(P,4,1), [0..(P-1)], [0..(P-1)]);
        elif c = "x ne 0, [x,y,z]~[ax,y,a^2z] if a^3=1" then
            vals := Cartesian(ValsFunction5(P,3,1), [0..(P-1)], [0..(P-1)]);
        elif c = "x ne 0, y=w^2,w^3,w^4,w^5, [x,y,z]~[ax,y,a^2z] if a^6=1" then
            vals := Cartesian(ValsFunction5(P,6,1),List([2,3,4,5], X -> ((W^X) mod P)), [0..P-1]);
        elif c = "x ne 0, z=w,w^2,w^3,w^4, [x,y,z]~[ax,a^3y,z] if a^5=1" then
            vals := Cartesian(ValsFunction5(P,5,1), [0..P-1], List([1,2,3,4], X -> ((W^X) mod P)));
        elif c = "x,z ne 0, [x,y,z]~[-x,y,-z]~[x,-y,z]" then
            vals := Cartesian([1..(P-1)],[0..(P-1)/2],[1..(P-1)/2]);
        elif c = "x=0,1" then
            vals := Cartesian([0,1],[0..(P-1)],[0..(P-1)]);
        elif c = "x=0,1, y=0,1" then
            vals := Cartesian([0,1],[0,1],[0..(P-1)]);
        elif c = "x=0,1, y=0,1, z ne 0,-1" then
            vals := Cartesian([0,1],[0,1],[1..(P-2)]);
        elif c = "x=0,1, y=0,1, z=0,1" then
            vals := Cartesian([0,1],[0,1],[0,1]);
        elif c = "x=0,1, y=0,1, z=0,1,w" then
            vals := Cartesian([0,1],[0,1],[0,1,W]);
        elif c = "x=1,w,...,w^(2gcd(p-1,3)-1), y ne 0, [x,y,z]~[x,ay,±a^2z] if a^3=1" then
            vals := Cartesian(List([0..2*Gcd(P-1,3)-1], X -> ((W^X) mod P)),
                              ValsFunction5(P,3,1), [0..(P-1)/2]);
        elif c = "x=w,w^2,...,w^((p-3)/2), y=1,w" then
            vals := Cartesian(List([1..(P-3)/2], X -> ((W^X) mod P)), [1,W], [0..(P-1)]);
        elif c = "y ne 0, [x,y,z]~[x,-y,-z]" then
            vals := Cartesian([0..(P-1)],[1..(P-1)/2],[0..(P-1)]);
        elif c = "y ne 0, [x,y,z]~[zy,y,x/y]" then
            vals := ValsFunction18(P); 
        elif c = "y ne 0,1, (x+y)(1+z)=1, [x,y,z]~[zy,y,x/y]" then
            vals := ValsFunction18a(P); 
        elif c = "z=1,w,w^2,w^3,w^4, [x,y,z]~[x,ay,z] if a^5=1" then
            vals := Cartesian([0..(P-1)],ValsFunction5(P,5,0),([1,W,W^2,W^3,W^4] mod P));
        fi;
        if P>3 and vals = [] then Error("conditions not found - para 3"); fi;
        return TranslatedLiePRings( L, P, para, vals, flag );
    fi;

    # 4 parameters
    if Length(para) = 4 then 

        if c = "[x,y,z,t]~[t+1,z+1,y-1,x-1]" then
            vals := ValsFunction16(P);
        elif c = "[x,y,z,t]~[x,-y,z,-t]" then
            vals := Cartesian([0..(P-1)],[1..(P-1)/2],[0..P-1],[0..P-1]);
            vals := Concatenation( vals, Cartesian([0..(P-1)],[0],[0..P-1],[0..(P-1)/2]));
        elif c = "x ne 0, [x,y,z,t]~[-x,y,t,z]" then
            vals := Cartesian([1..(P-1)/2],[0..P-1],[0..P-1],[0..P-1]);
        elif c = "x ne 0, [x,y,z,t]~[-x,y,z,-t]" then
            vals := Cartesian([1..(P-1)/2],[0..P-1],[0..P-1],[0..P-1]);
        elif c = "x ne 0, [x,y,z,t]~[ax,a^3y,a^4z,at] if a^5=1" then
            vals := Cartesian(ValsFunction5(P,5,1),[0..P-1],[0..P-1],[0..P-1]);
        elif c = "x ne 0, [x,y,z,t]~[ax,a^3y,a^4z,t] if a^5=1" then
            vals := Cartesian(ValsFunction5(P,5,1),[0..P-1],[0..P-1],[0..P-1]);
        fi;
        if P>3 and vals = [] then Error("conditions not found - para 4"); fi;
        return TranslatedLiePRings( L, P, para, vals, flag );
    fi;
  
end );

