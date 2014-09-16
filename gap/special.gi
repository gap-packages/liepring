
##
## Specialize Lie ring with given values
##
## P : a prime or false
## para : a list of parameters (not p or w)
## vals : the values of the parameters
##
SpecialiseLiePRingNC := function( L, P, para, vals )
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
        T.tab[i] := WordByExps(LRReduceExp(T,ExpsByWord(T.dim,T.tab[i])));
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
end;

##
## check if P satisfies conditions of L and if cl(L) < P
##
IsValidPrime := function(L, P)
    local c;
    if not IsPrimeInt(P) then return false; fi;
    if not HasLibraryConditions(L) then return true; fi;
    if P=2 then return false; fi;
    if P=3 and PClassOfLiePRing(L) > 2 then return false; fi;
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
    fi;
    return true;
end;

##
## Specialize with checking 
##
SpecialiseLiePRing := function( L, P, para, vals )
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
end;

##
## Specialize prime only
##
SpecialisePrimeOfLiePRing := function( L, P )
    return SpecialiseLiePRing( L, P, [], [] );
end;
SpecialisePrimeOfLiePRingNC := function( L, P )
    return SpecialiseLiePRingNC( L, P, [], [] );
end;

## 
## Lie rings in a family
##

LiePRingsInFamily := function( arg )
    local L, S, d, P, W, c, para, vals, fals, X, Y, Z,
          p, w, x, y, z, t, j, k, m, n, r, s, u, v;

    # get Lie ring
    L := arg[1];
    S := SCTable(Zero(L));
    d := DimensionOfLiePRing(L);

    # get prime
    if Length(arg) = 2 then 
        P := arg[2];
        if IsInt(S.prime) and P <> S.prime then 
            Print("cannot specialise prime to ",P," \n");
            return fail;
        elif not IsValidPrime(L,P) then 
            return fail;
        fi;
    else
        P := S.prime;
    fi;
    if not IsInt(P) then 
        Print("need to specify a prime \n");
        return fail;
    fi;

    if not IsBound(S.param) and IsInt(S.prime) then 
        return [L]; 
    elif not IsBound(S.param) then
        return [SpecialiseLiePRingNC( L, P, [], [] )];
    fi;

    # assign variable names
    p := IndeterminateByName("p");
    w := IndeterminateByName("w");
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

    para := Difference( S.param, [w] );
    vals := [];
    fals := 0;
    c := L!.LibraryConditions[1];
    W := PrimitiveRootMod(P);

    # 0 parameters
    if Length(para) = 0 then 
        return [SpecialiseLiePRingNC( L, P, [], [] )];
    fi;

    # 1 parameter
    if Length(para) = 1 then 

        # para = [x]
        if (Length(c) = 0) or (c = "all x") then 
            vals := List([0..P-1], X -> [X]);
        elif c = "x ne -1" then 
            vals := List([0..P-2], X -> [X]);
        elif c = "x ne -1,-1/2" then 
            vals := List(Difference([0..P-2], [(P-1)/2]), X -> [X]);
        elif c = "x ne -2" then 
            vals := List(Difference([0..P-1], [P-2]), X -> [X]);
        elif c = "x ne 0" then 
            vals := List([1..P-1], X -> [X]);
        elif c = "x ne 0,-1" then 
            vals := List([1..P-2], X -> [X]);
        elif c = "x ne 0,-1,2,1/2" then 
            vals := List(Difference([1..P-2], [2, (P+1)/2]), X -> [X]);
        elif c = "x ne 0,-1/4" then 
            vals := List(Difference([1..P-1], [((P-1)/4) mod P]), X -> [X]);
        elif c = "x ne 0,-w,2w,w/2" then 
            vals := List(Difference([1..P-1], ([-W,2*W,W/2] mod P)), X -> [X]);
        elif c = "x ne 0,1" then 
            vals := List([2..P-1], X -> [X]);
        elif c = "x ne 1" then 
            vals := List(Difference( [0..P-1], [1] ), X -> [X]);
        elif c = "x=0,1,w" then 
            vals := List([0,1,W], X -> [X]);
        elif c = "x=1/2" then 
            vals := List([(P+1)/2], X -> [X]);
        elif c = "x=1/4" then 
            vals := List([(P+1)/4], X -> [X]);
        elif c = "x=0,1,w,w^2,w^3" then 
            vals := List(([0,1,W,W^2,W^3] mod P), X -> [X]);
        elif c = "x=w,w^2,w^3,w^4" or c = "x=w^i, i=1,2,3,4" then 
            vals := List(([W, W^2, W^3, W^4] mod P), X -> [X]);
        elif c = "x=w,w^2,w^3,w^4,w^5,w^6" or c = "x=w^i, i=1,2,3,4,5,6" then 
            vals := List(([W, W^2, W^3, W^4, W^5, W^6] mod P), X -> [X]);
        elif c = "x=w^2,w^3,w^4,w^5" or c = "x=w^i, i=2,3,4,5" then 
            vals := List(([W^2, W^3, W^4, W^5] mod P), X -> [X]);
        elif c = "x=w^4,w^5,w^6,w^7" then 
            vals := List(([W^4, W^5, W^6, W^7] mod P), X -> [X]);
        elif c = "x=w/2" then 
            vals := [[W/2 mod P]];
        elif c = "x=w^i, 0 le i lt 2gcd[p-1,3]" then 
            vals := List([0..2*Gcd(P-1,3)-1], X -> [W^X mod P]);
        elif c = "x=w^i, 0 le i lt gcd[p-1,8]" then 
            vals := List([0..Gcd(P-1,8)-1], X -> [W^X mod P]);
        elif c = "all x, x~-x" or c = "all x,x~-x" or c = "x~-x" then 
            vals := List( [0..(P-1)/2], X -> [X] );
        elif c = "x ne 0, x~-x" then 
            vals := List( [1..(P-1)/2], X -> [X] );
        elif c = "x ne 0, x~x^-1" or c = "x ne 0,x~x^-1" then 
            vals := List( ValsFunction1(P, [1..P-1]), X -> [X]);
        elif c = "x ne 0,1,-1, x~x^-1" then 
            vals := List( ValsFunction1(P, [2..P-2]), X -> [X]);
        elif c = "all x, x~1-x" then 
            vals := List([1..(P+1)/2], X -> [X]);
        elif c = "x ne 0,-2, x~-x-2" then 
            vals := List(Union([1..(P-3)/2], [P-1]), X -> [X]);
        elif c = "all x, x~w-x" then 
            vals := List(Filtered([0..P-1], Y -> ((W-Y) mod P) >= Y), X -> [X]);
        elif c = "x ne 0,-w, x~-w-x" then 
            vals := List(Filtered([1..P-1], Y -> ((-W-Y) mod P) >= Y), X->[X]);
        elif c = "1+4x not a square mod p" or
             c = "all x, 1+4x not a square mod p" then 
            vals := List( ValsFunction2(P), X -> [X]);
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p" then 
            vals := List(ValsFunction5(P, 5, 1), X -> [X]);
        elif c = "x ne 0,  x~x' if x^4=x'^4 mod p" or 
             c = "x ne 0, x~x' if x^4=x'^4 mod p" then 
            vals := List(ValsFunction5(P, 4, 1), X -> [X]);
        elif c = "x ne 0, x~x' if x^3=x'^3  mod  p" or
             c = "x ne 0, x~x' if x^3=x'^3 mod p" or
             c = "x ne 0,x~x' if x^3=x'^3 mod p" then 
            vals := List(ValsFunction5(P, 3, 1), X -> [X]);
        elif c = "x~x' if x^4=x'^4 mod p" or
             c = "all x, x~x' if x^4=x'^4 mod p" then 
            vals := List(ValsFunction5(P, 4, 0), X -> [X]);
        elif c = "all x, x~x' if x^3=x'^3 mod p" then 
            vals := List(ValsFunction5(P, 3, 0), X -> [X]);
        elif c = "x ne 0, x~x' if x^6=x'^6 mod p" then 
            vals := List(ValsFunction5(P, 6, 1), X -> [X]);
        elif c = "x ne 0, x~x' if x^7=x'^7 mod p" then 
            vals := List(ValsFunction5(P, 7, 1), X -> [X]);
        elif c = "x FIXED with x^2-w not a square mod p" then 
            vals := [[ValsFunction10a(P,W)]];
        elif c = "x ne 0, equivalence classes {x,-x,1/x,-1/x}" then 
            vals := List(ValsFunction13(P), X -> [X]);
    
        # para = [y]
        elif c = "all y" then 
            vals := List([0..P-1], X -> [X]);
        elif c = "y=0,1" then 
            vals := List([0,1], X -> [X]);
        elif c = "y=0,1,w" then 
            vals := List([0,1,W], X -> [X]);
        elif c = "y=w,w^2" then 
            vals := List(([W,W^2] mod P), X -> [X]);
        elif c = "y=w^3,w^4,w^5,w^6,w^7,w^8" then 
            vals := List([3..8], X -> [W^X mod P]);
        elif c = "y=w^i, i=1,2,3,4" then 
            vals := List([1..4], X -> [W^X mod P]);
        elif c = "y=w^i, i=2,3,4,5" then
            vals := List([2..5], X -> [W^X mod P]);
        elif c = "y=w,w^2" then
            vals := List(([W,W^2] mod P), X -> [X]);
        elif c = "y=w^3,w^4,w^5,w^6,w^7,w^8" then 
            vals := List(([W^3,W^4,W^5,W^6,W^7,W^8] mod P), X -> [X]);
        elif c = "y FIXED such that 1-wy^2 is not a square mod p" then
            vals := [[ValsFunction10(P,W)]];
        elif c = "all y, y~-y" then 
            vals := List([0..(P-1)/2], X -> [X]);
        elif c = "all y, y~-y-2" then 
            vals := List(Union([0..(P-3)/2], [P-1]), X -> [X]);
        elif c = "all y, y~y' if y^4=y'^4 mod p" then 
            vals := List(ValsFunction5(P, 4, 0), X -> [X]);
        elif c = "y ne 0, y~y' if y^3=y'^3 mod p" then 
            vals := List(ValsFunction5(P, 3, 1), X -> [X]);
        elif c = "y ne 0, y~y' if y^4=y'^4 mod p" then 
            vals := List(ValsFunction5(P, 4, 1), X -> [X]);

        # para = [z]
        elif c = "all z" then 
            vals := List([0..P-1], X -> [X]);
        elif c = "z ne 0" then 
            vals := List([1..P-1], X -> [X]);
        elif c = "z ne 0, z~-z" then 
            vals := List([1..(P-1)/2], X -> [X]);
        elif c = "z=0,1" then 
            vals := List([0,1], X -> [X]);
        elif c = "z=1/2" then 
            vals := List([(P+1)/2], X -> [X]);
        elif c = "z=w^i, i=1,2,3,4" then 
            vals := List([1,2,3,4], X -> [W^X mod P]);

        # para = [t]
        elif c = "all t" then 
            vals := List([0..P-1], X -> [X]);
        elif c = "t=0,1,w" then 
            vals := List([0,1,W], X -> [X]);
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;
  
    # 2 parameter
    elif Length(para) = 2 then 

        # sort parameters correctly
        if x in para then # [x,y], [x,z], [x,t]
            para := [x, Difference(para,[x])[1]]; 
        elif y in para then # [y,z], [y,t]
            para := [y, Difference(para,[y])[1]]; 
        elif z in para then # [z, t]
            para := [z, Difference(para,[z])[1]]; 
        else
            Error("parameters out of range");
        fi;

        # para = [x,y];
        if Length(c) = 0 or c = "all x,y" then 
            vals := Cartesian( [0..P-1], [0..P-1]);
        elif c = "all x, y ne 0" then 
            vals := Cartesian([0..P-1], [1..P-1]);
        elif c = "x ne 0, all y" then 
            vals := Cartesian([1..P-1], [0..P-1]);
        elif c = "x ne 0, all y, y~-y" then 
            vals := Cartesian([1..P-1], [0..(P-1)/2]);
        elif c = "all x,y, x~-x" or c = "all x, x~-x, all y" then 
            vals := Cartesian([0..(P-1)/2], [0..P-1]);
        elif c = "all x, x~-x, y ne 0" then 
            vals := Cartesian([0..(P-1)/2], [1..P-1]);
        elif c = "x ne 0, all y, x~-x" or c = "x ne 0, x~-x, all y" then 
            vals := Cartesian([1..(P-1)/2], [0..P-1]);
        elif c = "x ne 0, all y, x~-x, y~-y" or 
             c = "x ne 0, x~-x, all y, y~-y" then 
            vals := Cartesian([1..(P-1)/2], [0..(P-1)/2]);
        elif c = "all x,y, y~-y" or c = "all x, all y, y~-y" then 
            vals := Cartesian([0..P-1], [0..(P-1)/2]);
        elif c = "all x, y ne 0, y~-y" then 
            vals := Cartesian([0..P-1], [1..(P-1)/2]);
        elif c = "y=1,w, all x, x~-x" then 
            vals := Cartesian([0..(P-1)/2], [1,W]);
        elif c = "x ne 0, x~-x, y=1,w" then 
            vals := Cartesian([1..(P-1)/2], [1,W]);
        elif c = "all x, x~-x, y=1/2" then 
            vals := Cartesian([0..(P-1)/2], [(P+1)/2]);
        elif c = "x ne 0, x~-x, y=w^-1" then 
            vals := Cartesian([1..(P-1)/2], [W^-1 mod P]);
        elif c = "x ne 0, y ne 0,1,-1, y~y^-1" then 
            vals := Cartesian([1..P-1], ValsFunction1(P,[2..P-2]));
        elif c = "x ne 0, y ne 0 ,y~y^-1" then 
            vals := Cartesian([1..P-1], ValsFunction1(P,[1..P-1]));
        elif c = "all x, y=w^2,w^3,w^4,w^5" then 
            vals := Cartesian([0..P-1], ([W^2, W^3, W^4, W^5] mod P));
        elif c = "all x, y=0,1" then 
            vals := Cartesian([0..P-1], [0,1]);
        elif c = "all x, y ne 1/2, y~1-y" then 
            vals := Cartesian([0..P-1], [1..(P-1)/2]);
        elif c = "x ne -1, y=[1+x]^-1 mod p" then 
            vals := List([0..P-2], X -> [X, (1+X)^-1 mod P]);
        elif c = "x ne -2, all y" then 
            vals := Cartesian( Difference([0..P-1],[P-2]), [0..P-1] );
        elif c = "x ne -2w, all y" then 
            vals := Cartesian( Difference([0..P-1], [(-2*W) mod P]),[0..P-1] );
        elif c = "y=1,w,w^2, all x, x1~x2 if x1^3=x2^3 mod p" or
             c = "y=1,w,w^2, all x, x~x' if x^3=x'^3 mod p" then 
            vals := Cartesian( ValsFunction5(P,3,0), [1,W,W^2]);
        elif c = "x ne 0, all y, x~x' if x^3=x'^3 mod p" then 
            vals := Cartesian( ValsFunction5(P,3,1), [0..P-1]);
        elif c = "all x,y, x~x' if x^3=x'^3 mod p" then 
            vals := Cartesian(ValsFunction5(P,3,0),[0..P-1]);
        elif c = "y=1,w, all x, x~x' if x^8=x'^8 mod p" then 
            vals := Cartesian(ValsFunction5(P,8,0),[1,W]);
        elif c = "y=1,w,w^2, all x, x1~x2 if x1^3=x2^3 mod p" then 
            vals := Cartesian(ValsFunction5(P,3,0),([1,W,W^2] mod P));
        elif c = "y=1,w,w^2,w^3,w^4,w^5, all x, x~x' if x^6=x'^6 mod p" then 
            vals := Cartesian(ValsFunction5(P,6,0), 
                              ([1,W,W^2,W^3,W^4,W^5] mod P));
        elif c = "y=w^2,w^3, all x, x~x' if x^8=x'^8 mod p" then 
            vals := Cartesian(ValsFunction5(P,8,0), ([W^2,W^3] mod P));
        elif c = "y=w^4,w^5,w^6,w^7, all x, x~x' if x^8=x'^8 mod p" then 
            vals := Cartesian(ValsFunction5(P,8,0), 
                              ([W^4,W^5,W^6,W^7] mod P));
        elif c = "x ne 0, all y, x~x' if x^3=x'^3 mod p" or 
             c = "x ne 0, x~x' if x^3=x'^3 mod p, all y" then 
            vals := Cartesian(ValsFunction5(P,3,1),[0..P-1]);
        elif c = "x ne 0, x~x' if x^6=x'^6 mod p, all y" then 
            vals := Cartesian(ValsFunction5(P,6,1),[0..P-1]);
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p, all y" then 
            vals := Cartesian(ValsFunction5(P,5,1),[0..P-1]);
        elif c = "y=w^i, i=1,2,3,4, x ne 0, x~x' if x^5=x'^5 mod p" then 
            vals := Cartesian(ValsFunction5(P,5,1),
                              ([W, W^2, W^3, W^4] mod P));
        elif c = "y=w^i, i=1,2,3,4,5,6, x ne 0, x~x' if x^7=x'^7 mod p" then 
            vals := Cartesian(ValsFunction5(P,7,1),
                              List([1..6],Y->(W^Y mod P)));
        elif c = "y=w^i, i=2,3,4,5, x ne 0, x~x' if x^3=x'^3 mod p" then
            vals := Cartesian(ValsFunction5(P,3,1),
                              List([2..5],Y->(W^Y mod P)));
        elif c = "x ne 1, all y" then 
            vals := Cartesian(Union([0],[2..P-1]), [0..P-1]);
        elif c = "x=0,-1, y=0,1" then 
            vals := Cartesian([0,P-1], [0,1]);
        elif c = "x=0,1, all y" then 
            vals := Cartesian([0,1], [0..P-1]);
        elif c = "x=0,1, y=0,1" then 
            vals := Cartesian([0,1], [0,1]);
        elif c = "x ne 1-w, x~-x+2[1-w], all y, y~-y" then 
            vals := Cartesian(ValsFunction14(P,W,1), [0..(P-1)/2]);
        elif c = "x ne 1-w^2, x~-x+2[1-w^2], all y, y~-y" then 
            vals := Cartesian(ValsFunction14(P,W,2), [0..(P-1)/2]);
        elif c = "x ne 1-w^3, x~-x+2[1-w^3], all y, y~-y" then 
            vals := Cartesian(ValsFunction14(P,W,3), [0..(P-1)/2]);
        elif c = "all x, y ne 0, y~y' if y^4=y'^4 mod p" then 
            vals := Cartesian([0..P-1], ValsFunction5(P,4,1));
        elif c = "all x, y ne 0, y~y' if y^5=y'^5 mod p" or
             c = "all x, y ne 0, y~y' if y^5=y'^5 mod p" then 
            vals := Cartesian([0..P-1], ValsFunction5(P,5,1));
        elif c = "all x,y, x~-x, y~y' if y^3 eq y'^3  mod  p" or
             c = "all x,y, x~-x, y~y' if y^3=y'^3 mod p" then 
            vals := Cartesian([0..(P-1)/2], ValsFunction5(P,3,0));
        elif c = "all x, y ne 0, y~y' if y^3=y'^3 mod p" or 
             c = "y ne 0, y~y' if y^3=y'^3 mod p, all x" then 
            vals := Cartesian([0..P-1], ValsFunction5(P,3,1));
        elif c = "all x, y ne 0, y~y' if y^4=y'^4 mod p" then 
            vals := Cartesian([0..P-1], ValsFunction5(P,4,1));
        elif c = "all x,y, y~y' if y^4=y'^4 mod p" then 
            vals := Cartesian([0..P-1], ValsFunction5(P,4,0));
        elif c = "x=w^i, 0 le i lt 2gcd[p-1,3], y ne 0, y~y' if y^6=y'^6 mod p"
        then 
            vals := Cartesian(List([0..2*Gcd(P-1,3)-1], Y->((W^Y) mod P)),
                              ValsFunction5(P,6,1));
        elif c = "x=w^i, 0 le i lt gcd[p-1,8], y ne 0, y~y' if y^8=y'^8 mod p"
        then 
            vals := Cartesian(List([0..Gcd(P-1,8)-1], Y->((W^Y) mod P)),
                              ValsFunction5(P,8,1));
        elif c = "all x,y, [x,y]~[y,x]" then 
            vals := Concatenation(List([0..P-1], X -> List([0..X],Y->[X,Y])));
        elif c = "x,y ne 0, x ne y, [x,y]~[y,x]" then 
            vals := Concatenation(List([1..P-1], X -> List([1..X-1],Y->[X,Y])));
        elif c = "all x,y such that xy ne 1, [x,y]~[y,x]" then 
            vals := Concatenation(List([0..P-1], X -> List([0..X],Y->[X,Y])));
            vals := Filtered(vals, X -> (X[1]*X[2] mod P)<>1); 
        elif c = "all x,y, [x,y]~[x',y'] if y^2-wx^2=y'^2-wx'^2 mod p" then 
            vals := ValsFunction6(P);
        elif c ="x ne 0, x~-x, y FIXED such that 1-wy^2 is not a square mod p"
        then 
            vals := Cartesian([1..(P-1)/2], [ValsFunction10(P,W)]);
        elif c = "x ne 0, all y, 4x+y^2 not a square mod p" then 
            vals := ValsFunction15(P);
        elif c = "If 2 is a square mod p and y^2=2, x ne 0, x~-x, y~-y" then 
            Y := IsSquareModP(P,2);
            if IsInt(Y) then vals := Cartesian([1..(P-1)/2], [Y]); fi;
        elif c = "If 2 is a square mod p and y^2=2w^2, x ne 0, x~-x, y~-y" then
            Y := IsSquareModP(P,2*W^2);
            if IsInt(Y) then vals := Cartesian([1..(P-1)/2], [Y]); fi;
        elif c = "If 2 is not a square mod p and wy^2=2, x ne 0, x~-x, y~-y" 
        then
            Y := IsSquareModP(P,2/W);
            if IsInt(Y) then vals := Cartesian([1..(P-1)/2], [Y]); fi;
        elif c = "If 2 is not a square mod p and y^2=2w^3, x ne 0, x~-x, y~-y"
        then 
            Y := IsSquareModP(P,2*W^3);
            if IsInt(Y) then vals := Cartesian([1..(P-1)/2], [Y]); fi;
        elif d = 6 and c = "See note6.62" then 
            vals := ValsFunction9(P);
        elif d = 7 and c = "See Notes6.178" then 
            fals := ValsFunction24(P); vals := fals[1];
        elif c = "See note1dec5.1" then 
            vals := ValsFunction9(P);

        # para = [x,z]
        elif c = "all x,z" then 
            vals := Cartesian([0..P-1], [0..P-1]);
        elif c = "x ne 0, all z" then 
            vals := Cartesian([1..P-1], [0..P-1]);
        elif c = "all x, z=1,w" then
            vals := Cartesian([0..P-1], [1,W]);
        elif c = "x,z ne 0, [x,z]~[z,x]" then 
            vals := Concatenation(List([1..P-1], X -> List([1..X],Y->[X,Y])));
        elif c = "x ne 0, x~-x, z=w^-1" then 
            vals := Cartesian([1..(P-1)/2], [W^-1]);
        elif c = "x ne 0, x~x' if x^4=x'^4 mod p, z=w^-1" then 
            vals := Cartesian(ValsFunction5(P,4,1),[W^-1]);
        elif c = "x ne 0, x~x' if x^6=x'^6 mod p, all z" then 
            vals := Cartesian(ValsFunction5(P,6,1),[0..P-1]);
        elif c = "z=w^i, i=1,2,3,4, x ne 0, x~x' if x^5=x'^5 mod p" then 
            vals := Cartesian(ValsFunction5(P,5,1),
                              List([1..4],X ->(W^X mod P)));
        elif c = "x ne -1,3, See Notes6.114" then 
            vals := ValsFunction20(P);

        # para = [x,t]
        elif c = "all x, t=0,1" then 
            vals := Cartesian([0..P-1], [0,1]);
        elif c = "t=1,-1, x ne 0" then 
            vals := Cartesian([1..P-1], [1, P-1]);
        elif c = "t=1,-1, x ne 0, x~-x" then 
            vals := Cartesian([1..(P-1)/2], [1, P-1]);
        elif c = "all x,t, [x,t]~[t,x]" then 
            vals := UnorderedTuples([0..P-1],2);

        # para = [y,z]
        elif c = "all y, z ne 0, y~-y, z~-z" then 
            vals := Cartesian( [0..(P-1)/2], [1..(P-1)/2]);
        elif c = "all y, z=0,1,w" then 
            vals := Cartesian( [0..P-1], [0,1,W] );
        elif c = "all y,z" then 
            vals := Cartesian( [0..P-1], [0..P-1]);
        elif c = "y,z ne 0, y~-y" then 
            vals := Cartesian( [1..(P-1)/2], [1..P-1]);
        elif c = "y,z ne 0, y~y' if y^4=y'^4 mod p" then 
            vals := Cartesian( ValsFunction5(P, 4, 1), [1..P-1]);
        elif c = "y=0,1, all z" then 
            vals := Cartesian( [0,1], [0..P-1]);
        elif c = "y=0,1, z=0,1,w" then 
            vals := Cartesian( [0,1], [0,1,W] );
        elif c = "y=1,w, all z" then 
            vals := Cartesian( [1,W], [0..P-1] );

        # para = [y,t]
        elif c = "all t, y=0,1" then 
            vals := Cartesian( [0,1], [0..P-1] );
        elif c = "all t,y" then 
            vals := Cartesian( [0..P-1], [0..P-1] );
        elif c = "t=0,1,w, y=0,1" then 
            vals := Cartesian( [0,1], [0,1,W] );
        elif c = "t=1,w, all y" then
            vals := Cartesian( [0..P-1], [1,W] );

        # para = [z,t]
        elif c = "z=0,1,w, t ne 0,1,-1, t~t^-1" then 
            vals := Cartesian( [0,1,W], ValsFunction1(P, [2..P-2] ));
        elif c = "t=1,-1, z=0,1,w" then 
            vals := Cartesian( [0,1,W], [1, P-1] );
        elif c = "all z,t" then 
            vals := Cartesian( [0..P-1], [0..P-1] );
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 3 then 
    
        # sort para
        if x in para and y in para and z in para then 
            para := [x,y,z];
        elif x in para and y in para and t in para then 
            para := [x,y,t];
        elif x in para and z in para and t in para then 
            para := [x,z,t];
        elif y in para and z in para and t in para then 
            para := [y,z,t];
        fi;

        # para = [x,y,z]
        if Length(c) = 0 or c = "all x,y,z" then 
            vals := Tuples( [0..P-1], Length(para) );
        elif c="all x,y, [x,y]~[y,x], z FIXED so that z^2-4 not a square mod p"
        then 
            Z := ValsFunction11(P);
            vals := Concatenation(List([0..P-1], X -> List([0..X],Y->[X,Y,Z])));
        elif c = "all x,y, z=1,w,w^2,w^3,w^4, y~y' if y^5=y'^5 mod p" then 
            vals := Cartesian([0..P-1],ValsFunction5(P,5,0),[1,W,W^2,W^3,W^4]);
        elif c = "x ne 0,all y,z, x~-x" then 
            vals := Cartesian([1..(P-1)/2],[0..P-1],[0..P-1]);
        elif c = "all x,y,z, z~-z" then 
            vals := Cartesian([0..P-1],[0..P-1],[0..(P-1)/2]);
        elif c = "all x,z, y ne 0, y~-y" then 
            vals := Cartesian([0..P-1],[1..(P-1)/2],[0..P-1]);
        elif c = "x ne -2, all y,z, z~-z" then 
            vals := Cartesian(Difference([0..P-1],[P-2]),[0..P-1],[0..(P-1)/2]);
        elif c = "x ne -2w, all y,z, z~-z" then 
            vals := Cartesian(Difference([0..P-1],[(-2*W) mod P]),[0..P-1],
                    [0..(P-1)/2]);
        elif c = "x ne 0, all y,z" then 
            vals := Cartesian([1..P-1],[0..P-1],[0..P-1]);
        elif c = "x ne 0, all y,z, x~-x" then 
            vals := Cartesian([1..(P-1)/2],[0..P-1],[0..P-1]);
        elif c = "x,z ne 0, all y, y~-y, z~-z" then 
            vals := Cartesian([1..P-1],[0..(P-1)/2],[1..(P-1)/2]);
        elif c = "x=0,1, y=0,1, z ne 0,-1" then 
            vals := Cartesian([0,1],[0,1],[1..P-2]);
        elif c = "all x,y,z, [x,y,z]~[-x,-y,-z]" then 
            vals := ValsFunction17(P);
        elif c = "x ne 0, all y,z, x~x' if x^4=x'^4 mod p" then 
            vals := Cartesian(ValsFunction5(P,4,1),[0..P-1],[0..P-1]);
        elif c = "x ne 0, x~x' if x^3=x'^3 mod p, all y,z" then 
            vals := Cartesian(ValsFunction5(P,3,1),[0..P-1],[0..P-1]);
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p, all y,z" then 
            vals := Cartesian(ValsFunction5(P,5,1),[0..P-1],[0..P-1]);
        elif c = "x=w^i, 0 le i lt 2gcd[p-1,3], y ne 0, y~y' if y^3=y'^3 mod p, all z, z~-z" then 
            vals := Cartesian(List([0..2*Gcd(P-1,3)-1], X -> ((W^X) mod P)),
                    ValsFunction5(P,3,1), [0..(P-1)/2]);
        elif c = "y=w^i, i=2,3,4,5, x ne 0, x~x' if x^6=x'^6 mod p, all z" 
        then 
            vals := Cartesian(ValsFunction5(P,6,1),List([2,3,4,5], X -> 
                    W^X mod P), [0..P-1]);
        elif c = "z=w^i, i=1,2,3,4, x ne 0, x~x' if x^5=x'^5 mod p, all y" 
        then
            vals := Cartesian(ValsFunction5(P,5,1), [0..P-1], 
                    List([1,2,3,4], X -> W^X mod P));

        elif c = "See Notes6.150, Case 3" then 
            vals := ValsFunction23(P);
        elif c = "See Notes6.150, Case 4" then 
            vals := ValsFunction23a(P);
        elif c = "See Notes6.150, Case 5" then 
            vals := ValsFunction23b(P);

        # para = [x,y,t]
        elif c = "all t, x=0,1, y=0,1" then 
            vals := Cartesian([0,1],[0,1],[0..P-1]);
        elif c = "t=0,1,w, x=0,1, y=0,1" then
            vals := Cartesian([0,1],[0,1],[0,1,W]);

        # para = [x,z,t]
        elif c = "t ne 0, all x,z, [x,z]~[tz,x/t]" then 
            vals := ValsFunction18(P);
        elif c="t ne 0,1, all x,z such that [x+t][1+z]=1 mod p, [x,z]~[tz,x/t]" 
        then 
            vals := ValsFunction18a(P);

        # para = [y,z,t]
        elif c = "all y,z,t" then 
            vals := Tuples( [0..P-1], Length(para) );
        elif c = "all y,z,t, t~-t" then 
            vals := Cartesian([0..P-1], [0..P-1], [0..(P-1)/2]);
        elif c = "y=0,1, all z,t" then 
            vals := Cartesian([0,1], [0..P-1], [0..P-1]);
        elif c = "y=0,1, z=0,1, t=0,1" then 
            vals := Cartesian([0,1], [0,1], [0,1]);
        elif c = "y=1,w, all z, t ne 0,1,-1, t~t^-1" then
            vals := Cartesian([1,W], [0..P-1], ValsFunction1(P,[2..P-2]));
        elif c = "all y,z,t, (z,t)~(t,z)" then 
            vals := Cartesian([0..P-1],[0..P-1],[0..P-1]);
            vals := Filtered(vals, X -> X[3]<=X[2]); 
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 4 then 
 
        # sort para
        para := [x,y,z,t];

        if Length(c) = 0 or c = "all x,y,z,t" then 
            vals := Cartesian([0..P-1],[0..P-1],[0..P-1],[0..P-1]);
        elif c = "all x,y,z,t, y~-y, if y=0 then t~-t" then 
            vals := Concatenation( 
                      Cartesian([0..P-1],[1..(P-1)/2],[0..P-1],[0..P-1]),
                      Cartesian([0..P-1],[0],[0..P-1],[0..(P-1)/2]));
        elif c = "all y,z,t, x ne 0, x~-x" or
             c = "x ne 0, x~-x, all y,z,t" then 
            vals := Cartesian([1..(P-1)/2],[0..P-1],[0..P-1],[0..P-1]);
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p, all y,z,t" then
            vals := Cartesian(ValsFunction5(P,5,1),[0..P-1],[0..P-1],[0..P-1]);
        elif c = "all x,y,z,t, [x,y,z,t]~[t+1,z+1,y-1,x-1]" then 
            vals := ValsFunction16(P);
        elif d = 6 and c = "See note6.178" then 
            vals := ValsFunction8(P);
        elif d = 7 and c = "See note5.38" and not w in S.param then 
            vals := ValsFunction12(P, 1);
        elif d = 7 and c = "See note5.38" and w in S.param then 
            vals := ValsFunction12(P, 2);
        elif d = 7 and c = "See Notes5.3, Case 4" then 
            vals := ValsFunction8a(P);
        elif c = "See Notes5.12" and w in S.param then 
            vals := ValsFunction22a(P);
        elif c = "See Notes5.12" then 
            vals := ValsFunction22(P);
        elif c = "See Notes6.150, Case 8" then 
            vals := ValsFunction23c(P);
        elif c = "See note2dec5.1" then 
            vals := ValsFunction25(P);
        elif c = "See Notes4.1 Case 4" then 
            vals := ValsFunction27(P);
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 5 then 

        # sort para
        if x in para then
            para := [x,y,z,t,u];
        else
            para := [y,z,t,u,v];
        fi;

        # para = [x,y,z,t,u];
        if Length(c) = 0 then 
            vals := Tuples([0..P-1],5);
        elif c = "See Notes6.163a" then  
            vals := ValsFunction26(P);
        elif c = "See Notes6.163b" then  
            vals := ValsFunction26a(P);

        # para = [y,z,t,u,v];
        elif c = "See Notes6.173" then 
            vals := ValsFunction21(P);
        elif d = 7 and c = "See Notes6.178a" then 
            if fals = 0 then fals := ValsFunction24(P); fi; vals := fals[2]; 
        elif d = 7 and c = "See Notes6.178b" then 
            if fals = 0 then fals := ValsFunction24(P); fi; vals := fals[3]; 
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 6 then 

        # sort para
        para := [x,y,z,t,u,v];

        # para = [x,y,z,t,u,v];
        if Length(c) = 0 then 
            vals := Tuples([0..P-1],6);
        elif c = "See Notes5.14, Case 1" then 
            vals := ValsFunction28(P, 1);
        elif c = "See Notes5.14, Case 2" then 
            vals := ValsFunction28(P, 2);
        elif c = "See Notes5.14, Case 3" then 
            vals := ValsFunction28(P, 3);
        elif c = "See Notes5.14, Case 4" then 
            vals := ValsFunction28(P, 4);
        elif c = "See Notes5.14, Case 5" then 
            vals := ValsFunction28(P, 5);
        elif c = "See Notes5.14, Case 6" then 
            vals := ValsFunction28(P, 6);
        elif c = "See Notes5.14, Case 7" then 
            vals := ValsFunction28(P, 7);
        elif c = "See Notes5.14, Case 8" then 
            vals := ValsFunction28(P, 8);
        elif c = "See Notes5.14, Case 9" then
            vals := ValsFunction28(P, 9);
        elif c = "See Notes5.14, Case 10" then 
            vals := ValsFunction28(P, 10);
        elif c = "See Notes5.14, Case 11" then
            vals := ValsFunction28(P, 11);
        elif c = "See Notes5.14, Case 12" then 
            vals := ValsFunction28(P, 12);
        elif c = "See Notes5.14, Case 13" then 
            vals := ValsFunction28(P, 13);
        elif c = "See Notes5.14, Case 14" then 
            vals := ValsFunction28(P, 14);
        elif c = "See Notes5.14, Case 15" then 
            vals := ValsFunction28(P, 15);
        elif c = "See Notes5.14, Case 16" then 
            vals := ValsFunction28(P, 16);
        elif c = "See Notes5.14, Case 17" then  
            vals := ValsFunction28(P, 17);
        elif c = "See Notes5.14, Case 18" then 
            vals := ValsFunction28(P, 18);
        elif c = "See Notes5.14, Case 19" then 
            vals := ValsFunction28(P, 19);
        elif c = "See Notes5.14, Case 20" then 
            vals := ValsFunction28(P, 20);
        elif c = "See Notes5.14, Case 21" then 
            vals := ValsFunction28(P, 21);
        elif c = "See Notes5.14, Case 22" then 
            vals := ValsFunction28(P, 22);
        elif c = "See Notes5.14, Case 23" then 
            vals := ValsFunction28(P, 23);
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 7 then 

        # sort para
        para := [x,y,z,t,u,v,s];

        # para = [x,y,z,t,u,v,s];
        if c = "See Notes5.14, Case 24" then 
            vals := ValsFunction28(P, 24);
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 8 then 

        # sort para
        para := [x,y,z,t,r,s,u,v];

        # para = [x,y,z,t,r,s,u,v];
        if Length(c) = 0 then 
            vals := Tuples([0..P-1],8);
        elif c = "See Notes5.3, Case 6" then 
            vals := ValsFunction19(P);
        elif c = "See Notes5.3, Case 7" then
            vals := ValsFunction19a(P);
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;

    elif Length(para) = 12 then 

        # sort para
        para := [x,y,z,t,j,k,m,n,r,s,u,v];

        # para = [x,y,z,t,j,k,m,n,r,s,u,v];
        if Length(c) = 0 then 
            vals := Tuples([0..P-1],12);
        elif c = "See Notes4.1, Case 5" then 
            vals := ValsFunction27a(P);
        elif c = "See Notes4.1, Case 6" then
            vals := ValsFunction27b(P);
        else
            Print("parameter ",para," case fails on condition ",c," \n");
        fi;
    else 
        Error("number of parameters is not in list");
    fi;

    if ForAny(vals, x -> Length(x) <> Length(para)) then 
        Error("parameters and values do not match");
    fi;

    # got vals and para
    return List(vals, X -> SpecialiseLiePRingNC( L, P, para, X ));
end;

