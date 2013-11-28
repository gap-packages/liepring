
DensityLPR := function( list )
    local i;
    for i in Reversed([1..Length(list)]) do
        if list[i] = true then return i+1; fi;
    od;
    return 1;
end;

NormedLPR := function( g, l )
    local p, e;

    if g = 0*g then return g; fi;
    p := SCTable(g).prime;
    e := Exponents(g)[l];

    if IsPolynomial(e) or IsRationalFunction(e) then 
        Print("WARNING: Dividing by ",e,"\n");
    fi;
    e := e^-1; 

    if IsInt(p) then 
        e := e mod p; 
    elif e <> 1 and e <> -1 then 
        Print("WARNING: Multiplying by ",e,"\n");
    fi;

    return e*g;
end; 

InsertLPR := function( list, g, k )
    local e, l;
    repeat
        e := Exponents(g);
        l := PositionNonZero(e);
        if l >= k then 
            return false;
        elif list[l] = true then 
            list[l] := NormedLPR( g, l );
            return l;
        else
            g := g - e[l]*list[l];
        fi;
    until false;
end;

StripLPR := function( list )
    local r, d, i, e, j, k;
    r := Length(list);
    d := List(list, x -> PositionNonZero(Exponents(x)));
    for i in [1..r] do
        e := Exponents(list[i]);
        for j in [d[i]+1..r] do
            if e[j] <> 0*e[j] then 
                k := Position( d, j );
                if IsInt(k) then 
                    list[i] := list[i] - e[j]*list[k];
                fi;
            fi;
        od;
    od;
    return list;
end;

IsIntLPR := function( g )
    return ForAll( Exponents(g), IsInt );
end;

BasisByGens := function( L, part, gens )
    local d, p, f, i, a, k, t, s, j, g, b, h; 

    # set up
    if not IsParentLiePRing(L) then L := Parent(L); fi;
    d := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    f := List([1..d], x -> true);

    # fill in part
    for i in [1..Length(part)] do
        a := PositionNonZero(Exponents(part[i]));
        f[a] := part[i];
    od;
    k := DensityLPR(f);

    # init
    t := Filtered( gens, x -> x <> Zero(L) and IsIntLPR(x) );
    s := Filtered( gens, x -> x <> Zero(L) and not IsIntLPR(x) );
    i := 1;
    j := Length(t);

    # step 1: process t (integral elements)
    while i <= j do
        g := t[i];
        a := InsertLPR( f, g, k );
        if IsInt(a) then 

            # reset density
            if a = k-1 then k := DensityLPR( f ); fi;

            # add powers 
            b := p*f[a]; 
            if b <> Zero(L) then 
                if IsIntLPR(b) then 
                    Add(t, b); j := j+1; 
                else 
                    Add(s, b);
                fi;
            fi;

            for h in [1..d] do
                if f[h] <> true and h <> a then  
                    b := f[h]*f[a];
                    if b <> Zero(L) then 
                        if IsIntLPR(b) then 
                            Add(t, b); j := j+1; 
                        else 
                            Add(s, b);
                        fi;
                    fi;
                fi;
            od;
        fi;
        i := i + 1;
    od;

    # step 2: process s (non-integral elements)
    i := 1;
    j := Length(s);
    while i <= j do
        g := s[i];
        a := InsertLPR( f, g, k );
        if IsInt(a) then 

            # reset density
            if a = k-1 then k := DensityLPR( f ); fi;

            # add powers 
            b := p*f[a]; 
            if b <> Zero(L) then Add(s, b); j := j+1; fi;

            for h in [1..d] do
                if f[h] <> true and h <> a then  
                    b := f[h]*f[a];
                    if b <> Zero(L) then Add(s, b); j := j+1; fi;
                fi;
            od;
        fi;
        i := i + 1;
    od;

    # strip and return
    f := Filtered(f, x -> x <> true);
    return StripLPR(f);
end;

LiePSubringByBasis := function( L, basis )
    local U, K;

    # get basis and parent
    K := Parent(L);
    if Length(basis) = DimensionOfLiePRing(K) then return K; fi;

    # compute
    if Length(basis) = 0 then 
        U := TrivialSubalgebra(L);
    else
        U := RingByGenerators( basis );
    fi;

    # add info
    SetBasisOfLiePRing(U, basis);
    SetDimensionOfLiePRing(U, Length(basis));
    SetPrimeOfLiePRing(U, PrimeOfLiePRing(K));
    SetParent(U, K);
    SetIsLiePRing(U, true);
    SetIsParentLiePRing(U, false);

    # return
    return U;
end;

LiePSubring := function( L, gens )
    return LiePSubringByBasis( L, BasisByGens(L, [], gens) );
end;

LiePClosure := function( L, U, gens )
    return LiePSubringByBasis( L, BasisByGens(L, BasisOfLiePRing(U), gens) );
end;

LiePIdeal := function( L, gens )
    local K, b, w, v, c;
    K := Parent(L);
    b := BasisByGens( K, [], gens );
    w := BasisOfLiePRing(L);
    repeat
        v := Flat(List( b, x -> List(w, y -> x*y) ));
        c := BasisByGens( K, b, v );
        if Length(c) = DimensionOfLiePRing(K) then return K; fi;
        if Length(c) = Length(b) then return LiePSubringByBasis( K, b ); fi;
        b := c;
    until false;
end;

LiePCommutator := function( U, V )
    local bU, bV, bC;
    bU := BasisOfLiePRing(U);
    bV := BasisOfLiePRing(V);
    bC := Flat(List(bU, x -> List(bV, y -> x*y)));
    return LiePSubring(U, bC);
end;

LiePPCommutator := function( U, L )
    local C, g;
    C := LiePCommutator( U, L );
    g := List(BasisOfLiePRing(U), x -> PrimeOfLiePRing(U)*x);
    return LiePClosure(L, C, g);
end;

LiePRump := function( L )
    return LiePPCommutator(L,L);
end;

LiePMinimalGeneratingSet := function( L )
    local U, bL, bU, eL, eU, v;
    U := LiePRump(L);
    bL := BasisOfLiePRing(L);
    bU := BasisOfLiePRing(U);
    eU := List(bU, Exponents);
    eL := List(bL, Exponents);
    v := BaseSteinitzVectors(eL, eU).factorspace;
    return List(v, x -> x*BasisOfLiePRing(Parent(L)));
end;

LiePLowerPCentralSeries := function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePPCommutator( s[Length(s)], L );
        Add(s, U);
    od;
    return s;
end;

LiePLowerCentralSeries := function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePCommutator( L, s[Length(s)] );
        Add(s, U);
    od;
    return s;
end;

LiePDerivedSeries := function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePCommutator( s[Length(s)], s[Length(s)] );
        Add(s, U);
    od;
    return s;
end;

LiePIsIdeal := function(L, U)
    local bL, bU, a, b;
    bL := BasisOfLiePRing(L);
    bU := BasisOfLiePRing(U);
    if ForAny(bU, x -> not x in L) then return false; fi;
    for a in bL do
        for b in bU do
            if not a*b in U then return false; fi;
        od;
    od;
    return true;
end;

LiePQuotientNC := function(L, U)
    local K, p, bL, eL, bU, eU, eQ, bQ, T, i, j, e;

    # check for parameters - this case is not supported
    K := Parent(L);
    if IsBound(SCTable(Zero(K)).param) then return fail; fi;
    p := PrimeOfLiePRing(K);

    # proceed
    bL := BasisOfLiePRing(L);
    eL := List(bL, Exponents);
    bU := BasisOfLiePRing(U);
    eU := List(bU, Exponents);
    eQ := BaseSteinitzVectors(eL, eU).factorspace;
    bQ := List(eQ, x -> x*bL);

    eL := Concatenation( eQ, eU );

    # set up new table
    T := rec( dim := Length(bQ),
              prime := PrimeOfLiePRing(U),
              tab := [],
              param := [] );

    for i in [1..T.dim] do
        for j in [1..i-1] do
            e := Exponents( bQ[i]*bQ[j] );
            e := SolutionIntMat( eL, e ){[1..T.dim]};
            Add(T.tab, WordByExps(e));
        od;
        e := Exponents( p*bQ[i] );
        e := SolutionIntMat( eL, e ){[1..T.dim]};
        Add(T.tab, WordByExps(e));
    od;
    return LiePRingBySCTableNC(T);
end;

LiePQuotient := function(L,U)
    if not LiePIsIdeal(L,U) then return false; fi;
    return LiePQuotientNC(L,U);
end;
   
