
DensityLPR := function( list )
    local i;
    for i in Reversed([1..Length(list)]) do
        if list[i] = true then return i+1; fi;
    od;
    return 1;
end;

NormedLPR := function( L, g, l )
    local p, e, f;

    if g = 0*g then return g; fi;
    p := SCTable(g).prime;
    e := Exponents(g)[l];

    if IsInt(e) and IsInt(p) then 
        f := e^-1 mod p;
        return f*g;
    elif IsInt(e) and AbsInt(e) in [1,2,3,1/2,1/3]  then  
        f := e^-1;
        return f*g;
    elif IsLaurentPolynomial(e) and IsRootPower(e) then 
        f := e^-1;
        return f*g;
    elif IsBound(L!.inv) and ((e in L!.inv) or (-e in L!.inv)) then 
        f := e^-1;
        return f*g;
    else
        Print("WARNING: Dividing by ",e,"\n");
        return fail;
    fi;
end; 

MyDepth := function(vec)
    local i;
    for i in [1..Length(vec)] do
        if vec[i] <> 0*vec[i] then return i; fi;
    od;
    return Length(vec)+1;
end;

InsertLPR := function( L, list, g, k )
    local e, l;
    repeat
        e := Exponents(g);
        l := MyDepth(e);
        if l >= k then 
            return false;
        elif list[l] = true then 
            list[l] := NormedLPR( L, g, l );
            if list[l] = fail then return fail; fi;
            return l;
        else
            g := g - e[l]*list[l];
        fi;
    until false;
end;

StripLPR := function( list )
    local r, d, i, e, j, k, f;
    r := Length(list);
    for i in [1..r] do
        e := Exponents(list[i]);
        d := DepthVector(e);
        for j in [1..i-1] do
            f := Exponents(list[j]);
            if f[d] <> 0*f[d] then 
                list[j] := list[j] - f[d]*list[i];
            fi;
        od;
    od;
    return list;
end;

IsIntLPR := function( g )
    local e, i;
    e := Exponents(g);
    for i in [1..Length(e)] do
        if not IsInt(e[i]) and not IsRootPower(e[i]) then return false; fi;
    od;
    return true;
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
        a := DepthVector(Exponents(part[i]));
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
        a := InsertLPR( L, f, g, k );
        if a = fail then 
            Error(" should not happen ");
            return fail;
        elif IsInt(a) then 

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
        a := InsertLPR( L, f, g, k );
        if a = fail then 
            return fail;
        elif IsInt(a) then 

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
    local U;

    # get basis and parent
    if Length(basis) = DimensionOfLiePRing(L) then return L; fi;

    # compute
    if Length(basis) = 0 then 
        U := TrivialSubalgebra(L);
    else
        U := RingByGenerators( basis );
    fi;

    # add info
    if IsBound(L!.inv) then U!.inv := L!.inv; fi;
    SetBasisOfLiePRing(U, basis);
    SetDimensionOfLiePRing(U, Length(basis));
    SetPrimeOfLiePRing(U, PrimeOfLiePRing(L));
    SetParent(U, Parent(L));
    SetIsLiePRing(U, true);
    SetIsParentLiePRing(U, false);

    # return
    return U;
end;

LiePSubring := function( L, gens )
    local b;
    b := BasisByGens(L, [], gens);
    if b = fail then return fail; fi;
    return LiePSubringByBasis( L, b );
end;

LiePClosure := function( L, U, gens )
    local b;
    b := BasisByGens(L, BasisOfLiePRing(U), gens);
    if b = fail then return fail; fi;
    return LiePSubringByBasis( L, b);
end;

LiePIdeal := function( L, gens )
    local K, b, w, v, c;
    b := BasisByGens( L, [], gens );
    if b = fail then return fail; fi;
    w := BasisOfLiePRing(L);
    repeat
        v := Flat(List( b, x -> List(w, y -> x*y) ));
        c := BasisByGens( L, b, v );
        if c = fail then return fail; fi;
        if Length(c) = DimensionOfLiePRing(L) then return L; fi;
        if Length(c) = Length(b) then return LiePSubringByBasis( L, b ); fi;
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
        if U = fail then return fail; fi;
        Add(s, U);
    od;
    return s;
end;

LiePLowerCentralSeries := function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePCommutator( L, s[Length(s)] );
        if U = fail then return fail; fi;
        Add(s, U);
    od;
    return s;
end;

LiePDerivedSeries := function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePCommutator( s[Length(s)], s[Length(s)] );
        if U = fail then return fail; fi;
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

LiePQuotientByTable := function( T, U )
    local V, b, c, Q, i, j, e;

    # the trivial case
    if Length(U) = 0 then return LiePRingBySCTableNC(T); fi;

    # get factor
    V := FactorSpace(Length(U[1]), U);
    b := Concatenation(V,U);
    b := b*IndeterminateByName("w")^0;
    c := MakeInt(b^-1);

    # set up new table
    Q := rec( dim := Length(V), prime := T.prime, tab := [], param := []);
    if IsBound(T.param) then Q.param := T.param; fi;
    for i in [1..Length(V)] do
        for j in [1..i-1] do
            e := LRMultiply( T, V[i], V[j] );
            e := e*c;
            e := e{[1..Q.dim]};
            Add(Q.tab, WordByExps(e));
        od;
        e := LRReduceExp( T, T.prime*V[i] );
        e := e*c;
        e := e{[1..Q.dim]};
        Add(Q.tab, WordByExps(e));
    od;
    return LiePRingBySCTableNC(Q);
end;

LiePQuotientNC := function(L, U)
    local S, u;
    S := SCTable(Zero(L));
    u := List(BasisOfLiePRing(U), Exponents);
    return LiePQuotientByTable(S, u);
end;

LiePQuotient := function(L,U)
    if not LiePIsIdeal(L,U) then return false; fi;
    return LiePQuotientNC(L,U);
end;
   
