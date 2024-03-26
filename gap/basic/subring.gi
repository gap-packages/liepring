
##
## g is an element of L in reduced form c1 b1 + ... + cn bn
## l is the depth of g, that is, cl <> 0
## cl is an element of Q[x_1, ..., x_m, w]
##
BindGlobal( "NormedLPR", function( L, g, l )
    local p, e, f;

    if g = 0*g then return g; fi;
    p := SCTable(g).prime;
    e := Exponents(g)[l];

    if e = e^0 then 
        return g;
    elif IsRat(e) and IsInt(p) then 
        f := e^-1 mod p;
    elif IsRat(e) or IsLiePUnit(SCTable(Zero(L)).ring.units, e) then  
        f := e^-1;
    else
        return fail;
    fi;
    return f*g;
end ); 

BindGlobal( "DensityLPR", function( list )
    local i;
    for i in Reversed([1..Length(list)]) do
        if list[i] = true then return i+1; fi;
    od;
    return 1;
end );

BindGlobal( "InsertLPR", function( L, list, g, k )
    local e, l, a;
    repeat
        e := Exponents(g);
        l := DepthVector(e);
        if l >= k then 
            return true;
        elif list[l] <> true then 
            g := g - e[l]*list[l];
        else
            a := NormedLPR(L, g, l);
            if a = fail then 
                return [l, e[l]];
            else
                list[l] := a; return l;
            fi;
        fi;
    until false;
end );

BindGlobal( "StripLPR", function( list )
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
end );

BindGlobal( "IsIntLPR", function( g )
    local e, i;
    e := Exponents(g);
    for i in [1..Length(e)] do
        if not IsInt(e[i]) then return false; fi;
    od;
    return true;
end );

BindGlobal( "BasisByGens", function( L, part, gens )
    local d, p, f, i, a, t, s, g, b, h, k;

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

    # sort nicely
    t := Filtered( gens, x -> x <> Zero(L) and IsIntLPR(x) );
    s := Filtered( gens, x -> x <> Zero(L) and not IsIntLPR(x) );
    t := Concatenation(t,s);

    # loop
    for g in t do
        a := InsertLPR( L, f, g, k );
        if IsInt(a) then 

            # reset density
            if a = k-1 then k := DensityLPR( f ); fi;

            # close under powers
            b := p*f[a]; 
            if b <> Zero(L) then Add(t, b); fi;

            # close under mult
            for h in [1..d] do
                if f[h] <> true and h <> a then  
                    b := f[h]*f[a];
                    if b <> Zero(L) then Add(t, b); fi;
                fi;
            od;
        elif IsList(a) then 
            return a[2];
        fi;
    od;

    # strip and return
    f := Filtered(f, x -> x <> true);
    return StripLPR(f);
end );

BindGlobal( "LiePSubringByBasis", function( L, basis )
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
    SetBasisOfLiePRing(U, basis);
    SetDimensionOfLiePRing(U, Length(basis));
    SetPrimeOfLiePRing(U, PrimeOfLiePRing(L));
    SetParent(U, Parent(L));
    SetIsLiePRing(U, true);
    SetIsParentLiePRing(U, false);

    # return
    return U;
end );

InstallMethod( BasisOfLiePRing, true, [IsLiePRing], 0, function(L)
    if IsParentLiePRing(L) then return GeneratorsOfRing(L); fi;
    return BasisByGens( Parent(L), GeneratorsOfRing(L) );
end );

BindGlobal( "LiePSubring", function( L, gens )
    local b;
    b := BasisByGens(L, [], gens);
    if not IsList(b) then return fail; fi;
    return LiePSubringByBasis( L, b );
end );

BindGlobal( "LiePClosure", function( L, U, gens )
    local b;
    b := BasisByGens(L, BasisOfLiePRing(U), gens);
    if not IsList(b) then return fail; fi;
    return LiePSubringByBasis(L, b);
end );

BindGlobal( "LiePRecSubring", function( arg )
    local L, gens, base, T, R, 
          t, b, A, U, Z, U1, L1, B1, g1, Z1, L2, B2, g2, b1, b2;

    # get arguments
    L := arg[1];
    gens := arg[2];
    if Length(arg) = 3 then base := arg[3]; else base := []; fi;

    # init
    T := [[L,gens, base]];
    R := [];

    # toop
    while Length(T) > 0 do
        t := T[Length(T)];
        Unbind(T[Length(T)]);
        b := BasisByGens( t[1], t[3], t[2] );
        if IsList(b) then 
            Add( R, LiePSubringByBasis(t[1], b) );

        elif not IsBool(b) then 

            # get info
            U := ShallowCopy(SCTable(Zero(L)).ring.units);
            Z := ShallowCopy(SCTable(Zero(L)).ring.zeros);

            # case 1
            U1 := Concatenation(U, [b]);
            L1 := LiePRingCopy(L, U1, Z);
            B1 := BasisOfLiePRing(L1);
            g1 := List(gens, x -> LiePImageByBasis(B1, x));
            b1 := List(base, x -> LiePImageByBasis(B1, x));
            Add(T, [L1, g1, b1]);

            # case 2
            Z1 := Concatenation(Z, [b]);
            L2 := LiePRingCopy(L, U, Z1);
            B2 := BasisOfLiePRing(L2);
            g2 := List(gens, x -> LiePImageByBasis(B2, x));
            b2 := List(base, x -> LiePImageByBasis(B2, x));
            Add(T, [L2, g2, b2]);
        fi;
    od;

    return R;
end );

BindGlobal( "LiePIdeal", function( L, gens )
    local K, b, w, v, c;
    b := BasisByGens( L, [], gens );
    if not IsList(b) then return fail; fi;
    w := BasisOfLiePRing(L);
    repeat
        v := Flat(List( b, x -> List(w, y -> x*y) ));
        c := BasisByGens( L, b, v );
        if not IsList(c) then return fail; fi;
        if Length(c) = DimensionOfLiePRing(L) then return L; fi;
        if Length(c) = Length(b) then return LiePSubringByBasis(L,b); fi;
        b := c;
    until false;
end );

BindGlobal( "LiePIsIdeal", function(L, U)
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
end );

BindGlobal( "LiePQuotientByTable", function( T, U )
    local V, b, c, Q, i, j, e;

    # the trivial case
    if Length(U) = 0 then return LiePRingBySCTableNC(T); fi;

    # get factor
    V := FactorSpace@(Length(U[1]), U);
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
            Add(Q.tab, WordByExps@(e));
        od;
        e := LRReduceExp( T, T.prime*V[i] );
        e := e*c;
        e := e{[1..Q.dim]};
        Add(Q.tab, WordByExps@(e));
    od;
    return LiePRingBySCTableNC(Q);
end );

BindGlobal( "LiePQuotientNC", function(L, U)
    local S, u;
    S := SCTable(Zero(L));
    u := List(BasisOfLiePRing(U), Exponents);
    return LiePQuotientByTable(S, u);
end );

BindGlobal( "LiePQuotient", function(L,U)
    if not LiePIsIdeal(L,U) then return false; fi;
    return LiePQuotientNC(L,U);
end );
   
