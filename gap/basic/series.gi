
BindGlobal( "LiePCommutator", function( U, V )
    local bU, bV, bC, sr;
    bU := BasisOfLiePRing(U);
    bV := BasisOfLiePRing(V);
    bC := Flat(List(bU, x -> List(bV, y -> x*y)));
    return LiePRecSubring(Parent(U), bC);
end );

BindGlobal( "LiePPCommutator", function( U, V )
    local bU, bV, b, g;
    bU := BasisOfLiePRing(U);
    bV := BasisOfLiePRing(V);
    b := Flat(List(bU, x -> List(bV, y -> x*y)));
    g := List(BasisOfLiePRing(V), x -> PrimeOfLiePRing(V)*x);
    return LiePRecSubring(Parent(U), Concatenation(b,g));
end );

BindGlobal( "LiePRump", function( L )
    return LiePPCommutator(L,L);
end );

BindGlobal( "LiePMinimalGeneratingSet", function( L )
    local U, bL, bU, eL, eU, v;
    U := LiePRump(L);
    bL := BasisOfLiePRing(L);
    bU := BasisOfLiePRing(U);
    eU := List(bU, Exponents);
    eL := List(bL, Exponents);
    v := BaseSteinitzVectors(eL, eU).factorspace;
    return List(v, x -> x*BasisOfLiePRing(Parent(L)));
end );

BindGlobal( "LiePLowerPCentralSeries", function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePPCommutator(L, s[Length(s)]);
        if Length(U) = 1 then Add(s, U[1]); else return fail; fi;
    od;
    return s;
end );

BindGlobal( "LiePLowerCentralSeries", function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePCommutator(L, s[Length(s)]);
        if Length(U) = 1 then Add(s, U[1]); else return fail; fi;
    od;
    return s;
end );

BindGlobal( "LiePDerivedSeries", function( L )
    local s, U;
    s := [L];
    while DimensionOfLiePRing(s[Length(s)]) > 0 do
        U := LiePCommutator( s[Length(s)], s[Length(s)] );
        if Length(U) = 1 then Add(s, U[1]); else return fail; fi;
    od;
    return s;
end );


