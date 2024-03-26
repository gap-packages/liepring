
InstallGlobalFunction( LPRElementConstruction, function( SC, list )
    local i, elm;
    for i in [1..SC.dim] do
        if IsPolynomial(list[i]) and IsConstantRationalFunction(list[i]) then 
            list[i] := Value(list[i], 1);
        fi;
    od;
    elm := rec( sctable := SC, exponents := Immutable(list), name := "l" );
    return Objectify( SC.fam!.type, elm );
end );

InstallGlobalFunction( LPRElementByExponentsNC, function( SC, list )
    return LPRElementConstruction( SC, list );
end );

InstallGlobalFunction( LPRElementByExponents, function( SC, list )
    return LPRElementConstruction( SC, LRReduceExp( SC, list ) );
end );

InstallGlobalFunction( LPRElementByWordList, function( SC, list )
    return LPRElementConstruction( SC, LRCollectWord( SC, list ) );
end );

InstallMethod( SCTable, [ IsLPRElementRep ], l -> l!.sctable );
InstallMethod( Exponents, [ IsLPRElementRep ], l -> l!.exponents );
InstallMethod( NameTag, [ IsLPRElementRep ], l -> l!.name );

InstallMethod( PrintObj, true, [ IsLPRElement ], 0, function(elm)
    local l, e, n, o, z, i, j;
    l := NameTag(elm);
    e := Exponents(elm);
    n := Length(e);
    i := First([1..n], x -> e[x] <> 0*e[x]);

    # the trivial case
    if i = fail then Print("0"); return; fi;

    # the first entry
    if e[i] = e[i]^0 then 
        Print(Concatenation(l,String(i)));
    else
        Print(e[i],"*",Concatenation(l,String(i)));
    fi;

    # the other entries
    for j in [i+1..n] do
        if e[j] = e[j]^0 then
            Print(" + ",Concatenation(l,String(j)));
        elif e[j] <> 0*e[j] then
            Print(" + ",e[j],"*",Concatenation(l,String(j)));
        fi;
    od;
end);

InstallMethod( \=, IsIdenticalObj, [IsLPRElement, IsLPRElement], 0,
function( elm1, elm2)
    return Exponents(elm1) = Exponents(elm2);
end );

InstallMethod( \<, IsIdenticalObj, [IsLPRElement, IsLPRElement], 0,
function( elm1, elm2 )
    return Exponents(elm1) < Exponents(elm2);
end );

InstallMethod( \+, IsIdenticalObj, [IsLPRElement, IsLPRElement], 0,
function( elm1, elm2 )
    return LPRElementByExponents( SCTable(elm1), 
                                 Exponents(elm1)+Exponents(elm2) );
end );
    
InstallMethod( \-, IsIdenticalObj, [IsLPRElement, IsLPRElement], 0,
function( elm1, elm2 )
    return LPRElementByExponents( SCTable(elm1), 
                                 Exponents(elm1)-Exponents(elm2) );
end );
    
InstallMethod( \*, true, [IsRat, IsLPRElement], 0,
function( a, elm )
    return LPRElementByExponents( SCTable(elm), a*Exponents(elm) );
end );

InstallMethod( \*, true, [IsRationalFunction, IsLPRElement], 0,
function( a, elm )
    return LPRElementByExponents( SCTable(elm), a*Exponents(elm) );
end );

InstallMethod( \*, IsIdenticalObj, [IsLPRElement, IsLPRElement], 0,
function( elm1, elm2 )
    return LPRElementByExponents( SCTable(elm1),
       LRMultiply(SCTable(elm1), Exponents(elm1), Exponents(elm2)));
end );

InstallMethod( ZeroOp, [IsLPRElement], 0, function(elm)
    return LPRElementByExponents( SCTable(elm), 0*Exponents(elm) );
end );

InstallMethod( AdditiveInverseOp, [IsLPRElement], 0, function(elm)
    return LPRElementByExponents( SCTable(elm), -Exponents(elm) );
end );

BindGlobal( "ExponentsByBasis", function( list, elm )
    local dep, exp, d, j;
    if Length(list)=0 and elm = 0*elm then return []; fi;
    if Length(list)=0 and elm <> 0*elm then return fail; fi;
    dep := List(list, x -> DepthVector(Exponents(x)));
    exp := Exponents(elm);
    while exp <> 0*exp do
        d := DepthVector(exp);
        j := Position(dep, d);
        if j = fail then return fail; fi;
        elm := elm - exp[d]*list[j];
        exp := Exponents(elm);
    od;
    return exp;
end );

InstallMethod( \in, true, [IsLPRElement, IsLiePRing], 0,
function(elm, L)
    if FamilyObj(Zero(L)) <> FamilyObj(elm) then return false; fi;
    if IsParentLiePRing(L) then return true; fi;
    return ExponentsByBasis( BasisOfLiePRing(L), elm ) <> fail;
end );
    
BindGlobal( "ExponentsLPR", function( L, elm )
    if IsParentLiePRing(L) then return Exponents(elm); fi;
    return ExponentsByBasis(BasisOfLiePRing(L), elm);
end );

BindGlobal( "LiePOrder", function( elm )
    local p, i;
    p := SCTable(elm).prime;
    i := 0;
    repeat
        if p^i*elm = 0*elm then return p^i; fi;
        i := i+1;
    until false;
end );

