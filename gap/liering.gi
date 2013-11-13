
InstallGlobalFunction( CreateLiePRing, function(SC) 
    local fam, I, b, R;

    # check
    if Length(SC.param) = 0 then Unbind(SC.param); fi;

    # create family with type 
    fam := NewFamily( "LPRElementsFamily", IsLPRElement, IsLPRElement );
    fam!.type := NewType( fam, IsLPRElementRep );

    # add coefficientsfamily
    SetCoefficientsFamily(fam, ElementsFamily( FamilyObj( Integers )));
    SC.fam := fam;

    # create basis and ring
    if SC.dim > 0 then 
        I := IdentityMat(SC.dim);
        b := List(I, x -> LPRElementByExponentsNC( SC, x ) );
        R := RingByGenerators( b );
        SetBasisOfLiePRing(R, b);
    else
        b := [LPRElementByExponentsNC(SC, [])];
        R := RingByGenerators( b );
        SetBasisOfLiePRing(R, []);
    fi;

    # try for size
    if IsInt(SC.prime) then SetSize(R, SC.prime^SC.dim); fi;
    SetDimensionOfLiePRing(R, SC.dim);
    SetPrimeOfLiePRing(R, SC.prime);
    SetIsParentLiePRing(R, true);

    # return
    return R;
end );

InstallGlobalFunction( LiePRingBySCTableNC, function( SC )
    local R;
    R := CreateLiePRing(SC);
    SetIsLiePRing(R, true);
    return R;
end );

InstallGlobalFunction( LiePRingBySCTable, function( SC )
    local R;
    if SC = fail then return fail; fi;
    R := CreateLiePRing(SC);
    if not IsLiePRing(R) then return fail; fi;
    return R;
end );

CheckIsLiePRing := function(L)
    local l, p, d, a, i, j, k, c;

    # get gens and table
    l := GeneratorsOfRing(L);
    p := PrimeOfLiePRing(L);
    d := DimensionOfLiePRing(L);

    # precompute pairs
    a := List([1..d], x -> List([1..d], y -> l[x]*l[y]));

    # test biadditivity
    for i in [1..d] do
        if (p*l[i])*l[i] <> Zero(L) then 
            Print("p* ",i," \n");
            return false;
        fi;
        for j in [1..d] do
            c := p*a[i][j];
            if c <> l[i]*(p*l[j]) then 
                Print(i," * (p * ",j,")\n");
                return false;
            fi;
            if c <> (p*l[i])*l[j] then 
                Print("(p * ",i,") * ",j," \n");
                return false;
            fi;
        od;
    od;

    # test Jacobi
    for i in [1..d] do
        for j in [1..d] do
            for k in [1..d] do
                if l[i]*a[j][k]+l[j]*a[k][i]+l[k]*a[i][j] <> Zero(L) then 
                    Print(i,"*",j,"*",k,"\n");
                    return false; 
                fi;
            od;
        od;
    od;
    
    return true;
end;

InstallMethod( IsLiePRing, true, [IsRing], 0, function(L)
    return CheckIsLiePRing(L);
end );

InstallMethod( PrintObj, true, [IsLiePRing], SUM_FLAGS, function(L)
    Print("<Lie ring with ",Length(GeneratorsOfRing(L))," generators> \n");
end );

InstallMethod( ViewObj, true, [IsLiePRing], SUM_FLAGS, function(L)
    if IsBound(SCTable(Zero(L)).param) then 
        Print("<Lie ring of dimension ",DimensionOfLiePRing(L),
              " over prime ",PrimeOfLiePRing(L),
              " with parameters ",SCTable(Zero(L)).param,">");
    else
        Print("<Lie ring of dimension ",DimensionOfLiePRing(L),
              " over prime ",PrimeOfLiePRing(L),">");
    fi;
end );

PrintExp := function( n, e )
    local d, i;
    d := First([1..n], x -> e[x] <> 0*e[x]);
    if d = fail then Print(0,"\n"); return; fi;

    if e[d] = e[d]^0 then 
        Print("l",d);
    elif e[d] = -e[d]^0 then 
        Print("-l",d);
    else
        Print(e[d],"*l",d);
    fi;
    for i in [d+1..n] do
        if e[i] <> 0*e[i] then 
            if e[i] = e[i]^0 then 
                Print(" + l",i);
            elif e[i] = -e[i]^0 then 
                Print(" - l",i);
            else
                Print(" + ",e[i],"*l",i);
            fi;
        fi;
    od;
    Print("\n");
end;

InstallGlobalFunction(ViewShortPresentation, function(L)
    if not IsBound(L!.ShortPresentation) then return; fi;
    Print(L!.ShortPresentation);
end );

InstallGlobalFunction(ViewPCPresentation, function(L)
    local S, n, p, l, i, e, j;
    S := SCTable(Zero(L));
    n := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    l := BasisOfLiePRing(L);

    for i in [1..n] do
        e := ExponentsLPR(L, p*l[i]);
        if e <> 0*e then 
           Print(p,"*l",i," = "); PrintExp(n,e);
        fi;
    od;

    for i in [1..n] do
        for j in [1..i-1] do
            e := ExponentsLPR(L, l[i]*l[j]);
            if e <> 0*e then 
               Print("[l",i,",l",j,"] = "); PrintExp(n,e);
            fi;
        od;
    od;
end ); 

InstallGlobalFunction( ParametersOfLiePRing, function(L)
    local S;
    S := SCTable(Zero(L));
    if not IsBound(S.param) then 
        return [];
    else
        return S.param;
    fi;
end );

InstallMethod( PrimeOfLiePRing, true, [IsLiePRing], 0, function(L)
    return SCTable(Zero(L)).prime;
end );

InstallMethod( BasisOfLiePRing, true, [IsLiePRing], 0, function(L)
    if IsParentLiePRing(L) then return GeneratorsOfRing(L); fi;
    return BasisByGens( Parent(L), GeneratorsOfRing(L) );
end );

InstallMethod( DimensionOfLiePRing, true, [IsLiePRing], 0, function(L)
    return Length(BasisOfLiePRing(L));
end );


