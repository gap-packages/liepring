
#    [ dimension, m, c, 
#      short presentation, 
#      structure constants of pcp,
#      conditions,,name]

FindIndets := function( list )
    local l, r, i, e, j, k;
    l := Flat(list);
    if ForAll(l, x -> IsInt(x)) then return []; fi;
    r := [];
    for i in [1..Length(l)] do
        if IsPolynomial(l[i]) then 
            e := ExtRepPolynomialRatFun(l[i]);
            for j in [1,3..Length(e)-1] do
                for k in [1,3..Length(e[j])-1] do
                    AddSet(r, e[j][k]);
                od;
            od;
        fi;
    od;
    return List(r, x -> Indeterminate(Integers, x));
end;

LiePRingByData := function( dim, l )
    local r, S, L, cond;

    # create a Lie p-ring
    r := FindIndets( l[5] );
    S := rec( dim := dim, prime := Indeterminate(Integers,"p"), 
              tab := l[5], param := r );
    L := CreateLiePRing(S);
    SetIsLiePRing(L, true);

    # generate some info
    cond := [];
    if IsBound(l[6]) then cond[1]  := l[6]; else cond[1] := ""; fi;
    if IsBound(l[7]) then cond[2]  := l[7]; else cond[2] := ""; fi;

    # add some attributes
    SetPClassOfLiePRing( L, l[3] );
    SetMinimalGeneratorNumberOfLiePRing( L, l[2] );
    SetDimensionOfLiePRing( L, l[1] );
    SetLibraryName( L, l[8] );
    SetShortPresentation( L, l[4] );
    SetLibraryConditions( L, cond );
    
    return L;
end;

NumberOfLiePRings := function( arg )
    local dim, P;
    dim := arg[1];

    # if a prime is given
    if Length(arg) = 2 then 
        P := arg[2];
        if not IsPrimeInt(P) or P=2 then return fail; fi;
        if P > 5 then 
            return NumberSmallGroups( arg[2]^dim );
        elif P = 5 then 
            return 34317;
        elif P = 3 then 
            return 1566;
        fi;
    fi;

    # otherwise return the number of generic Lie $p$-rings
    if dim = 1 then 
        return 1;
    elif dim = 2 then 
        return 2;
    elif dim = 3 then 
        return 5;
    elif dim = 4 then 
        return 15;
    elif dim = 5 then 
        return 75;
    elif dim = 6 then 
        return 542;
    elif dim = 7 then 
        return 4773;
    fi;
end;

##
## dim
## dim, P
## dim, gen, cl
## dim, P, gen, cl
##
InstallGlobalFunction( LiePRingsByLibrary, function( arg )
    local dim, lie, spe;

    # get the dimension and the Lie p-rings
    dim := arg[1];
    lie := List([1..NumberOfLiePRings(dim)], 
                 n -> LiePRingByData( dim, LIE_DATA[dim][n] ) );

    # consider cases
    if Length(arg) = 1 then 
        return lie; 
    elif Length(arg) = 2 then  
        spe := List( lie, L -> LiePRingsInFamily( L, arg[2] ) );
        spe := Filtered(spe, l -> l <> fail );
        return Flat(spe);
    elif Length(arg) = 3 then 
        lie := Filtered(lie, x -> MinimalGeneratorNumberOfLiePRing(x)=arg[2]);
        lie := Filtered(lie, x -> PClassOfLiePRing(x)=arg[3]);
        return lie;
    elif Length(arg) = 4 then 
        lie := Filtered(lie, x -> MinimalGeneratorNumberOfLiePRing(x)=arg[3]);
        lie := Filtered(lie, x -> PClassOfLiePRing(x)=arg[4]);
        spe := List( lie, L -> LiePRingsInFamily( L, arg[2] ) );
        spe := Filtered(spe, l -> l <> fail );
        return Flat(spe);
    fi;

end ); 

    
