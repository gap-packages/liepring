
#    [ dimension, m, c, 
#      short presentation, 
#      structure constants of pcp,
#      conditions,,name]

BindGlobal( "LiePRingByData", function( dim, l )
    local r, S, L, cond, i;

    # generate some info
    cond := [];
    if IsBound(l[6]) then cond[1] := l[6]; else cond[1] := ""; fi;
    if IsBound(l[7]) then cond[2] := l[7]; else cond[2] := ""; fi;
    if IsBound(l[9]) then cond[3] := l[9]; fi;

    # create a Lie p-ring
    r := VarsOfSCTab( l[5] );
    i := RingInvariantsByData( cond[1], l[1], r ); 
    S := rec( dim := dim, prime := Indeterminate(Integers,"p"), 
              tab := l[5], param := r, ring := i );
    L := CreateLiePRing(S);
    SetIsLiePRing(L, true);

    # add some attributes
    SetPClassOfLiePRing( L, l[3] );
    SetMinimalGeneratorNumberOfLiePRing( L, l[2] );
    SetDimensionOfLiePRing( L, l[1] );
    SetLibraryName( L, l[8] );
    SetShortPresentation( L, l[4] );
    SetLibraryConditions( L, cond );
    
    return L;
end );

#############################################################################
##
#F  LIEPRING_NumberOfGroupsOfPrimePower( <p>, <k> )
##
##  Returns the number of groups of order p^k for k at most 7.
##
##  The counts for k at most 5 are classical, see
##      H. U. Besche, B. Eick and E. A. O'Brien, "A millennium project:
##      constructing small groups", Internat. J. Algebra Comput. 12 (2002),
##      623-644, doi:10.1142/S0218196702001115
##  the case k = 6 is due to
##      M. F. Newman, E. A. O'Brien and M. R. Vaughan-Lee, "Groups and
##      nilpotent Lie rings whose order is the sixth power of a prime",
##      J. Algebra 278 (2004), 383-401, doi:10.1016/j.jalgebra.2003.11.012
##  and the case k = 7 is due to
##      E. A. O'Brien and M. R. Vaughan-Lee, "The groups of order p^7 for odd
##      prime p", J. Algebra 292 (2005), 243-258,
##      doi:10.1016/j.jalgebra.2005.01.019
##
##  This avoids having to consult the SmallGrp package, and also covers cases
##  for which the small groups library has no data, such as p = 13, k = 7.
##
BindGlobal( "LIEPRING_NumberOfGroupsOfPrimePower", function( p, k )
    local g;
    g := d -> Gcd(p-1, d);
    if k = 1 then return 1;
    elif k = 2 then return 2;
    elif k = 3 then return 5;
    elif k = 4 then if p = 2 then return 14; else return 15; fi;
    elif k = 5 then
      if p = 2 then return 51; elif p = 3 then return 67; fi;
      return 61 + 2*p + 2*g(3) + g(4);
    elif k = 6 then
      if p = 2 then return 267; elif p = 3 then return 504; fi;
      return 3*p^2 + 39*p + 344 + 24*g(3) + 11*g(4) + 2*g(5);
    elif k = 7 then
      if p = 2 then return 2328; elif p = 3 then return 9310;
      elif p = 5 then return 34297; fi;
      return 3*p^5 + 12*p^4 + 44*p^3 + 170*p^2 + 707*p + 2455
             + (4*p^2 + 44*p + 291)*g(3)
             + (p^2 + 19*p + 135)*g(4)
             + (3*p + 31)*g(5) + 4*g(7) + 5*g(8) + g(9);
    fi;
    return fail;
end );

BindGlobal( "NumberOfLiePRings", function( arg )
    local dim, P;
    dim := arg[1];

    # if a prime is given
    if Length(arg) = 2 then 
        P := arg[2];
        if not IsPrimeInt(P) or P=2 then return fail; fi;
        if P < 7 then 
            Print("prime must be at least 7\n");
        elif dim > 7 then 
            Print("dimension must be at most 7\n");
        else
            # the number of Lie p-rings of this dimension equals the number
            # of groups of order p^dim
            return LIEPRING_NumberOfGroupsOfPrimePower( P, dim );
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
end );

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

    
