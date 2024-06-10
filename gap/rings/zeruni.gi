#############################################################################
##
##  R = Q[w,x,y,z,t,...]
##  P = R[p]
##
##  g in P hat die Form g = g_0 + g_1 p + ... + g_t p^t
##    Einheit: g ist Einheit, falls g_0 Einheit ist
##    PDegree: t = PDegree(g)
##    PValue : falls g = p^l h mit h_0 <> 0, dann l = PValue(g)
##
##  U = {g in R | g Einheit}
##  W = Localisation von R bei U = {g/h | g in R, h in U}
##
##  g in W hat die Form g = p^l a/b mit a_0 <> 0, b in U
##    Einheit: g ist Einheit, falls l = 0 und a Einheit
##    PDegree: PDegree(g) = l + PDegree(a) - PDegree(b)
##    PValue : PValue(g) = l
##

#############################################################################
##
## Zeros 
##

BindGlobal( "ReduceCoeffsByZeros", function( zeros, elm )
    local p, a, c, b;

    if elm = 0*elm then return elm; fi;

    p := IndeterminateByName("p");

    a := NumeratorOfRationalFunction(elm);
    c := PolynomialCoefficientsOfPolynomial(a, p);
    c := List(c, x -> RedPol(zeros,x));
    a := Sum(List([1..Length(c)], x -> c[x]*p^(x-1)));

    b := DenominatorOfRationalFunction(elm);
    c := PolynomialCoefficientsOfPolynomial(b, p);
    c := List(c, x -> RedPol(zeros,x));
    b := Sum(List([1..Length(c)], x -> c[x]*p^(x-1)));

    return a/b;
end );

#############################################################################
##
## Units 
##

BindGlobal( "ReduceByUnits", function( units, elm )
    local w, U, u, t, f;

    # the denominator is always a unit
    elm := NumeratorOfRationalFunction(elm);

    # set up
    w := IndeterminateByName("w");
    U := Concatenation(units, [w]);
    Append(U, List(VarsOfPoly(elm), x -> x^2-w));

    # loop
    for u in U do
        t := true;
        while t = true do
            t := false;
            f := elm/u;
            if IsPolynomial(f) then elm := f; t := true; fi;
        od;
    od;

    return elm;
end );

BindGlobal( "IsLiePUnit", function( uni, elm )
    local a, b;
    if elm = 0*elm then return false; fi;
    if IsRat(elm) then return true; fi;
    if PDegree(elm)<>0 then return false; fi;
    a := NumeratorOfRationalFunction(elm);
    return IsCRF(ReduceByUnits(uni, a));
end );

BindGlobal( "LeadingUnit", function( elm )
    local p;
    if elm = 0*elm then return elm; fi;
    p := IndeterminateByName("p");
    return elm / p^PValue(elm);
end );

#############################################################################
##
## Expand input to a system of units and zeros
## Returns either [U,Z] or fail if an inconsistency is detected
##
BindGlobal( "SetupUZSystem", function( pp, U, Z )
    local done, zz, uu, x, y, z, w;

    # set up
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    w := IndeterminateByName("w");

    # Groebner basis for Z 
    if Length(Z)>0 then 
        Z := SquareFreeGB(Z); 
        if Length(Z) = 1 and (IsInt(Z[1]) or IsCRF(Z[1])) then return fail; fi;
    
        if w-1 in Z or -(w-1) in Z then return fail; fi;

        if w-x^2 in Z or -w+x^2 in Z then return fail; fi;
        if w-y^2 in Z or -w+y^2 in Z then return fail; fi;
        if w-z^2 in Z or -w+z^2 in Z then return fail; fi;

        if w^3-x^2 in Z or -w^3+x^2 in Z then return fail; fi;
        if w^3-y^2 in Z or -w^3+y^2 in Z then return fail; fi;
        if w^3-z^2 in Z or -w^3+z^2 in Z then return fail; fi;

        if w*x^2-1 in Z then return fail; fi;
        if w*y^2-1 in Z then return fail; fi;
        if w*z^2-1 in Z then return fail; fi;

        if w^3-1/9*x^2 in Z then return fail; fi;
        if 9*w^3-x^2 in Z then return fail; fi;
        if w*z^2-y^2 in Z then return fail; fi;
    fi;

    # reduction for U
    U := List(U, x -> RedPol(Z, x));
    if ForAny(U, x -> x = 0*x) then return fail; fi;

    U := Unique(Concatenation(List(U, x -> SQParts(pp,x))));
    U := Filtered(U, x -> not IsCRF(x));
    U := List(U, x -> x / LeadingCoefficient(x));
    U := Filtered(U, x -> x <> w);

    # now return result
    return [U,Z];
end );

#############################################################################
##
## Create units and zeros
##
BindGlobal( "CreateUnits", function( pp, A, B )
    return SetupUZSystem(pp, Unique(Concatenation(A[1],B)), A[2]);
end );

BindGlobal( "CreateZeros", function( pp, A, B )
    return SetupUZSystem( pp, A[1], Concatenation(A[2],[Product(B)]) );
end );

