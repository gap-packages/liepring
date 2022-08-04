
##
## determine variables in polynomial or list of polynomials
##
BindGlobal( "VarsOfPoly", function( poly )
    local e, v, i, j;
    e := ExtRepPolynomialRatFun(poly);
    v := [];
    for i in [1,3..Length(e)-1] do
        for j in [1,3..Length(e[i])-1] do
            Add(v, e[i][j]);
        od;
    od;
    v := Set(v);
    v := List(v, x -> Indeterminate(Rationals, x));
    v := Filtered(v, x -> x <> IndeterminateByName("w"));
    return v;
end );

BindGlobal( "VarsOfSCTab", function( list )
    local l;
    l := Filtered(Flat(list), x -> IsPolynomial(x) );
    l := Set(Flat(List(l, VarsOfPoly)));
    l := List(l, x -> IndeterminateNumberOfUnivariateRationalFunction(x));
    return List(Set(l), x -> Indeterminate(Rationals, x));
end );

##
## get coeffs
##
BindGlobal( "CoefficientsOfPolynomial", function( f )
    local e;
    e := ExtRepPolynomialRatFun(f);
    return e{[2,4..Length(e)]};
end );

##
## split rat fun a/b into a and b so that a and b have int coeffs
##
BindGlobal( "SplitRatFun", function( rf )
    local a, b, u, v;
    a := NumeratorOfRationalFunction(rf);
    b := DenominatorOfRationalFunction(rf);
    u := Lcm(List(CoefficientsOfPolynomial(a), DenominatorRat));
    v := Lcm(List(CoefficientsOfPolynomial(b), DenominatorRat));
    return Lcm(u,v)*[a,b];
end );

##
## reduce mod zeros
##
BindGlobal( "RedPol", function( zeros, elm )
    local u, v, a, b;

    # set up
    if Length(zeros) = 0 then return elm; fi;
    if elm = 0*elm or IsRat(elm) then return elm; fi;
    u := SplitRatFun(elm); v := u[2]; u := u[1];

    a := PolynomialDivisionAlgorithm( u, zeros, ORDER );
    b := PolynomialDivisionAlgorithm( v, zeros, ORDER );
    if a[1] = 0*a[1] then return 0; fi;
    if b[1] = 0*b[1] then Error("division by zero"); fi;

    return a[1]/b[1];
end );


##
## degree of multivariate poly
##
BindGlobal( "DegreeOfPoly", function( poly )
    local e, d, i, s, j, w;
    e := ExtRepPolynomialRatFun(poly);
    w := ExtRepPolynomialRatFun(IndeterminateByName("w"))[1][1];
    d := 0;
    for i in [1,3..Length(e)-1] do
        s := 0;
        for j in [2,4..Length(e[i])] do
            if e[i][j-1] <> w then s := s + e[i][j]; fi;
        od;
        d := Maximum(d,s);
    od;
    return d;
end );

##
## factors of poly
##
BindGlobal( "SQPart", function(pp, h)
    local f;
    f := Collected(MyFactors(pp, h));
    return Product(List(f, x -> x[1]));
end );

BindGlobal( "SQParts", function(pp, h)
    local f;
    f := Collected(MyFactors(pp, h));
    return List(f, x -> x[1]);
end );

##
## constant poly versus int
##
BindGlobal( "IsCRF", function ( elm )
    if IsRat( elm )  then return true; fi;
    return IsConstantRationalFunction( elm );
end );

DeclareGlobalFunction( "MakeInt" );
InstallGlobalFunction( "MakeInt", function(M)
    local i, j;
    if IsList(M) then
        return List(M, x -> MakeInt(x));
    elif IsInt(M) then
        return M;
    elif M = 0*M then
        return 0;
    elif IsUnivariatePolynomial(M) and Degree(M) = 0 then
        return CoefficientsOfUnivariatePolynomial(M)[1];
    else
        return M;
    fi;
end );

##
## groebner square free
##
BindGlobal( "SquareFreeGB", function( polys )
    local f, S, w, i;
    f := ShallowCopy(polys);
    w := IndeterminateByName("w");
    repeat
        S := ShallowCopy(f);
        f := List(f, pol -> Product(List(Collected(Factors(pol)),x ->x[1])));
        for i in [1..Length(f)] do
            while IsPolynomial(f[i]/w) do f[i] := f[i]/w; od;
        od;
        f := CallGroebner(f);
    until S = f; 
    return f;
end );

#############################################################################
##
## Valuation and Degree
##
BindGlobal( "PValuePol", function(elm)
    local p, c, j;
    p := IndeterminateByName("p");
    c := PolynomialCoefficientsOfPolynomial( elm*p^0, p );
    j := 1; while c[j] = 0*c[j] do j := j + 1; od;
    return j-1;
end );

BindGlobal( "PValue", function( elm )
    local a, b;
    if elm = 0*elm then return fail; fi;
    if IsCRF(elm) then return 0; fi;
    a := PValuePol(NumeratorOfRationalFunction(elm));
    b := PValuePol(DenominatorOfRationalFunction(elm));
    if b <> 0 then Error("check this out"); fi;
    return a-b;
end );

BindGlobal( "PDegree", function( elm )
    local p, a, b, c, d;
    if elm = 0*elm then return fail; fi;
    if IsCRF(elm) then return 0; fi;
    p := IndeterminateByName("p");
    a := NumeratorOfRationalFunction(elm)*p^0;
    b := DenominatorOfRationalFunction(elm)*p^0;
    c := Length(PolynomialCoefficientsOfPolynomial( a, p ))-1;
    d := Length(PolynomialCoefficientsOfPolynomial( b, p ))-1;
    return c-d;
end );

BindGlobal( "FDegree", function(elm)
    local p, a, b, e, m;
    if elm=0*elm then return -infinity; fi;
    if IsCRF(elm) then return 0; fi;
    p := IndeterminateByName("p");
    a := NumeratorOfRationalFunction(elm)*p^0;
    b := DenominatorOfRationalFunction(elm)*p^0;
    e := ExtRepPolynomialRatFun(a);
    e := e{[1,3..Length(e)-1]};
    m := List(e, x -> Sum(x{[2,4..Length(x)]}));
    return Maximum(m);
end );

BindGlobal( "FSummands", function(elm)
    local a, b, p, e;
    if elm=0*elm then return 0; fi;
    if IsCRF(elm) then return 1; fi;
    p := IndeterminateByName("p");
    a := NumeratorOfRationalFunction(elm)*p^0;
    b := DenominatorOfRationalFunction(elm)*p^0;
    e := ExtRepPolynomialRatFun(a);
    return Length(e)/2;
end );


