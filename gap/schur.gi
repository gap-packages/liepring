         
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
##
##
PValuePol := function(elm)
    local p, c, j; 
    p := IndeterminateByName("p");
    c := PolynomialCoefficientsOfPolynomial( elm*p^0, p );
    j := 1; while c[j] = 0*c[j] do j := j + 1; od; 
    return j-1;
end;

PValue := function( elm )
    local a, b;
    if elm = 0*elm then return fail; fi;
    if IsCRF(elm) then return 0; fi;
    a := PValuePol(NumeratorOfRationalFunction(elm));
    b := PValuePol(DenominatorOfRationalFunction(elm));
    if b <> 0 then Error("check this out"); fi;
    return a-b;
end;

PDegree := function( elm )
    local p, a, b, c, d;
    if elm = 0*elm then return fail; fi;
    if IsCRF(elm) then return 0; fi;
    p := IndeterminateByName("p");
    a := NumeratorOfRationalFunction(elm)*p^0;
    b := DenominatorOfRationalFunction(elm)*p^0;
    c := Length(PolynomialCoefficientsOfPolynomial( a, p ))-1;
    d := Length(PolynomialCoefficientsOfPolynomial( b, p ))-1;
    return c-d;
end;

FDegree := function(elm)
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
end;

#############################################################################
##
## elm, units in R
##
StripUnitPol := function( units, elm )
    local w, U, u, t;

    if not IsPolynomial(elm) then return fail; fi;

    # w is always a unit
    w := IndeterminateByName("w"); 
    U := Concatenation(units, [w]);
    U := Filtered(U, x -> not IsCRF(x));
 
    # loop
    for u in U do 
        t := true;
        while t = true do 
            t := false;
            if IsPolynomial(elm/u) then 
                elm := elm/u; 
                t := true;
            fi;
        od;
    od;

    return elm;
end;

IsMyUnit := function( units, elm )
    local a, b;
    if elm = 0*elm then return false; fi;
    if IsRat(elm) then return true; fi;
    if PDegree(elm)<>0 then return false; fi;
    a := NumeratorOfRationalFunction(elm);
    return IsCRF(StripUnitPol(units, a));
end;

LeadingUnit := function( elm )
    local p;
    if elm = 0*elm then return elm; fi;
    p := IndeterminateByName("p");
    return elm / p^PValue(elm);
end;

#############################################################################
##
## 
##
FindNiceElm := function( units, elms )
    local d, m, u, p, a;

    p := IndeterminateByName("p");
    while true do

        d := List(elms, PDegree);
        m := Minimum(d);
        elms := Filtered(elms, x -> PDegree(x)=m);

        d := List(elms, FDegree);
        m := Minimum(d);
        elms := Filtered(elms, x -> FDegree(x)=m);

        for u in elms do
            a := LeadingUnit(u);
            if IsMyUnit(units, a) then return [u, true]; fi;
        od;

        return [elms[1], false];
    od;

end;

#############################################################################
##
## 
##
RedPol := function( zeros, elm )
    if elm=0*elm or IsCRF(elm) then return elm; fi;
    if Length(zeros)=0 then return elm; fi;
    return PolynomialReducedRemainder(elm, zeros, MonomialLexOrdering());
end;

ReduceByZeros := function( zeros, elm )
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
end;

#############################################################################
##
## Next pivot in matrix
##
FindPivot := function( L, m, i )
    local j, k, f, d, g, n, u;

    # reduce entries in m with zeros and filter non-zero elms
    f := [];
    for j in [i..Length(m)] do
       for k in [i..Length(m[j])] do
           m[j][k] := ReduceByZeros(L!.zeros, m[j][k]);
           if m[j][k] <> 0*m[j][k] then Add(f, m[j][k]); fi;
       od;
    od;
    f := Unique(f);

    # check 
    if Length(f) = 0 then return 0; fi;

    # loop over degrees
    d := 0;
    while true do
        g := Filtered(f, x -> PValue(x)=d);
        if Length(g) > 0 then 

            # choose minimal element
            u := FindNiceElm( L!.units, g);
            if u[2] = true then 
                return MatPos(m,i,u[1]);
            else
                return u[1];
            fi;

        fi;
        d := d+1;
    od;
end;

#############################################################################
##
## SNF for the matrix determined by L - use units and zeros
##
GenericSNF := function(L, m)
    local n, d, p, u, i, e, j, k, w, a, A;

    # catch arguments
    p := IndeterminateByName("p");
    n := Length(m);
    d := Length(m[1]);
 
    # extend by 0-rows
    for i in [n+1..d] do Add(m, 0*[1..d]); od;

    for i in [1..d] do

        e := FindPivot(L, m, i );
#Print("  ",i," has pivot ",e,"\n");

        if e = 0*e then 
            m := List([1..d], x -> m[x][x]);
            m := Filtered(m, x -> x <> 0*x);
            m := Filtered(m, x -> x <> x^0);
            return rec( norm := m );
        elif not IsList(e) then 
            return e;
        else
            j := e[1]; k := e[2];

            # swap
            if j <> i then 
                w := m[i]; m[i] := m[j]; m[j] := w;
            fi;
            if k <> i then 
                m := TransposedMatMutable(m);
                w := m[i]; m[i] := m[k]; m[k] := w;
                m := TransposedMatMutable(m);
            fi;

            # norm
            a := LeadingUnit( m[i][i] );
            m[i] := m[i] * a^-1;

            # clear
            for j in [i+1..n] do
                m[j] := m[j] - (m[j][i]/m[i][i])*m[i];
            od;
            m := TransposedMatMutable(m);
            for k in [i+1..d] do
                m[k] := m[k] - (m[k][i]/m[i][i])*m[i];
            od;
            m := TransposedMatMutable(m);
        fi;
    od;
    m := List([1..Length(m[1])], x -> m[x][x]);
    m := Filtered(m, x -> x <> 0*x);
    m := Filtered(m, x -> x <> x^0);
    return rec(norm := m);
end;
    
#############################################################################
##
## Some helpers
##
Pos := function(i,j) return i*(i-1)/2 + j; end;

MakeIntPoly := function(L, vec)
    local p, i, u, e, f;

    # set up
    p := IndeterminateByName("p");
    vec := p^0 * vec;

    # eliminate poly denominators
    for i in [1..Length(vec)] do
        u := DenominatorOfRationalFunction(vec[i]);
        if not IsMyUnit(L!.units, u) then Error("no unit"); fi;
        vec := u*vec;
    od;

    # eliminate coeff denominators
    for i in [1..Length(vec)] do
        if vec[i] <> 0*vec[i] then 
            e := ExtRepPolynomialRatFun(vec[i]);
            e := Lcm(List([2,4..Length(e)], x -> DenominatorRat(x)));
            if e < 0 then e := (-1)*e; fi;
            f := Set(Factors(e));
            if ForAny(f, x -> x > 2) then Error("big units"); fi;
            vec := e*vec;
        fi;
    od;

    return vec;
end;

#############################################################################
##
## Create matrix from structure constants
##
SetUpSchurMultSystem := function(L)
    local d, p, l, n, a, b, i, j, h, k, v, w, u, s;

    # set up units
    if not IsParentLiePRing(L) then return fail; fi;
    if not IsBound(L!.units) and HasLibraryConditions(L) then 
        RingInvariants(L);
    elif not IsBound(L!.units) then 
        L!.units := []; L!.zeros := [];
    fi;

    # catch arguments
    d := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    l := BasisOfLiePRing(L);
    n := d * (d+1)/2;
    b := rec( ppps := [], vecs := [] );

    # a special case
    if d = 1 then return rec( norm := [], unit := []); fi;

    # precompute structure constants
    a := NullMat(d,d);
    for i in [1..d] do
        for j in [1..d] do
            if i = j then 
                a[i][i] := ExponentsLPR(L, p*l[i]);
            else
                a[i][j] := ExponentsLPR(L, l[i]*l[j]);
            fi;
        od;
    od;

    # evaluate p-th powers: p [li,lj] = [p li, lj] = [li, p lj]
    for i in [1..d] do
        for j in [1..i-1] do

            # p [li,lj]
            v := 0*[1..n];
            if i <> j then 
                v[Pos(i,j)] := p;
                for k in [1..d] do 
                    if a[i][j][k] <> 0 then 
                        v[Pos(k,k)] := a[i][j][k]; 
                    fi;
                od;
            fi;

            # [p li, lj]
            w := 0*[1..n];
            for k in [1..d] do 
                if a[i][i][k] <> 0 then 
                    w[Pos(k,j)] := a[i][i][k]; 
                fi;
            od;
            if w<>v then Add(b.ppps, MakeIntPoly(L,v-w)); fi;

            # [li, p lj]
            u := 0*[1..n];
            for k in [1..d] do 
                if a[j][j][k] <> 0 and k < i then 
                    u[Pos(i,k)] := a[j][j][k]; 
                elif a[j][j][k] <> 0 and i < k then 
                    u[Pos(k,i)] := -a[j][j][k]; 
                fi;
            od;
            if w<>u then Add(b.vecs, MakeIntPoly(L,w-u)); fi;
        od;

        # [p li,li] = p [li,li] = 0
        w := 0*[1..n];
        for k in [1..d] do
            if a[i][i][k] <> 0 then 
                w[Pos(k,i)] := a[i][i][k]; 
            fi;
        od;
        if w<>0*w then Add(b.vecs, MakeIntPoly(L,w)); fi;
    od;

    # evaluate jacobi
    for i in [1..d] do
        for j in [1..d] do
            for k in [1..d] do
                v := 0*[1..n];
                w := 0*[1..n];
                u := 0*[1..n];
                for h in [1..d] do
                    if j<>k and a[j][k][h] <> 0 and h < i then 
                        v[Pos(i,h)] := a[j][k][h];     
                    elif j<>k and a[j][k][h] <> 0 and i < h then 
                        v[Pos(h,i)] := -a[j][k][h];     
                    fi;
                od;
                for h in [1..d] do
                    if i<>j and a[i][j][h] <> 0 and h < k then 
                        u[Pos(k,h)] := a[i][j][h];     
                    elif i<>j and a[i][j][h] <> 0 and k < h then 
                        u[Pos(h,k)] := -a[i][j][h];     
                    fi;
                od;
                for h in [1..d] do
                    if k<>i and a[k][i][h] <> 0 and h < j then 
                        w[Pos(j,h)] := a[k][i][h];     
                    elif k<>i and a[k][i][h] <> 0 and j < h then 
                        w[Pos(h,j)] := -a[k][i][h];     
                    fi;
                od;
                w := v+w+u;
                if w<>0*w then Add(b.vecs, MakeIntPoly(L,w)); fi;
            od;
        od;
    od;

    # add and return
    return Concatenation(b.vecs, b.ppps);
end;

#############################################################################
##
## Expand input to a system of units and zeros
##
SetupUZSystem := function( U, Z )
    local done, zz, uu, x, y, z, w;

    done := false;
    while not done do
        done := true;

        if Length(Z) = 0 then return [U,Z]; fi;

        # reduce zeros with units
        zz := List(Z, x -> StripUnitPol(U,x));

        # expand zeros
        zz := ReducedGroebnerBasis(zz, MonomialLexOrdering());
        zz := List(zz, x -> Product(List(Collected(Factors(x)), y -> y[1])));
        zz := ReducedGroebnerBasis(zz, MonomialLexOrdering());
        if zz[1] = zz[1]^0 then return fail; fi;

        # reduce units with zeros
        uu := List(U, x -> RedPol(zz, x));
        if ForAny(uu, x -> x = 0*x) then return fail; fi;

        if zz <> Z then done := false; fi;
        Z := zz;

    od;
    
    # some special checks
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    w := IndeterminateByName("w");
 
    if w in Z then return fail; fi;
    if w^3-1/9*x^2 in Z then return fail; fi;
    if w-x^2 in Z or -w+x^2 in Z then return fail; fi;
    if w^-1-x^2 in Z or -w^-1+x^2 in Z then return fail; fi;
    if w-z^2 in Z or -w+z^2 in Z then return fail; fi;
    if w^3-z^2 in Z or -w^3+z^2 in Z then return fail; fi;
    if w^3-x^2 in Z or -w^3+x^2 in Z then return fail; fi;
    if w*y^2-1 in Z then return fail; fi;
    if w*x^2-1 in Z then return fail; fi;
    if -(x^2 - w*(w-1)^2) in Z then return fail; fi;
    if w-1 in Z or -(w-1) in Z then return fail; fi;
    if w * (w*y^2-1)^2 - (x*y)^2 in Z then return fail; fi;

    # now return result
    return [U,Z];

end;

#############################################################################
##
## Create units and zeros
##
CreateUnits := function( A, B )
    local b, u;

    # set up
    u := Unique( Concatenation( A[1], B ) );

    # create system
    return SetupUZSystem(u, A[2]);
end;

CreateZeros := function( A, B )
    local z;

    # set up
    z := Concatenation( A[2], [Product(B)] );

    # create system
    return SetupUZSystem( A[1], z );
end;

#############################################################################
##
## split up into factors
##
ReduceToSquareFree := function(h)
    local b, l, c, k;
    b := NumeratorOfRationalFunction(h);
    c := LeadingUnit(b);
    k := Collected(Factors(c));
    return List(k, x -> x[1]);
end;

#############################################################################
##
## Main function
##
LiePSchurMult := function(L)
    local b, c, R, U, Z, T, i, V, W, h, w;

    b := SetUpSchurMultSystem(L);
    R := [];
    U := L!.units;
    Z := L!.zeros;
    w := IndeterminateByName("w");

    T := [[U,Z]];
    i := 1;
    while i <= Length(T) do
#Print("start ",T[i],"\n");
        L!.units := T[i][1];
        L!.zeros := T[i][2];
        c := StructuralCopy(b);
        h := GenericSNF(L, c);
        if IsRecord(h) then 
            h.units := T[i][1];
            h.zeros := T[i][2];
            Add(R, h);
        else
            h := ReduceToSquareFree(h);
            V := CreateUnits( T[i], h );
            W := CreateZeros( T[i], h);
            if W <> fail then Add(T, W); fi;
            if V <> fail then Add(T, V); fi;
        fi;
        i := i+1;
    od;
          
    # reset
    L!.units := U;
    L!.zeros := Z;
            
    # that's it
    return R;
end;

#############################################################################
##
## Checking routines
##
LiePSchurMultByPrime := function(L, p)
    local F, G, c;
    F := LiePRingsInFamily(L,p);
    if F = fail then return F; fi;
    G := List(F, x -> PcGroupToPcpGroup(PGroupByLiePRing(x)));
    c := Collected(List(G, SchurMultPcpGroup));
    return List(c, x -> [List(x[1], y -> Length(Factors(y))),x[2]]);
end;

CheckLiePSM := function(L)
    local fix, gen, p;

    p := 7;
    fix := fail;
    while fix = fail do
        fix := LiePSchurMultByPrime(L,p);
        p := NextPrimeInt(p);
    od;
    fix := Set(List(fix, x -> x[1]));
    Print(fix,"\n");

    gen := LiePSchurMult(L);
    gen := List(gen, x -> List(x.norm, y -> Length(Factors(y))));
    gen := Set(gen);
    Print(gen,"\n");

    return fix = gen;
end;

CheckOrderPSM := function(d, k)
    local LL, res, i, t;
    LL := LiePRingsByLibrary(d);
    res := [];
    for i in [1..Length(LL)] do
        if Length(ParametersOfLiePRing(LL[i]))=k then 
            t := CheckLiePSM(LL[i]);
            if t = false then Add(res, [i,t]); fi;
        fi;
        Print(i, "  done \n");
    od;
    return res;
end;

