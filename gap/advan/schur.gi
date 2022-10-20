         
#############################################################################
##
## 
##
BindGlobal( "FindNiceElm", function( units, elms )
    local d, m, u, p, a;

    p := IndeterminateByName("p");

    # catch the easy case
    if 1 in elms or p^0 in elms then return [p^0, true]; fi;
    if -1 in elms or -p^0 in elms then return [-p^0, true]; fi;

    # go 
    while true do

        d := List(elms, PDegree);
        m := Minimum(d);
        elms := Filtered(elms, x -> PDegree(x)=m);

        d := List(elms, FDegree);
        m := Minimum(d);
        elms := Filtered(elms, x -> FDegree(x)=m);

        for u in elms do
            a := LeadingUnit(u);
            if IsLiePUnit(units, a) then return [u, true]; fi;
        od;

        d := List(elms, FSummands);
        m := Minimum(d);
        elms := Filtered(elms, x -> FSummands(x)=m);
        
        #return [Random(elms), false];
        return [elms[1], false];
    od;

end );

#############################################################################
##
## Next pivot in matrix
##
BindGlobal( "FindPivot", function( R, m, i )
    local j, k, f, d, g, n, u;

    # reduce entries in m with zeros and filter non-zero elms
    f := [];
    for j in [i..Length(m)] do
       for k in [i..Length(m[j])] do
           m[j][k] := ReduceCoeffsByZeros(R!.zeros, m[j][k]);
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
            u := FindNiceElm( R!.units, g);
            if u[2] = true then 
                return MatPos(m,i,u[1]);
            else
                return u[1];
            fi;

        fi;
        d := d+1;
    od;
end );

#############################################################################
##
## SNF for the matrix determined by L - use units and zeros
##
BindGlobal( "GenericSNF", function(R, m)
    local n, d, p, u, i, e, j, k, w, a, A;

    # catch arguments
    p := IndeterminateByName("p");
    n := Length(m);
    d := Length(m[1]);
 
    # extend by 0-rows
    for i in [n+1..d] do Add(m, 0*[1..d]); od;

    for i in [1..d] do

        e := FindPivot(R, m, i );
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
end );
    
#############################################################################
##
## Some helpers
##
BindGlobal( "MakeIntPoly", function(R, vec)
    local p, i, u, e, f, uni;

    # set up
    p := IndeterminateByName("p");
    vec := p^0 * vec;
    uni := R!.units;

    # eliminate poly denominators
    for i in [1..Length(vec)] do
        u := DenominatorOfRationalFunction(vec[i]);
        if not IsLiePUnit(uni, u) then Error("no unit"); fi;
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
end );

#############################################################################
##
## Create matrix from structure constants
##
BindGlobal( "SetUpSchurMultSystem", function(L)
    local R, d, p, l, n, a, b, i, j, h, k, v, w, u, s, Pos;

    # catch arguments
    R := RingInvariants(L);
    d := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    l := BasisOfLiePRing(L);
    n := d * (d+1)/2;
    b := rec( ppps := [], vecs := [] );

    # a special case
    if d = 1 then return rec( norm := [], unit := []); fi;

    Pos := {i,j} -> i*(i-1)/2 + j;

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
            if w<>v then Add(b.ppps, MakeIntPoly(R,v-w)); fi;

            # [li, p lj]
            u := 0*[1..n];
            for k in [1..d] do 
                if a[j][j][k] <> 0 and k < i then 
                    u[Pos(i,k)] := a[j][j][k]; 
                elif a[j][j][k] <> 0 and i < k then 
                    u[Pos(k,i)] := -a[j][j][k]; 
                fi;
            od;
            if w<>u then Add(b.vecs, MakeIntPoly(R,w-u)); fi;
        od;

        # [p li,li] = p [li,li] = 0
        w := 0*[1..n];
        for k in [1..d] do
            if a[i][i][k] <> 0 then 
                w[Pos(k,i)] := a[i][i][k]; 
            fi;
        od;
        if w<>0*w then Add(b.vecs, MakeIntPoly(R,w)); fi;
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
                if w<>0*w then Add(b.vecs, MakeIntPoly(R,w)); fi;
            od;
        od;
    od;

    # add and return
    return Concatenation(b.vecs, b.ppps);
end );
    
#############################################################################
##
## Main function
##
BindGlobal( "LiePSchurMult", function(L)
    local S, b, c, R, T, i, V, W, h, w, pp;

    # set up
    S := [];
    b := SetUpSchurMultSystem(L);
    R := RingInvariants(L);
    w := IndeterminateByName("w");
    pp := ParametersOfLiePRing(L);

    # loop
    T := [[R.units, R.zeros]];
    i := 1;
    while i <= Length(T) do
#Print("start ",T[i],"\n");
        c := StructuralCopy(b);
        R := rec(units := T[i][1], zeros := T[i][2]);
        h := GenericSNF(R, c);
        if IsRecord(h) then 
            h.units := R.units;
            h.zeros := R.zeros;
            Add(S, h);
        else
            h := SQParts(pp, LeadingUnit(NumeratorOfRationalFunction(h)));
            V := CreateUnits( pp, T[i], h );
            W := CreateZeros( pp, T[i], h);
            if W <> fail and not W in T then Add(T, W); fi;
            if V <> fail and not V in T then Add(T, V); fi;
        fi;
        i := i+1;
    od;
          
    # that's it
    return S;
end );
