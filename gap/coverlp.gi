
LPRWeights := function(L)
    local d, s, t, b, w, i;

    # set up and check
    if not IsParentLiePRing(L) then return fail; fi;
    d := DimensionOfLiePRing(L);
    s := LiePLowerPCentralSeries(L);
    t := List(s, BasisOfLiePRing);
    b := t[1];

    # get weights
    w := List(b, x -> true);
    for i in [1..d] do
        w[i] := First([1..Length(t)], x -> not (b[i] in t[x]))-1;
        if w[i] = fail then return fail; fi;
    od;
    return w;
end;


LPRDefs := function(L)
    local d, b, S, v, i, u, c, j, k;

    # set up and check
    if not IsParentLiePRing(L) then return fail; fi;
    d := DimensionOfLiePRing(L);
    b := BasisOfLiePRing(L);
    S := SCTable(Zero(L)).tab;

    # get defs
    v := List(b, x -> true);
    for i in [1..d] do
        u := [];
        c := [];
        for j in [1..Length(S)] do
            k := Length(S[j]);
            if k > 0 and S[j][k-1] = i then 
                Add(u, j); 
                Add(c, S[j][k]);
            fi;
        od;
        if Length(c) > 0 then  
            j := First([1..Length(c)], x -> IsInt(c[x]));
            if IsInt(j) then v[i] := u[j]; fi;
        fi;
    od;

    return v;
end;

LiePCover := function(L)
    local d, p, S, v, k, l, T, c, i, j, a, b, I, Q, J, w, u, e;

    # check
    if not IsParentLiePRing(L) then return fail; fi;

    # set up
    d := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    S := SCTable(Zero(L));
    v := LPRDefs(L);
    k := Length(Filtered(v, x -> IsInt(x)));
    if not IsInt(p) then return fail; fi;

    # init cover
    l := d + d*(d-1)/2;
    T := rec( dim := d+l-k,
              prime := p,
              tab := [] );
    if IsBound(S.param) then T.param := S.param; fi;

    # add powers and conjugates
    c := d+1;
    for i in [1..d] do
        k := Sum([1..i-1]);
        for j in [1..i] do
            l := k+j;
            if IsBound(S.tab[l]) then 
                a := S.tab[l];
            else
                a := [];
            fi;
            if l in v then
                b := [];
            else
                b := [c,1]; c := c+1;
            fi;
            Add( T.tab, Concatenation( a, b ) );
        od;
    od;

    # evaluate pairs
    I := IdentityMat(T.dim);
    Q := MutableNullMat( T.dim, T.dim );
    for i in [1..T.dim] do
        for j in [1..T.dim] do
            Q[i][j] := LRMultiply( T, I[i], I[j] );
        od;
    od;

    # evaluate Jacobi
    J := [];
    for i in [1..T.dim] do
        for j in [1..T.dim] do
            for k in [1..T.dim] do
                a := LRMultiply( T, I[i], Q[j][k] );
                b := LRMultiply( T, I[j], Q[k][i] );
                c := LRMultiply( T, I[k], Q[i][j] );
                w := LRReduceExp( T, a+b+c );
                if w <> 0*w and not w in J and not -w in J then 
                    AddSet(J, w); 
                fi;
            od;
        od;
    od;

    # get a basis
    u := TriangulizedMat(BaseMat(J * One(GF(p)) ));
    I := IdentityMat( T.dim, GF(p) );
    v := BaseSteinitzVectors(I, u).factorspace;
    b := Concatenation(v,u)^-1;
    v := List(v, IntVecFFE);

    # set up new table
    Q := rec( dim := Length(v), prime := p, tab := [], param := []);
    for i in [1..Q.dim] do
        for j in [1..i-1] do
            e := LRMultiply( T, v[i], v[j] );
            e := e*b; e := e{[1..Q.dim]};
            Add(Q.tab, WordByExps(IntVecFFE(e)));
        od;
        e := LRReduceExp( T, p*v[i] );
        e := e*b; e := e{[1..Q.dim]};
        Add(Q.tab, WordByExps(IntVecFFE(e)));
    od; 

    return LiePRingBySCTableNC(Q);
end;

