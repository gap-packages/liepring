         
LPRWeights := function(L)
    local d, s, t, b, w, i;

    # set up and check
    if not IsParentLiePRing(L) then return fail; fi;
    d := DimensionOfLiePRing(L);
    s := LiePLowerPCentralSeries(L);
    b := BasisOfLiePRing(L);

    # get weights
    w := List(b, x -> true);
    for i in [1..d] do
        w[i] := First([1..Length(s)], x -> not (b[i] in s[x]))-1;
    od;

    # return
    return w;
end;

LPRDefs := function(L)
    local d, b, S, v, c, i, j, k, e;

    # check
    if not IsParentLiePRing(L) then return fail; fi;
    if IsBound(L!.defs) then return L!.defs; fi;

    # set up
    d := DimensionOfLiePRing(L);
    b := BasisOfLiePRing(L);
    S := SCTable(Zero(L)).tab;
    v := List(b, x -> true);

    # get defs
    c := 0;
    for i in [1..d] do
        for j in [1..i] do
            c := c + 1;
            if IsBound(S[c]) and Length(S[c]) > 0 then 
                k := S[c][Length(S[c])-1];
                e := S[c][Length(S[c])];
                if k > i and v[k] = true and (e in [1,-1] or IsRootPower(e)) 
                    then 
                    v[k] := [i,j]; 
                elif k > i and e in [1,-1] then 
                    v[k] := [i,j];
                fi;
            fi;
        od;
    od;
  
    # add them and return
    L!.defs := v;
    return v;
end;

InsertVec := function( T, J, d, vec )
    local wec;
    wec := LRReduceExp( T, vec );
    if ForAny(wec{[1..d]}, x -> x <> 0*x) then return fail; fi;
    if wec <> 0*wec and not wec in J and not -wec in J then Add(J, wec); fi;
    return J;
end;

LiePCover := function(L)
    local d, p, S, v, T, c, r, i, j, k, I, Q, P, J, a, b, w, u, e, s, R;

    # check
    if not IsParentLiePRing(L) then return fail; fi;

    # set up
    d := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    S := SCTable(Zero(L));
    v := LPRDefs(L);
    w := LPRWeights(L);
    for i in [1..d] do
        for j in [1..i-1] do
            if w[i]+w[j] > w[d]+1 then Add(v, [i,j]); fi;
        od;
    od;

    # init cover
    T := rec( prime := p, tab := [] );

    # add tails
    c := 0; r := d;
    for i in [1..d] do
        for j in [1..i] do
            c := c+1;
            if IsBound(S.tab[c]) then 
                T.tab[c] := StructuralCopy(S.tab[c]);
            else
                T.tab[c] := [];
            fi;
            if not [i,j] in v then
                r := r+1;
                Append( T.tab[c], [r,1]);
            fi;
        od;
    od;
    T.dim := r;

    # evaluate powers and pairs
    I := IdentityMat(T.dim);
    Q := MutableNullMat( T.dim, T.dim );
    P := [];
    for i in [1..T.dim] do
        P[i] := LRReduceExp(T, p*I[i]);
        for j in [1..T.dim] do
            Q[i][j] := LRMultiply( T, I[i], I[j] );
        od;
    od;

    # evaluate relations
    J := [];

    # test biadditivity
    for i in [1..T.dim] do
        c := LRMultiply(T, P[i], I[i]);
        J := InsertVec(T, J, S.dim, c);
        if J = fail then Error("bi-add case 1"); fi;

        for j in [1..i-1] do
            a := LRReduceExp(T, p*Q[i][j]);
            b := LRMultiply(T, P[i], I[j]);
            c := LRMultiply(T, I[i], P[j]);
            J := InsertVec(T, J, S.dim, a-c);
            if J = fail then Error("bi-add case 2"); fi;
            J := InsertVec(T, J, S.dim, a-b);
            if J = fail then Error("bi-add case 3"); fi;
        od;
    od;

    # test Jacobi
    for i in [1..T.dim] do
        for j in [1..T.dim] do
            for k in [1..T.dim] do
                a := LRMultiply(T, I[i], Q[j][k]);
                b := LRMultiply(T, I[j], Q[k][i]);
                c := LRMultiply(T, I[k], Q[i][j]);
                J := InsertVec(T, J, S.dim, a+b+c);
                if J = fail then Error("Jacobi"); fi;
            od;
        od;
    od;

    # check
    if IsBound(L!.inv) then 
        u := MyBaseMat(J, L!.inv);
    else
        u := MyBaseMat(J, []);
    fi;
    if u = fail then return fail; fi;

    # get quotient
    if IsBound(S.param) then T.param := S.param; fi;
    R := LiePQuotientByTable(T, u);

    # add info
    if IsBound(L!.inv) then R!.inv := L!.inv; fi;
    b := BasisOfLiePRing(R);
    R!.mult := LiePSubringByBasis( R, b{[S.dim+1..Length(b)]});
    s := LiePLowerPCentralSeries(R);
    R!.nucl := s[PClassOfLiePRing(L)+1];
    return R;
end;

StepSize := function(L)
    local C;
    C := LiePCover(L);
    return DimensionOfLiePRing(C!.nucl);
end;

IsTerminal := function(L)
    return StepSize(L)=0;
end;


