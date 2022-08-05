BindGlobal( "InsertVec", function( T, J, d, vec )
    local wec;
    wec := LRReduceExp( T, vec );
    if ForAny(wec{[1..d]}, x -> x <> 0*x) then return fail; fi;
    if wec <> 0*wec and not wec in J and not -wec in J then Add(J, wec); fi;
    return J;
end );

BindGlobal( "LiePCover", function(L)
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
    T.ring := S.ring;

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
    Q := NullMat( T.dim, T.dim );
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
    u := MyBaseMat@(J, T.ring.units);
    if not IsList(u) then return u; fi;

    # get quotient
    if IsBound(S.param) then T.param := S.param; fi;
    R := LiePQuotientByTable(T, u);

    # add info
    b := BasisOfLiePRing(R);
    R!.mult := LiePSubringByBasis( R, b{[S.dim+1..Length(b)]});
    s := LiePLowerPCentralSeries(R);
    if s <> fail then R!.nucl := s[PClassOfLiePRing(L)+1]; fi;
    return R;
end );

BindGlobal( "IsTerminalLiePRing", function(L)
    local C;
    C := LiePCover(L);
    if IsLiePRing(C) and IsBound(C!.nucl) then 
        return DimensionOfLiePRing(C!.nucl)=0;
    else
        return fail;
    fi;
end );

BindGlobal( "StepSizeLiePRing", function(L)
    local C;
    C := LiePCover(L);
    if IsLiePRing(C) and IsBound(C!.nucl) then 
        return DimensionOfLiePRing(C!.nucl);
    else
        return fail;
    fi;
end );

BindGlobal( "AutoOnMult", function( C, mat )
    local dim, sml, new;
    dim := DimensionOfLiePRing(C);
    sml := dim - DimensionOfLiePRing(C!.mult);
    new := ExtendAuto(C,mat);
    return new{[sml+1..dim]}{[sml+1..dim]};
end );



