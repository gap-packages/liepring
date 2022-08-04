         
BindGlobal( "LPRWeights", function(L)
    local d, s, t, b, w, i;

    # set up and check
    if not IsParentLiePRing(L) then return fail; fi;
    d := DimensionOfLiePRing(L);
    s := LiePLowerPCentralSeries(L);
    if s = fail then return fail; fi;
    b := BasisOfLiePRing(L);

    # get weights
    w := List(b, x -> true);
    for i in [1..d] do
        w[i] := First([1..Length(s)], x -> not (b[i] in s[x]))-1;
    od;

    # return
    return w;
end );

BindGlobal( "LPRDefs", function(L)
    local d, b, T, S, U, v, c, i, j, k, e;

    # check
    if not IsParentLiePRing(L) then return fail; fi;
    if IsBound(L!.defs) then return L!.defs; fi;

    # set up
    d := DimensionOfLiePRing(L);
    b := BasisOfLiePRing(L);
    T := SCTable(Zero(L));
    U := T.ring.units;
    S := T.tab;
    v := List(b, x -> true);

    # get defs
    c := 0;
    for i in [1..d] do
        for j in [1..i] do
            c := c + 1;
            if IsBound(S[c]) and Length(S[c]) > 0 then 
                k := S[c][Length(S[c])-1];
                e := S[c][Length(S[c])];
                if k > i and v[k] = true and IsLiePUnit(U, e) then
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
end );

