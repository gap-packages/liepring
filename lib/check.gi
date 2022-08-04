
LoadPackage("anupq");

BindGlobal( "SCodeByGroup", function( G )
    local F, S, K; 
    F := FpGroupPcGroup(G);
    S := StandardPresentation( F : Prime := PrimePGroup(G) );
    K := PcGroupFpGroup(S);
    return CodePcGroup(K);
end );

BindGlobal( "CheckNumberLiePRings", function(dim, P)
    local L;
    L := LiePRingsByLibrary(dim,P);
    return Length(L) = NumberSmallGroups(P^dim);
end );

BindGlobal( "FingerPrintLCS", function( G )
    local s, t;
    s := LowerCentralSeries(G);
    t := List([1..Length(s)-1], x -> AbelianInvariants(s[x]/s[x+1]));
    return t;
end );

BindGlobal( "FingerPrintJen", function( G )
    local s, t;
    s := JenningsSeries(G);
    t := List([1..Length(s)-1], x -> AbelianInvariants(s[x]/s[x+1]));
    return t;
end );

BindGlobal( "FingerPrintPCS", function(G)
    local s, t;
    s := PCentralSeries(G);
    t := List([1..Length(s)-1], x -> Size(s[x]/s[x+1]));
    return t;
end );

BindGlobal( "FingerPrintAut", function(G)
    return AutomorphismGroupPGroup(G).size;
end );

BindGlobal( "FingerPrintCQ", function(G)
    local c, Q;
    c := Center(G);
    Q := G/c;
    if Length(Factors(Size(Q))) <= 4 then 
        return IdGroup(Q); 
    else
        return Size(Q);
    fi;
end );

BindGlobal( "FingerPrintPRump", function(G)
    local U;
    U := PRump(G, PrimePGroup(G));
    if Length(Factors(Size(U))) <= 4 then 
        return IdGroup(U); 
    else
        return Size(U);
    fi;
end );
      
BindGlobal( "FingerPrintDerived", function(G)
    local U;
    U := DerivedSubgroup(G);
    if Length(Factors(Size(U))) <= 4 then 
        return IdGroup(U); 
    else
        return Size(U);
    fi;
end );

BindGlobal( "FingerPrintClasses", function(G)
    local H, cl, rs, i, a, b, c, g, e;
    H := PcGroupCode(CodePcGroup(G),Size(G));
    cl := ConjugacyClasses(H);
    rs := [];
    for i in [1..Length(cl)] do
        g := Representative(cl[i]);
        a := Size(cl[i]);
        b := Order(g);
        rs[i] := [a,b];
    od;
    return rs;
end );

BindGlobal( "FingerPrintMaxSub", function(G)
    local n, t, i;
    n := Length(Factors(Size(G)));
    t := [G];
    for i in [n, n-1 .. 5] do
        t := Set(Flat(List(t, MaximalSubgroups)));
    od;
    return Collected(List(t, IdGroup));
end );

BindGlobal( "RefineBin", function(list, func)
    local val, spl;
    val := List(list, x -> func(x));
    spl := List(Set(val), x -> Filtered([1..Length(val)], y->val[y]=x));
    return List(spl, x -> list{x});
end );
   
BindGlobal( "BinsByGroups", function(grps)
    local bins, abel, invs, s; 

    # abelian groups
    abel := Filtered(grps, x -> IsAbelian(x));
    invs := List(abel, AbelianInvariants);
    if not Length(Set(invs)) = Length(invs) then 
        Error("abelian groups do not split");
    fi;
    Print("bins for ",Length(grps)-Length(abel)," groups \n");

    # other
    bins := [Filtered(grps, x -> not IsAbelian(x))];

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintLCS)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintPCS)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintDerived)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintJen)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintCQ)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintPRump)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintAut)));
    bins := Filtered(bins, x -> Length(x)>1);
    s := List(bins, Length);
    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
    Print("  ",SortedList(s),"\n");

#    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintMaxSub)));
#    bins := Filtered(bins, x -> Length(x)>1);
#    s := List(bins, Length);
#    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
#    Print("  ",SortedList(s),"\n");
#
#    bins := Concatenation(List(bins, x -> RefineBin(x, FingerPrintClasses)));
#    bins := Filtered(bins, x -> Length(x)>1);
#    s := List(bins, Length);
#    Print("  got ",Length(bins)," bins with ",Sum(s)," groups \n");
#    Print("  ",SortedList(s),"\n");

    return bins;
end );

BindGlobal( "CheckLiePRings", function( P, dim, gen, pcl )
    local L, G, A, S, i, T;

    # get groups
    L := LiePRingsByLibrary(dim, P);
    S := [];
    for i in [1..Length(L)] do
        G := PGroupByLiePRing(L[i]);
        Print("checking ",i);
        if not IsAbelian(G) then 
            if Index(G, PRump(G, P)) = P^gen then 
                if (not IsInt(pcl)) or 
                   (IsInt(pcl) and Length(PCentralSeries(G)) = pcl+1) 
                then 
                    S[i] := SCodeByGroup(G);
                    Print(" -- is relevant");
                fi;
            fi;
        fi;
        Print(" \n");
    od;

    T := Filtered(S, x -> x <> true);
    return Length(T) = Length(Set(T));
end );

BindGlobal( "CheckLiePRingsDim7ByFilePlus", function( nr, P )
    local L, spe, v, r, i, G, S;

    # read the desired lie p-rings
    LIE_DATA[7] := [];
    ReadPackage("liepring", Concatenation("lib/dim7/",LIE_TABLE[nr][1]));

    # convert to gap lie p-rings
    Print("  get lie rings \n");
    L := List(LIE_DATA[7], x -> LiePRingByData(7, x));

    # add primes
    Print("  specialize \n");
    spe := List( L, X -> LiePRingsInFamily( X, P ) );
    spe := Flat(Filtered(spe, l -> l <> fail ));

    # check number
    Print("  check number \n");
    v := EvaluatePorcPoly( LIE_TABLE[nr][2], P );
    if Length(spe) <> v then
        Print("wrong number: ",Length(spe)-v," off \n");
        return fail;
    fi;

    Print("  check isoms \n");
    r := [];
    for i in [1..Length(spe)] do
        G := PGroupByLiePRing(spe[i]);
        S := SCodeByGroup(G);
        Add(r, S);
        Print("    ",i," of ",Length(spe)," done \n");
    od;
    return Length(r) = Length(Set(r));
end );

BindGlobal( "CheckLiePRingsDim7Gen2", function( P )
    local t, i;
    t := [];
    for i in [2..59] do
        Print("starting file ",i," \n");
        t[i] := CheckLiePRingsDim7ByFile(i,P);
        Print("-- got ",t[i],"\n");
    od;
    return t;
end );

