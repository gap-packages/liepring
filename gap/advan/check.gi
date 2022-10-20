
#############################################################################
##
## Checking routines
##
BindGlobal( "LiePSchurMultByPrime", function(L, p)
    local F, G, c;
    F := LiePRingsInFamily(L,p);
    if F = fail then return F; fi;
    G := List(F, x -> PcGroupToPcpGroup(PGroupByLiePRing(x)));
    c := Collected(List(G, SchurMultPcpGroup));
    return List(c, x -> [List(x[1], y -> Length(Factors(y))),x[2]]);
end );

BindGlobal( "CheckLiePSM", function(L)
    local fix, gen, p;

    p := 7;
    fix := fail;
    while fix = fail do
        fix := LiePSchurMultByPrime(L,p);
        p := NextPrimeInt(p);
    od;
    fix := Set(List(fix, x -> x[1]));
    #Print(fix,"\n");

    gen := LiePSchurMult(L);
    gen := List(gen, x -> List(x.norm, y -> Length(Factors(y))));
    gen := Set(gen);
    #Print(gen,"\n");

    return fix = gen;
end );

BindGlobal( "CheckOrderPSM", function(d, k)
    local LL, res, i, t;
    LL := LiePRingsByLibrary(d);
    res := [];
    for i in [1..Length(LL)] do
        if Length(ParametersOfLiePRing(LL[i]))=k then 
            t := CheckLiePSM(LL[i]);
            if t = false then Add(res, [i,t]); fi;
        fi;
        #Print(i, "  done \n");
    od;
    return res;
end );
