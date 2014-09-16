
RequirePackage("liepring");
LL := List([1..7], x -> LiePRingsByLibrary(x));

CheckFamily := function(L, P)
    local f, a, b;
    if P = 3 and PClassOfLiePRing(L) > 2 then return true; fi;
    f := NumberOfLiePRingsInFamily(L);
    a := EvaluatePorcPoly(f,P);
    b := LiePRingsInFamily(L,P);
    if b = fail then b := 0; else b := Length(b); fi;
    return a=b;
end;

CheckPrime := function(P)
    local i, j, L;
    for i in [1..7] do
        for j in [1..Length(LL[i])] do
#            if not (i=7 and j in [1798,1799,772,773,774,776,779,780,3285]) then
            if not (i=7 and j in [1798,1799]) then 
                L := LL[i][j];
                if CheckFamily(L, P) = false then 
                    Error("dim ",i," num ",j," prime ",P); 
                fi;
                Print("dim ",i," num ",j," prime ",P," done \n");
            fi;
        od;
    od;
end;

CheckSpecial := function(P)
    local L1, L2, f1, f2, a, b1, b2;
    if P = 3 then return true; fi;
    L1 := LL[7][1798];
    L2 := LL[7][1799];
    f1 := NumberOfLiePRingsInFamily(L1);
    f2 := NumberOfLiePRingsInFamily(L2);
    a := EvaluatePorcPoly(f1+f2,P);
    b1 := Length(LiePRingsInFamily(L1,P));
    b2 := Length(LiePRingsInFamily(L2,P));
    return a = b1+b2;
end;
  


