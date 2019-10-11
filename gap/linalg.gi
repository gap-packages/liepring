
MakeInt := function(M)
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
end;

DepthVector := function(vec)
    local i;
    for i in [1..Length(vec)] do
        if vec[i] <> 0*vec[i] then return i; fi;
    od;
    return Length(vec)+1;
end;

Swapped := function( J, k, j, d )
    local t;
    t := J[k];
    J[k] := J[j];
    J[j] := t;
    J[k] := J[k][d]^-1 * J[k];
    return J;
end;

IsRootPower := function(v)
    local a, b, c;

    # trivial cases
    if v = 0*v then return false; fi;
    if not IsLaurentPolynomial(v) then return false; fi;

    # check if v is a polynomial over w
    a := IndeterminateNumberOfLaurentPolynomial(IndeterminateByName("w"));
    b := IndeterminateNumberOfLaurentPolynomial(v);
    if a <> b then return false; fi;
 
    # final check
    c := CoefficientsOfLaurentPolynomial(v);
    if Length(c[1]) <> 1 then return false; fi;
    return c[1][1] in [1,-1,2,-2];
end;

IsLiePUnit := function(L, v)
    local u;

    # the trivial case
    if v = 0*v then return false; fi;

    # make int if possible
    v := MakeInt(v);

    # the integer case
    if IsRat(v) then 
        v := NumeratorRat(v);
        if (v in [1,-1,2,-2] or v in L!.inv) then 
            return true;
        else
            return fail;
        fi;
    fi;

    # the polynomial case
    if IsRationalFunction(v) then 
        v := NumeratorOfRationalFunction(v);
        if v in L!.inv then return true; fi;
        u := v/IndeterminateByName("w");
        if IsPolynomial(u) then 
            return IsUnit(L, u);
        else
            return fail;
        fi;
    fi;

end;

Pivot := function(J, d, k, inv)
    local v, m, j;

    # cut out relevant part
    v := J{[k..Length(J)]}[d];

    # check for ints
    m := Filtered(v, x -> (x <> 0 and IsInt(x)));
    if Length(m) > 0 then 
        m := Minimum(m);
        j := Position(v,m);
        if not m in [1,-1,2,-2] then Print("inverting ",v[j],"\n"); fi;
        return Swapped(J, k, j+k-1, d); 
    fi;

    # then check for w
    j := First([1..Length(v)], x -> IsRootPower(v[x]));
    if IsInt(j) then return Swapped(J, k, j+k-1, d); fi;

    # now consider invertibles
    j := First([1..Length(v)], x -> (v[x] in inv) or (v[x] in -inv));
    if IsInt(j) then return Swapped(J, k, j+k-1, d); fi;

    # finally print a statement
    Print("need invertible in ", v,"\n");
    return fail;
end;

MyBaseMat := function(J, inv)
    local k, d, i, l, m;

    if Length(J) = 0 then return J; fi;

    J := StructuralCopy(J);
    l := Length(J);
    m := Length(J[1]);
    k := 0;
    repeat

        # get positions
        k := k+1;
        if k > Length(J) then return J; fi;
        d := Minimum(List(J{[k..l]}, DepthVector));
        if d = m+1 then return J{[1..k-1]}; fi;

        # move pivot
        J := Pivot(J, d, k, inv);
        if J = fail then return fail; fi;

        # clear out
        for i in [1..Length(J)] do
            if i <> k and J[i][d] <> 0*J[i][d] then
                J[i] := J[i] - J[i][d]*J[k];
            fi;
        od;
    until false;
end;

FactorSpace := function( d, sub )
    return IdentityMat(d){Difference([1..d], List(sub, DepthVector))};
end;


