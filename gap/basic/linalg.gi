# DepthVector differs from PositionNonZero in that it also works on
# inhomogeneous lists
DepthVector := function( vec )
    local k;
    k := PositionProperty(vec, k -> not IsZero(k));
    if k = fail then return Length(vec) + 1; fi;
    return k;
end;

BindGlobal( "Swapped", function( J, k, j, d )
    local t;
    t := J[k];
    J[k] := J[j];
    J[j] := t;
    J[k] := J[k][d]^-1 * J[k];
    return J;
end );

BindGlobal( "MyMinimum", function( list )
    local d, i;
    d := list[1];
    for i in [2..Length(list)] do
        if AbsInt(list[i]) < AbsInt(d) then 
            d := list[i];
        elif AbsInt(list[i])=AbsInt(d) and d < 0 and list[i] > 0 then 
            d := list[i];
        fi;
    od;
    return d;
end );

BindGlobal( "Pivot", function(J, d, k, units)
    local v, m, j, l;

    # cut out relevant part
    v := J{[k..Length(J)]}[d];
    m := Unique(Filtered(v, x -> x <> 0*x));

    # check for ints
    l := Filtered(m, IsInt);
    if Length(l) > 0 then 
        l := MyMinimum(l);
        j := Position(v,l);
        return Swapped(J, k, j+k-1, d); 
    fi;

    # then check for w
    l := Filtered(m, x -> IsLiePUnit(units, x));
    if Length(l) > 0 then 
        l := l[1];
        j := Position(v,l);
        return Swapped(J, k, j+k-1, d); 
    fi;

    # finally print a statement
    Print("need invertible in ", v,"\n");
    return Random(m);
end );

BindGlobal( "MyBaseMat@", function(J, units)
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
        J := Pivot(J, d, k, units);
        if not IsList(J) then return J; fi;

        # clear out
        for i in [1..Length(J)] do
            if i <> k and J[i][d] <> 0*J[i][d] then
                J[i] := J[i] - J[i][d]*J[k];
            fi;
        od;
    until false;
end );

BindGlobal( "FactorSpace@", function( d, sub )
    return IdentityMat(d){Difference([1..d], List(sub, DepthVector))};
end );
