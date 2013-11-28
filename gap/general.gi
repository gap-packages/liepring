
#DepthVector := function ( vec )
#    local  i;
#    for i  in [ 1 .. Length( vec ) ]  do
#        if vec[i] <> 0 * vec[i]  then
#            return i;
#        fi;
#    od;
#    return Length( vec ) + 1;
#end;

WordByExps := function ( exp )
    local  w, i;
    w := [  ];
    for i  in [ 1 .. Length( exp ) ]  do
        if exp[i] <> 0  then
            Add( w, i );
            Add( w, exp[i] );
        fi;
    od;
    return w;
end;

MyCutVector := function ( vec, l )
    local  d, new, i;
    if Length( vec ) = 0  then
        return [  ];
    fi;
    d := Length( vec ) / l;
    new := [  ];
    for i  in [ 1 .. l ]  do
        Add( new, vec{[ d * (i - 1) + 1 .. d * i ]} );
    od;
    return new;
end;

ExpsByWord := function( n, w )
    local v, i;
    v := List([1..n], X -> 0);
    for i in [1,3..Length(w)-1] do
        v[w[i]] := v[w[i]] + w[i+1];
    od;
    return v;
end;

ValueRatFun := function( f, para, vals )
    local a, b;
    a := NumeratorOfRationalFunction(f);
    b := DenominatorOfRationalFunction(f);
    return Value( a, para, vals )/ Value( b, para, vals );
end;



