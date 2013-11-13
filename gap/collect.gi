
QuotRemInt := function( a, p )
    local q, r;
    q := 0; r := a;
    while not r in [0..p-1] do
        if r < 0 then 
            q := q-1;
            r := r + p; 
        else
            q := q+1;
            r := r - p; 
        fi;
    od;
    return [q,r];
end;

ReduceExtRep := function( e )
    local I, i;
    I := [];
    for i in [2,4..Length(e)] do
        if e[i] <> 0 then Add(I, i-1); Add(I, i); fi;
    od;
    return e{I};
end;

QuotRemPoly := function( f, p )
    local F, e, a, r, i, l, k, c;

    # if p is a prime, then reduce coeffs of f
    if IsInt(p) then 
        F := FamilyObj(f);
        e := ExtRepPolynomialRatFun(f);
        a := ShallowCopy(e); 
        r := ShallowCopy(e);
        for i in [2,4..Length(e)] do
            l := QuotRemInt( e[i], p );
            a[i] := l[1]; r[i] := l[2];
        od;
        a := ReduceExtRep(a); r := ReduceExtRep(r);
        return [PolynomialByExtRep(F, a), PolynomialByExtRep(F, r)];

    # otherwise p is indeterminate 
    else
        c := PolynomialCoefficientsOfPolynomial( f, p );
        if Length(c) <= 1 then  # f is constant
            return [0*f, f];
        else
            a := Sum(List([2..Length(c)], x -> c[x]*p^(x-2)));
            return [a, c[1]];
        fi;
    fi;
end;

GetEntryMult := function( SC, i, j )
    local k, e;
    if i = j then return [1,[]]; fi;
    if i < j then 
        k := Sum([1..j-1])+i;
        e := -1;
    elif i > j then 
        k := Sum([1..i-1])+j;
        e := 1;
    fi;
    if IsBound(SC.tab[k]) then 
        return [e, SC.tab[k]]; 
    else
        return [e, []];
    fi;
end;

GetEntryPow := function( SC, i )
    local k;
    k := Sum([1..i]);
    if IsBound(SC.tab[k]) then 
        return SC.tab[k];
    else
        return [];
    fi;
end;

AddWordToExp := function( exp, a, word )
    local i; 
    if a = 0 then return exp; fi;
    for i in [1,3..Length(word)-1] do
        if a = 1 then 
            exp[word[i]] := exp[word[i]] + word[i+1];
        else
            exp[word[i]] := exp[word[i]] + a * word[i+1];
        fi;
    od;
end;

LRReduceExp := function( SC, exp )
    local p, n, a, b, i, new, u, v;
    p := SC.prime;
    n := SC.dim;
    new := ShallowCopy(exp);
    for i in [1..n] do
        if IsInt(new[i]) and IsInt(p) then 
            if new[i] < 0 or new[i] >= p then     
                a := QuotRemInt( new[i], p );
                b := GetEntryPow( SC, i );
                AddWordToExp( new, a[1], b );
                new[i] := a[2];
            fi;
        elif IsRat(new[i]) and IsInt(p) then 
            u := NumeratorRat(new[i]);
            v := DenominatorRat(new[i])^-1 mod p;
            if i = n or (u>=0 and u<p) then 
                new[i] := u*v mod p;
            else
                Error("check this");
                a := QuotRemInt( u, p );
                b := GetEntryPow( SC, i );
                AddWordToExp( new, (a[1] * v) mod p, b );
                new[i] := (a[2] * v) mod p;
            fi;
        elif IsPolynomial(new[i]) then 
            a := QuotRemPoly( new[i], p );
            if a[1] <> 0*a[1] then 
                b := GetEntryPow( SC, i );
                AddWordToExp( new, a[1], b );
                new[i] := a[2];
            fi;
        elif IsRationalFunction(new[i]) then 
            u := NumeratorOfRationalFunction(new[i]);
            v := DenominatorOfRationalFunction(new[i]);
            a := QuotRemPoly( u, p );
            if a[1] <> 0*a[1] then 
                b := GetEntryPow( SC, i );
                AddWordToExp( new, (a[1]/v), b );
                new[i] := a[2]/v;
            fi;
        fi;
    od;
    return new;
end;

LRCollectWord := function( SC, word )
    local exp;
    exp := List([1..SC.dim], x -> 0);
    AddWordToExp( exp, 1, word );
    return LRReduceExp( SC, exp );
end;

LRMultiply := function( SC, exp1, exp2 )
    local exp, k, i, j, a, b;
  
    exp := List([1..SC.dim], x -> 0);
    for i in [1..SC.dim] do
        for j in [1..SC.dim] do
            if exp1[i] <> 0 and exp2[j] <> 0 then 
                a := GetEntryMult( SC, i, j );
                b := a[1]*exp1[i]*exp2[j];
                AddWordToExp( exp, b, a[2] );
            fi;
        od;
    od;

    return LRReduceExp( SC, exp );
end;


