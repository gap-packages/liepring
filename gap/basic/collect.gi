##
## quotient remainder case 1: 
##      a is a rat fun in Z[w, x_1, ..., x_m] and p is prime
## return a list [quotient, remainder] with a = quotient*p + remainder
##      and remainder is a polynomial with coeffs in [0..p-1]
##
BindGlobal( "QuotRemInt", function( a, p, flag )
    local u, v, e, b, r, i, l;
    
    if IsRat(a) then 

        # catch info
        u := NumeratorRat(a);
        v := DenominatorRat(a);

        # compute
        if flag=true then 
            return [0*a, a mod p];
        elif v = 1 then 
            return QuotientRemainder(u,p);
        else
            Error("QuotRem: division case 1");
        fi;

    elif IsRationalFunction(a) then 

        # catch info
        u := SplitRatFun(a); v := u[2]; u := u[1];

        # compute
        if flag=true then 
            return [0*a, a mod p];
        elif v = 1 then 
            e := ExtRepPolynomialRatFun(u);
            b := ShallowCopy(e); 
            r := ShallowCopy(e);
            for i in [2,4..Length(e)] do
                l := QuotientRemainder( e[i], p );
                b[i] := l[1]; r[i] := l[2];
            od;
            b := PolynomialByExtRep(FamilyObj(a), b);
            r := PolynomialByExtRep(FamilyObj(a), r);
            return [b,r];
         else
            Error("QuotRem: division case 2");
         fi;

    else
        Error("case not found"); 
    fi;
end );

##
## quotient remainder case 2: 
##      a is a rat fun in Q[x_1, ..., x_m, w][p] and p is indeterminate
## return a list [quotient, remainder] with a = quotient*p + remainder
##      and remainder is a polynomial with leading coeff in [0..p-1]
##
BindGlobal( "QuotRemPoly", function( f, p, flag )
    local u, v, c;

    if IsRat(f) or 0*f=f then return [0*f, f]; fi;
    u := SplitRatFun(f); v := u[2]; u := u[1];

    if flag=true then 
        c := PolynomialCoefficientsOfPolynomial(u,p)[1];
        return [0*c, c/v];
    elif v=v^0 then 
        c := PolynomialCoefficientsOfPolynomial(u,p)[1];
        return [(u-c)/p, c];
    else 
        c := PolynomialCoefficientsOfPolynomial(u,p)[1];
        return [(u-c)/p, c]/v;
    fi;

end );

##
## extract entries from SC table
##

BindGlobal( "GetEntryMult", function( SC, i, j )
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
end );

BindGlobal( "GetEntryPow", function( SC, i )
    local k;
    k := Sum([1..i]);
    if IsBound(SC.tab[k]) then 
        return SC.tab[k];
    else
        return [];
    fi;
end );

##
## exp is a list of coefficients [c1, ..., cn]
## a is rational
## word is a list [i1, c1, i2, c2, ...] with i1 an index in [1..n] and
##    c1 is typically an entry from the SCTable
##
BindGlobal( "AddWordToExp", function( exp, a, word )
    local i; 
    if a = 0 then return exp; fi;
    for i in [1,3..Length(word)-1] do
        if a = 1 then 
            exp[word[i]] := exp[word[i]] + word[i+1];
        else
            exp[word[i]] := exp[word[i]] + a * word[i+1];
        fi;
    od;
end );

## 
## apply zeros
##
BindGlobal( "LRApplyZeros", function( SC, exp )
    return List(exp, x -> RedPol(SC.ring.zeros, x));
end );

## 
## apply the relations p bi and zeros
##
BindGlobal( "LRReduceExp", function( SC, exp )
    local p, n, a, b, i, new, u, v;
    p := SC.prime;
    n := SC.dim;
    new := ShallowCopy(exp);
    for i in [1..n] do
        b := GetEntryPow( SC, i );
        if IsInt(p) then 
           a := QuotRemInt( new[i], p, (Length(b)=0) );
        else
           a := QuotRemPoly( new[i], p, (Length(b)=0) );
        fi;
        if a[1] <> 0 * a[1] then AddWordToExp( new, a[1], b ); fi;
        new[i] := RedPol( SC.ring.zeros, a[2] );
    od;

    return MakeInt(new);
end );

##
## determine a reduced word
##
BindGlobal( "LRCollectWord", function( SC, word )
    local exp;
    exp := List([1..SC.dim], x -> 0);
    AddWordToExp( exp, 1, word );
    return LRReduceExp( SC, exp );
end );

##
## multiply
##
BindGlobal( "LRMultiply", function( SC, exp1, exp2 )
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
end );


