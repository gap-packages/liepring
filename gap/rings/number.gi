
BindGlobal( "ElmNumberForm", function( pp, polys )
    local f, e, p, a, b, i;

    f := polys[1];
    e := ExtRepPolynomialRatFun(f);
    p := IndeterminateByName("p");

    # case ac (xy) = 0
    if Length(e) = 2 then return 2*p-1; fi;

    # case ac (xy) + e = 0 with e <> 0
    if Length(e) = 4 and [] in e then return p-1; fi;

    # case ac (xy) + ad (x) = 0 
    if Length(e) = 4 and not ([] in e) then return 2*p-1; fi;

    # case ac (xy) + ad (x) + e = 0 with e <> 0 then 
    if Length(e) = 6 and [] in e then return p-1; fi;

    # general case
    a := 1; b := 0;
    for i in [1,3..Length(e)-1] do
        if Length(e[i])=0 then 
            b := e[i+1];
        elif Length(e[i])=2 then 
            a := a * e[i+1];
        elif Length(e[i])=4 then 
            a := a / e[i+1];
        fi;
    od;
 
    if b=a then return 2*p-1; else return p-1; fi;
end );
    
BindGlobal( "ElmNumberLin", function(pp, polys) 
    local r, s, v, A, t, i, e, j, b, p, k, a, w, W;

    # set up
    r := Length(pp);
    s := Length(polys);
    v := List(pp, x -> IndeterminateNumberOfLaurentPolynomial(x));
    A := NullMat(r,s);
    t := List([1..s], x -> 0);
    w := ExtRepPolynomialRatFun(IndeterminateByName("w"))[1][1];
    W := IndeterminateByName("w");

    # a trivial case
    if r = 0 and s = 0 then return 1; fi;
    if r = 0 and s > 0 then return 0; fi;

    # determine matrix A and vector t
    for i in [1..s] do
        e := ExtRepPolynomialRatFun(polys[i]);
        for j in [1,3..Length(e)-1] do
            if Length(e[j]) = 0 then 
                t[i] := t[i] - e[j+1];
            elif Length(e[j])=2 and e[j][1] = w then 
                t[i] := t[i] - W^e[j][2]*e[j+1];
            elif Length(e[j])=4 and e[j][1] = w then
                k := Position(v, e[j][3]);
                A[k][i] := A[k][i] + W^e[j][2]*e[j+1]; 
            elif Length(e[j])=4 and e[j][3] = w then
                k := Position(v, e[j][1]);
                A[k][i] := A[k][i] + W^e[j][4]*e[j+1]; 
            elif Length(e[j])=2 then 
                k := Position(v, e[j][1]);
                A[k][i] := A[k][i] + e[j+1]; 
            else
                Error("should not happen");
            fi;
        od;
    od;

    # determine number of solutions
    A := A*W^0; t := t*W^0;
    a := SolutionMat(A, t);
    if IsBool(a) then return 0; fi;
    b := Length(TriangulizedNullspaceMat(A));
    p := IndeterminateByName("p");
    return p^b;
end );

BindGlobal( "ElmNumberUni", function(pp, polys)
    local r, s, f, i, sub, res, p, w, d, x, g4, s16, y, c, e, z, v;

    # set up
    r := Length(pp);
    s := Length(polys);
    p := IndeterminateByName("p");
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    v := IndeterminateByName("v");
    w := IndeterminateByName("w");
    g4 := IndeterminateByName("(p-1,4)");
    s16 := IndeterminateByName("(p^2-1,16)");

    # a trivial case
    if r = 0 and s = 0 then return 1; fi;
    if r = 0 and s > 0 then return 0; fi;

    # split polys according to vars
    sub := List(pp, x -> Filtered(polys, f -> VarsOfPoly(f)=[x])); 
    res := 1;
    for i in [1..Length(sub)] do
        if Length(sub[i]) = 0 then 
            res := res * p; 
        else
            f := Gcd(sub[i]);
            c := List(Collected(Factors(f)), x -> x[1]);
            e := List(c, DegreeOfPoly);
            if Length(c)=1 and c[1] = c[1]^0 then 
                d := 0;
            elif ForAll(e, x -> x <= 1) then 
                d := Length(c);
            elif f = x^2-2 or f = w^2-1/2*x^2 then 
                d := (s16-8)/4;
            elif f = w*x^2-2 or f = w^3-1/2*x^2 then 
                d := 2-(s16-8)/4;
            elif f = x^2+1 then 
                d := g4-2;
            elif f = w^3-x^2 then 
                d := 0;
            elif f = x^2+w then 
                d := 4-g4;
            elif f = z^3-w*z or f = w^3*z-z^3 then 
                d := 1;
            else
                Error("cannot find number of roots of poly");
            fi;
            res := res * d;
        fi;
    od;
    return res;
end );

DeclareGlobalFunction( "ElmNumberMon" );
InstallGlobalFunction( "ElmNumberMon", function(pp, polys)
    local f, a, b, c, r, p;

    p := IndeterminateByName("p");

    # catch trivial case
    if Length(polys) = 0 then return p^Length(pp); fi;

    f := SquareFreeGB(polys);

    # these cases should not happen
    if Length(pp) = 0 then return 0; fi;
    if f[1] = f[1]^0 then return 0; fi;


    # split 
    a := pp[1];
    b := pp{[2..Length(pp)]};
    c := Filtered(polys, x -> IsPolynomial(x/a));
    r := Difference(polys, c);

    # this is the base case
    if Length(pp) = 1 and c = [a] then 
        return 1;
    elif Length(pp) = 1 then 
        return p;
    fi;

    # recuse
    if Length(c) = 0 then 
        return p*ElmNumberMon(b, r);
    elif c = [a] then 
        return ElmNumberMon(b,r);
    else
        return ElmNumberMon(b,r)+(p-1)*ElmNumberMon(b, Concatenation(r,c/a));
    fi;
end );

BindGlobal( "IsLinearSystem", function( pp, polys )
    local ee;
    ee := List(polys, x -> DegreeOfPoly(x));
    return ForAll(ee, x -> x <= 1);
end );

BindGlobal( "IsUnivarSystem", function( pp, polys )
    local ee;
    ee := List(polys, x -> VarsOfPoly(x));
    return ForAll(ee, x -> Length(x) = 1);
end );

BindGlobal( "IsFormSystem", function( pp, polys )
    local e, t, i, j;
    if Length(polys) <> 1 or Length(pp) <> 2 then return false; fi;
    e := ExtRepPolynomialRatFun(polys[1]);
    t := false;
    for i in [1,3..Length(e)-1] do
        if Length(e[i]) > 4 then return false; fi;
        for j in [2,4..Length(e[i])] do
            if e[i][j] <> 1 then return false; fi;
        od;
        if Length(e[i]) = 4 then t := true; fi;
    od;
    return t;
end );

BindGlobal( "IsMonomSystem", function(pp, polys)
    local f, e;
    for f in polys do
        e := ExtRepPolynomialRatFun(f);
        if Length(e)>2 then return false; fi;
        if e[2] <> 1 then return false; fi;
        if Length(e[1])=0 then return false; fi;
        if ForAny(e[1]{[2,4..Length(e[1])]},x->x<>1) then return false; fi;
    od;
    return true;
end ); 

BindGlobal( "NumberOfZeros", function( pp, polys )
    local w, f, S, p, x, y, z, t, u, g3, g4, s16, v, r, s;

    # check
    p := IndeterminateByName("p");
    w := IndeterminateByName("w");
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    t := IndeterminateByName("t");
    u := IndeterminateByName("u");
    v := IndeterminateByName("v");
    s := IndeterminateByName("s");
    r := IndeterminateByName("r");
    g3 := IndeterminateByName("(p-1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    s16 := IndeterminateByName("(p^2-1,16)");
    if Length(polys) = 0 then return p^Length(pp); fi;

    # 1. step : groebner
    f := SquareFreeGB(polys);

    # 2. step : check cases

    # no zeros
    if 1 in f or p^0 in f or w in f or w-1 in f or w+1 in f or
       4*w-1 in f or -x^2+w in f or -v^2+w in f or v^2-w in f or
       -z^2+w in f or z^2-w in f then 
        return 0; 
    fi;

    # single case
    if Length(pp) = 1 then 
        if f = [x^2-2,w-2] then return 0; fi;
    fi;

    # quadratic cases
    if Length(pp) = 2 then
        if f = [y^2-y, x^2-x*y-x+y, w-y] then return 0; fi;
        if f = [x*y-x, x^2-x, w-x] then return 0; fi;
        if f = [x*y+x, x^2-x, w-x] then return 0; fi;
        if f = [x^2+y+1, w+y+1] then return 0; fi;
        if f = [x^2+y^2+y, x^2+w+y] then return 0; fi;
        if f = [2*y^2+y, 2*x^2+y, 2*w+y] then return 0; fi;
        if f = [2*x^2+y, 2*w+y] then return 0; fi; 
        if f = [y^2-3, x+y-1, w+2*y-4] then return 0; fi;
        if f = [x+y-1, -y^2+w+2*y-1] then return 0; fi;
        if f = [2*y^2+x, 2*w+x] then return 0; fi;
        if f = [x*(y-y^2)-1] or f = [-x*(y-y^2)+1] then return p-2; fi; 
        if f = [x^2-x*y+y] then return p-1; fi; 
        if f = [x*(y-1), x*(x-1)] then return p+1; fi;
        if f = [x*(y+1), x*(x-1)] then return p+1; fi;
        if f = [y^2-1,x-y] then return 2; fi;
        if f = [y^2-1,w^2*x^2-y] then return g4; fi;
        if f = [x^2+y+1] then return p; fi;
        if f = [x^2+w+y] then return p; fi;
        if f = [y^2-3*y+3,x+y-4] then return g3-1; fi;
        if f = [y^2-3,x+y-1] then return 2 - (g3-g4+1)^2/2; fi;
        if f = [(x-y)*(y+1)] then return 2*p-1; fi;
        if f = [(w*x-y)*(y+1)] then return 2*p-1; fi;
    fi;

    # triple cases
    if Length(pp) = 3 then 
        if f = [-2*z^3+2*w*z-x] then return p^2; fi;
        if f = [2*w^3*z-w^2*x-2*z^3] then return p^2; fi;
        if f = [z^2-z, y*z-y, y^2-y, x, w-y] then return 0; fi;
        if f = [-z^2+y, -y*z+x, w-y] then return 0; fi;
        if f = [z^4-1, -z^2+y, -y*z+x, w-y] then return 0; fi;
        if f = [y*z^2-1, -y*z+x, w-y] then return 0; fi;
        if f = [-z^2+y, x, w-y] then return 0; fi;
        if f = [-z^2+y, w-y] then return 0; fi;
        if f = [z^5-z, -z^2+y, -y*z+x, w-y] then return 0; fi;
        if f = [-z^2+y, x^2-y, w-y] then return 0; fi;
        if f = [y^2*z^2-y, -y*z+x, w-y] then return 0; fi;
        if f = [x^2-y, w-y] then return 0; fi;
        if f = [-y*z+x, -(y-1)*(y-w)] then return 2*p; fi;
        if f = [y^3*z^2-y^2*z^2-y^2+y, -y*z+x, -x^2+w] then return 0; fi;
        if f = [y*z, x, y^2*z+w*z-y*z, w*y-y^2-w+y] then return 2; fi;
        if f = [z, x+y-1, -y^2+w+2*y-1] then return 0; fi;
        if f = [y*z^2-y+1, -y*z^2+y*z+x+y-1] then return p-2; fi;
        if f = [y*z^2-y+1, -y*z^2+y*z+x+y-1, -y^2+w+y] then return 0; fi;
        if f = [y*z^3-y*z+z, -y*z^2+y*z+x+y-1] then return 2*p-3; fi;
        if f = [(z+1)*(x*y-z)] then return 2*p^2-p+1; fi;
        if f = [z+1, x*y+1] then return p-1; fi;
        if f = [x*y-z] then return p^2; fi;
        if f = [-z^3+w*z-1/2*x] then return p^2; fi;
        if f = [w^3*z-1/2*w^2*x-z^3] then return p^2; fi;
        if f = [-y*z+x, w-y] then return p; fi;
        if f = [y*z, x, w-y] then return 1; fi;
        if f = [(y-1)*((x-y)*(z-1)-1)] then return p^2 + (p-1)^2; fi;
        if f = [y, x*z-x-1] then return p-1; fi;
        if f = [z, x*y-y^2-x+2*y-1] then return 2*p-1; fi;
        if f = [ y-1, x*z-x-z ] then return p-1; fi;
        if f = [(x-y)*(z-1)-1] then return p*(p-1); fi;
        if f = [y*z+y-1, x] then return p-1; fi;
        if f = [(z+1)*(x+y)-1] then return p*(p-1); fi;
        if f = [y, x*(z+1)-1] then return p-1; fi;
        if f = [z^2+2*z+2,y,x+z+1] then return g4-2; fi;
        if f = [y-1, (x+1)*(z+1)-1] then return p-1; fi;
        if f = [z*(y*z+y-1), (y-1)*(y*z+y-1), y*z+x+y-1] then return p-1; fi;
        if f = [y*z^2-y+1, y*z+x] then return p-2; fi;
        if f = [(y*z+y-1)*(y*z^2-y+1), 
                -y^2*z^2+y^2+y*z+x-y] then return 2*p-4; fi;
        if f = [(y*z+y-1)*(y*z-1/2*z^2+y-z-1),  
                2*y^2*z-y*z^2+2*y^2-2*y*z+x-3*y+z+1] then return 2*p-4; fi;
    fi;

    # cases with 4 params
    if Length(pp) = 4 then 
        if f = [y,x*t+t^2] then return 2*p^2-p; fi;
        if f = [z,x*t+t^2] then return 2*p^2-p; fi;
        if f = [z,y,x*t+t^2] then return 2*p-1; fi;
        if f = [z,y,x*t*(x+t)] then return 3*p-2; fi;
        if f = [z,y,x*t*(x+t)] then return 3*p-2; fi;
        if f = [t*(z-t)] then return p^2*(2*p-1); fi;
        if f = [t*(z-t),t*(x*t-y*t-z+t),z*(x-1)+t*(1-y)] then 
               return 3*p^2-p; fi;
        if f = [t*(z-t), t*(x*t-y*t-z+t), x*z-y*t-z+t, 
           t*(x*y-y^2-x+y), x^2-x*y-x+y] then return 2*p^2-1; fi; 
        if f = [t*(z-t), t*(y-1), x-1] then return p^2+p-1; fi;
        if f = [t*(z-t), y-1, x-1] then return 2*p-1; fi;
        if f = [t*(z-t), (y-1)*(z-t), x-y] then return p^2+p-1; fi;

        if f = [x*y-z*t] then return p^3+p^2-p; fi;
        if f = [x*t-y*z] then return p^3+p^2-p; fi;
        if f = [w*y*t-x*z] then return p^3+p^2-p; fi;
        if f = [(x-1)*z-(y-1)*t] then return p^3+p^2-p; fi;

        if f = [t,z*(x-1)] then return p*(2*p-1); fi;
        if f = [z-t,(x-y)*t] then return p*(2*p-1); fi;
        if f = [z,x*(x+t)] then return p*(2*p-1); fi;
        if f = [z,y,x*(x+t)] then return (2*p-1); fi;

        if f = [y,x*z,x*(x+t)] then return p^2+p-1; fi;

        if f = [z+t,x*y+t^2] then return p^2; fi;
        if f = [y*z+t^2,x+t] then return p^2; fi;
        
        if f = [(x+t)*(x*t-y*z)] then return 2*p^3-p; fi;

        if f = [z,t*x*(x+t)] then return p*(3*p-2); fi;
        if f = [y,x*t*(x+t)] then return p*(3*p-2); fi;

        if f = [x*t-y*z,w*y*t-x*z,z*(w*y^2-x^2)] then return 2*p^2-1; fi;
        if f = [(x+t)*(x*t-y*z),w*y*t-x*z,(w*y^2-x^2)-(x*t-y*z)] then 
            return 2*p^2-1; fi;
        if f = [t*(z-t), t*(x-y), (x-1)*z-(y-1)*t, (x-1)*(x-y)] then 
            return 2*p^2-1; fi;

        if f = [y*z+t^2,x+t,w*t^2-z^2,w*y+z] then return 1; fi; 
        if f = [x*t-y*z,w*y*t-x*z,w*y^2-x^2] then return p^2; fi;
    fi;

    if Length(pp) = 7 then 
        if f = [w*x-s*v] then return p^6; fi;
        if f = [-z*s+y,w*x-s*v] then return p^5; fi;
        if f = [u,-z*s+y,w*x-s*v] then return p^4; fi;
        if f = [u,-z*s+y,x*s+2*t-v,1/2*s^2*v+w*t-1/2*w*v,w*x-s*v] then
            return p^3; fi;
        if f = [u,z+v,s*v+y,x*s+2*t-v,1/2*s^2*v+w*t-1/2*w*v,w*x-s*v] then
             return p^2; fi;
        if f = [v,u,z,y,x*s+2*t,w*t,w*x] then return p; fi;
        if f = [u,t,z+v,s*v+y,x*s-v,-s^2*v+w*v,w*x-s*v] then return p; fi;
    fi;

#############################################################################
##
## Bis zu 4 Parametern ist die Liste der Polynome vollstaendig und geprueft.
## Ebenso fuer 7 Parameter. Alle anderen Faelle 5,6,8,12 sind unvollstaendig
## und auch nicht unbedingt strikt ueberprueft.
##
#############################################################################

    if Length(pp) = 5 then
        if f = [t, x*u-1] then return p^2*(p-1); fi;
        if f = [t, (x*u-1)*(x*y+z)] then return 2*p^3-2*p^2+p; fi;
        if f = [t,z*u+y,x*u-1] then return p*(p-1); fi;
        if f = [t,x*y+z] then return p^3; fi;
        if f = [u,t,x*y+z] then return p^2; fi;
        if f = [x*u+t*u-1] then return p^3*(p-1); fi;
        if f = [t*(y-z*u), u*(x+t)-1] then return p*(p-1)*(2*p-1); fi;
        if f = [x*u+1/2*t*u-1] then return p^3*(p-1); fi;
        if f = [t*u-2,x] then return p^2*(p-1); fi;
        if f = [1/2*z*t*u^2+(x*u+1/2*t*u-1)*y] then 
                       return p^2*(p-1)^2 + p^3 + p*(2*p-1)*(p-1); fi;
        if f = [(-z*t*u^2+y*t*u-2*z*u-2*y)*t, x*t*u^2+t^2*u^2+2*x*u-t*u-2, 
                 1/2*z*t*u^2+x*y*u+1/2*y*t*u-y, (z*t*u+x*y+z)*t ] then 
                       return 2*p*(p-1)^2; fi;
        if f = [y*u*v-z^2*t] then return p*(p^3+2*p^2-3*p+1); fi;
        if f = [y*u*v-z^2*t, w*z*t] then return (2*p-1)*(3*p^2-3*p+1); fi;
        if f = [u,t,z+v,y,v^2+w] then return 0; fi;
        if f = [w*y*z*t-1/2*y*u*v+1/2*z^2*t] then return 6*p^3-9*p^2+5*p-1; fi;
        if f = [u,t,y,2*z^2-v^2+w] then return 0; fi;
        if f = [v,u,z*t,w*y+1/2*z] then return 2*p-1; fi;
        if f = [v,z*t,w*y*u+2*z*u] then return (2*p-1)^2; fi;
    fi;

    if Length(pp) = 6 then 
        if f = [y*u-z*t,z*(x*t-u*v),(x+t)*(x*t-u*v)] then 
                                return p*(2*p^3+2*p^2-4*p+1); fi;
        if f = [u,z*t,x+t] then return p^2*(2*p-1); fi;
        if f = [y*u-z*t] then return p^2*(p^3+p^2-p); fi;
        if f = [y*u-z*t,x+t] then return p*(p^3+p^2-p); fi;
        if f = [u,t,x*y-z*v] then return p^3+p^2-p; fi;
    fi;

    if Length(pp) = 8 then 
        if f = [x*r-s*v] then return p^4*(p^3+p^2-p); fi;
        if f = [v,r+u,y*s+z*u,x*u,x*t] then return p*(p^3+2*p^2-4*p+1); fi;
        if f = [s,r,x*t-u*v] then return p^2*(p^3+p^2-p); fi;
        if f = [u,s,r,t,x*y-z*v] then return (p^3+p^2-p); fi;
        if f = [v,r+u,y*x+z*u,x*u,x*t] then return p^3*(2*p-1)+p^2*(p-1); fi;
    fi;

    # generic cases
    if IsFormSystem(pp, f) then return ElmNumberForm(pp,f); fi;
    if IsLinearSystem(pp, f) then return ElmNumberLin(pp,f); fi;
    if IsUnivarSystem(pp, f) then return ElmNumberUni(pp,f); fi;
    if IsMonomSystem(pp, f) then return ElmNumberMon(pp,f); fi;

    # bad case
    #Print("cannot solve ",f,"\n");

#    Error("NumberOfZeros cannot find the answer");

    # cannot find the answer
    return fail;
end );

BindGlobal( "ElementNumber", function( pp, units, zeros )
    local d, x, y, u, t, z, r, U, v, s, W, R, T, i, w, S, p;

    d := Length(pp);
    if d = 0 then return 1; fi;

    w := IndeterminateByName("w");
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    t := IndeterminateByName("t");
    v := IndeterminateByName("v");
    u := IndeterminateByName("u");
    s := IndeterminateByName("s");
    r := IndeterminateByName("r");

    U := units;
    U := Filtered(U, a -> a <> x^2-w); 
    U := Filtered(U, a -> a <> -x^2+w); 
    U := Filtered(U, a -> a <> y^2-w); 
    U := Filtered(U, a -> a <> -y^2+w); 
    U := Filtered(U, a -> a <> z^2-w); 
    U := Filtered(U, a -> a <> -z^2+w); 
    U := Filtered(U, a -> a <> t^2-w); 
    U := Filtered(U, a -> a <> -t^2+w); 
    U := Filtered(U, a -> a <> v^2-w); 
    U := Filtered(U, a -> a <> -v^2+w); 
    U := Filtered(U, a -> a <> u^2-w); 
    U := Filtered(U, a -> a <> -u^2+w); 
    U := Filtered(U, a -> a <> s^2-w); 
    U := Filtered(U, a -> a <> -s^2+w); 
    U := Filtered(U, a -> a <> r^2-w); 
    U := Filtered(U, a -> a <> -r^2+w); 

    t := Length(U);
    s := Length(zeros);
    p := IndeterminateByName("p");
    if t=0 and s=0 then return p^d; fi;

    W := Combinations([1..t]);
    R := List(W, x -> 0);
    T := 0;
    for i in [1..Length(W)] do
        w := W[i];
        S := Concatenation(zeros, U{w});
        R[i] := NumberOfZeros( pp, S );
        #Print(W[i]," yields ",R[i],"\n");
        if R[i]=fail then return fail; fi;
        T := T + (-1)^Length(w) * R[i];
    od;
    
    return T;
end );

BindGlobal( "ElementNumbers", function( pp, s )
    local t, e, p, j, k, f;

    # set up and catch trivial case
    t := Set(List(s, x -> x.norm));
    p := IndeterminateByName("p");
    e := List(t, x -> 0*p^0);
    if Length(t) = 1 then 
        return rec( norms := t, numbs := [p^Length(pp)]); 
    fi;

    # go through cases
    for j in [1..Length(s)] do
        k := Position(t, s[j].norm);
        f := ElementNumber(pp, s[j].units, s[j].zeros);
        if f = fail or e = fail then
            return fail;
        else
            e[k] := e[k] + f;
        fi;
    od;

    # filter
    j := Filtered([1..Length(e)], x -> e[x] <> 0*e[x]);
    e := e{j}; t := t{j};

    # sort 
    SortParallel(e, t);
   
    return rec( norms := t, numbs := e);
end );

