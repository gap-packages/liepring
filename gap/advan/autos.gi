BindGlobal( "ExtendAuto", function(L, mat)
    local p, dim, def, bas, img, i, j, k, h, e, new;
    p := PrimeOfLiePRing(L);
    dim := DimensionOfLiePRing(L);
    def := LPRDefs(L);
    bas := BasisOfLiePRing(L);
    img := [];
    for i in [1..dim] do
        if IsBool(def[i]) then 
           img[i] := mat[i]*bas;
        elif def[i][1] = def[i][2] then
           j := def[i][1]; 
           e := Exponents(p*bas[j]);
           img[i] := p*img[j];
           for k in [1..dim] do
               if e[k] <> 0*e[k] and k < i then 
                   img[i] := img[i] - e[k]*img[k];
               fi;
           od;
           img[i] := e[i]^-1 * img[i];
        else
           j := def[i][1];
           h := def[i][2];
           e := Exponents(bas[j]*bas[h]);
           img[i] := img[j]*img[h];
           for k in [1..dim] do
               if e[k] <> 0*e[k] and k < i then 
                   img[i] := img[i] - e[k]*img[k];
               fi;
           od;
           img[i] := e[i]^-1 * img[i];
        fi;
    od;
    return List(img, Exponents);
end );

BindGlobal( "MyVals", function( pol, vars, vals )
    if IsBool(pol) or IsInt(pol) then return pol; fi;
    return Value(pol, vars, vals);
end );

BindGlobal( "EliminateVars", function( pol, eli )
    local i;
    for i in [1..Length(eli)] do
        pol := MyVals( pol, [eli[i][1]], [eli[i][2]]);
    od;
    return pol;
end );
   
BindGlobal( "TryEliminate", function( pol )
    local x, q;
    if pol = 0*pol then return fail; fi;
    for x in VarsOfPoly(pol) do
        q := PolynomialCoefficientsOfPolynomial(pol, x);
        if Length(q) = 2 and q[2]^2 = q[2]^0 then 
            return [x, -q[1]/q[2]];
        fi;
    od;
    return fail;
end );
 
BindGlobal( "CheckForElimination", function( tod, new, h )
    local rel, var, eli, don, i, j, e;
    
    rel := new;
    var := tod.var;
    eli := tod.eli;

    # loop
    don := false;
    while not don do 
        don := true;
        for i in [1..Length(rel)] do
            e := TryEliminate(rel[i]);
            if not IsBool(e) then 
                e := [e[1], e[2], h];
                don := false;
                rel[i] := false;
                Add(eli, e);
                var := Filtered(var, x -> x <> e[1]);
                for j in [1..Length(rel)] do
                    rel[j] := MyVals(rel[j], [e[1]], [e[2]]);
                od;
                for j in [1..Length(eli)] do
                    eli[j][2] := MyVals(eli[j][2], [e[1]], [e[2]]);
                od;
            fi;
        od;
        if Length(var) = 0 then don := true; fi;
        rel := Filtered(rel, x -> x <> false );
        rel := Filtered(rel, x -> x <> 0*x);
        if Length(rel) = 0 then don := true; fi;
    od;

    # result
    tod.eli := eli; return rel;
end );

BindGlobal( "AddByExpo", function (k, rel, c, r)
    local i, j, f, d;
    for i in [1..Length(r)] do
        if r[i] <> 0*r[i] then 
            j := c[i];
            f := NumeratorOfRationalFunction(r[i]);
            d := Length(rel[j]);
            if not f in rel[j] then rel[j][d+1] := f; fi;
        fi;
    od;
end );

BindGlobal( "FindWord", function( tod, elm )
    local w, i, v;

    # eliminate
    w := EliminateVars( elm, tod.eli );
    if w = 0*w then return false; fi;

    # check
    for i in Reversed([1..Length(tod.rel)]) do
        v := PolynomialDivisionAlgorithm(w, tod.rel[i], ORDER);
        if v[1] = 0*v[1] then return v[2]; fi;
    od;

    # no solution
    return fail;
end );

BindGlobal( "AutGroupDescription", function(arg)
    local L, flag, p, b, d, a, q, c, k, D, A, i, j, img, rel, var, e, r, U, 
          don, new, tod, old, uni, wd, w, v, try, M, B, m, l,
          lname, auto1, auto2, auto3, auto4, eqns1, eqns2, eqns3, eqns4,
          comment1, comment2, comment3, comment4, x, y;

    # catch args
    L := arg[1];
    if Length(arg)=2 then 
        flag := arg[2];
    else
        flag := 0;
    fi;

    # set up
    p := PrimeOfLiePRing(L);
    b := BasisOfLiePRing(L);
    d := DimensionOfLiePRing(L);
    a := MinimalGeneratorNumberOfLiePRing(L);
    q := ParametersOfLiePRing(L);
    c := List(b, x -> Length(Factors(LiePOrder(x))));
    lname:=LibraryName(L);
# mrvl addition
# The defining generators never appear in the relations Sum r_i.b_i
    for i in [1..a] do
      c[i] := 0;
    od;
    k := Maximum(c);

    # variables
    w:=IndeterminateByName("w");
    x:=IndeterminateByName("x");
    y:=IndeterminateByName("y");
    D := Indeterminate(Rationals, "D");
    A := NullMat(a,d);
    m := NullMat(a,d);
    for i in [1..a] do
        for j in [1..d] do
            A[i][j] := Indeterminate(Rationals, Concatenation("A",
                        String(i), String(j)));
            m[i][j] := A[i][j];
            if flag=1 and i<>j and j<=a then A[i][j] := 0; fi;
            if flag=1 and i=j then A[i][j] := 1; fi;
            if flag=2 and j>a then A[i][j] := 0; fi;
        od;
    od;
    var := Concatenation([D], Flat(A));

#Difficult LiePRings

#1 7.487
if lname = "7.487" then

auto1:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns1:=[
        [-m[1][3]^2-m[1][3]*m[2][3]+m[1][4]-2*m[1][5]-m[2][5],
         m[1][3]*m[2][3]-m[1][4]+m[2][3]^2-2*m[2][4]+m[2][5]]];
comment1:="p^8 automorphisms";
#Determinant 1

auto2:=[
        [-1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,-1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns2:=[
        [m[1][3]^2+m[1][3]*m[2][3]+m[1][4]-2*m[1][5]-m[2][5],
         -m[1][3]*m[2][3]-m[1][4]-m[2][3]^2-2*m[2][4]+m[2][5]]];
comment2:="p^8 automorphisms";
#Determinant 1

auto3:=[
        [0,1,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [1,0,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns3:=[
        [m[1][3]^2+m[1][3]*m[2][3]-2*m[1][4]+m[1][5]-m[2][4],
         -m[1][3]*m[2][3]-m[1][5]-m[2][3]^2+m[2][4]-2*m[2][5]]];
comment3:="p^8 automorphisms";
#Determinant -1

auto4:=[
        [0,-1,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-1,0,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns4:=[
        [-m[1][3]^2-m[1][3]*m[2][3]-2*m[1][4]+m[1][5]-m[2][4],
         m[1][3]*m[2][3]-m[1][5]+m[2][3]^2+m[2][4]-2*m[2][5]]];
comment4:="p^8 automorphisms";
#Determinant -1

#Automorphism group has order 4*p^8

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2),
        rec(auto:=auto3,eqns:=eqns3,comment:=comment3),
        rec(auto:=auto4,eqns:=eqns4,comment:=comment4)];

fi;
#---------------------------------------

#2 7.488
if lname = "7.488" then

auto1:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns1:=[
        [-w*m[1][3]^2-2*w*m[1][5]-m[1][3]*m[2][3]+m[1][4]-m[2][5],
         w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5]+m[2][3]^2-2*m[2][4]]];
comment1:="p^8 automorphisms";
#Determinant 1

auto2:=[
        [-1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,-1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns2:=[
        [w*m[1][3]^2-2*w*m[1][5]+m[1][3]*m[2][3]+m[1][4]-m[2][5],
         -w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5]-m[2][3]^2-2*m[2][4]]];
comment2:="p^8 automorphisms";
#Determinant 1

auto3:=[
        [0,m[1][2],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [w,0,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns3:=[
        [2*w^2*m[1][2]*m[1][5]-w*m[1][2]^2*m[2][4]+w*m[1][2]*m[1][3]*m[2][3]-w*m[1][2]*m[1][4]+w*m[1][3]^2-w*m[1][5]-m[1][4],
         -w^2*m[1][2]*m[2][5]-w^2*m[1][5]-w*m[1][2]*m[2][3]^2+2*w*m[1][2]*m[2][4]-w*m[1][3]*m[2][3]-w*m[2][5]-m[2][4]],
        [w*m[1][2]-1]];
comment3:="p^8 automorphisms";
#Determinant -1

auto4:=[
        [0,m[1][2],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-w,0,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns4:=[
        [-2*w^2*m[1][2]*m[1][5]-w*m[1][2]^2*m[2][4]+w*m[1][2]*m[1][3]*m[2][3]+w*m[1][2]*m[1][4]-w*m[1][3]^2-w*m[1][5]-m[1][4],
         w^2*m[1][2]*m[2][5]-w^2*m[1][5]-w*m[1][2]*m[2][3]^2-2*w*m[1][2]*m[2][4]+w*m[1][3]*m[2][3]-w*m[2][5]-m[2][4]],
        [w*m[1][2]+1]];
comment4:="p^8 automorphisms";
#Determinant -1

#Automorphism group has order 4*p^8

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2),
        rec(auto:=auto3,eqns:=eqns3,comment:=comment3),
        rec(auto:=auto4,eqns:=eqns4,comment:=comment4)];

fi;
#---------------------------------

#3 7.489
if lname = "7.489" then

auto1:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns1:=[
        [x*m[1][3]^2+2*x*m[1][5]-m[1][3]*m[2][3]+m[1][4]-m[2][5],
         -x*m[1][3]*m[2][3]+x*m[1][4]-x*m[2][5]+m[2][3]^2-2*m[2][4]]];
comment1:="p^8 automorphisms";

auto2:=[
        [0,m[1][2],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-x,0,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns2:=[
       [x^2*m[1][2]*m[1][3]+2*x^2*m[1][2]*m[1][5]+x*m[1][2]^2*m[2][4]-x*m[1][2]*m[1][3]*m[2][3]+x*m[1][2]*m[1][4]-x*m[1][3]^2+x*m[1][3]+x*m[1][5]-m[1][4],
        -x^3*m[1][2]^2*m[2][3]+x^3*m[1][2]*m[1][3]-x^2*m[1][2]*m[2][5]+x^2*m[1][3]-x^2*m[1][5]+x*m[1][2]*m[2][3]^2-2*x*m[1][2]*m[2][4]+x*m[1][3]*m[2][3]+x*m[2][3]+x*m[2][5]-m[2][4]],
       [x*m[1][2]+1]];
comment2:="p^8 automorphisms when x <> 0 mod p";

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#------------------------------------------

#4 7.490
if lname = "7.490" then

auto1:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns1:=[
        [w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5],
         w^2*m[1][3]^2+2*w^2*m[1][5]]];
comment1:="p^8 automorphisms";
#Determinant 1

auto2:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,-1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns2:=[
        [-w*m[1][3]*m[2][3]-w*m[1][4]-w*m[2][5],
         -w^2*m[1][3]^2-2*w^2*m[1][5]]];
comment2:="p^8 automorphisms";
#Determinant -1

auto3:=[
        [-1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns3:=[
        [w*m[1][3]*m[2][3]-w*m[1][4]-w*m[2][5],
         w^2*m[1][3]^2-2*w^2*m[1][5]+2*x*m[2][3]]];
comment3:="p^8 automorphisms if x=0 mod p";
#Determinant -1

auto4:=[
        [-1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,-1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns4:=[
        [-w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5],
         -w^2*m[1][3]^2+2*w^2*m[1][5]+2*x*m[2][3]]];
comment4:="p^8 automorphisms if x=0 mod p";
#Determinant 1

#Automorphism group has order 4*p^8 if x=0 and 2*p^8 if x ne 0

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2),
        rec(auto:=auto3,eqns:=eqns3,comment:=comment3),
        rec(auto:=auto4,eqns:=eqns4,comment:=comment4)];

fi;
#-------------------------------------

#5 7.491
if lname = "7.491" then

auto1:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns1:=[
        [-w*y*m[2][3]^2+2*w*y*m[2][4]+w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5],
         -w^2*y*m[1][3]*m[2][3]+w^2*y*m[1][4]-w^2*y*m[2][5]+w^2*m[1][3]^2+2*w^2*m[1][5]]];
comment1:="p^8 automorphisms, 1-w*y^2 not square mod p";
#Determinant 1

auto2:=[
        [m[1][1],-y*m[1][1]-y,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [w*(y*m[1][1]+y),-m[1][1],m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns2:=[
        [-w^2*y*m[1][1]*m[1][5]-w^2*y*m[1][5]+w*y*m[1][1]*m[2][4]+w*y*m[2][3]^2-w*y*m[2][4]-w*m[1][1]*m[1][4]+w*m[1][1]*m[2][5]-w*m[1][3]*m[2][3]-2*w*m[2][5],
         -w^2*y*m[1][1]*m[1][4]+w^2*y*m[1][1]*m[2][5]+w^2*y*m[1][3]*m[2][3]-2*w^2*y*m[1][4]-w^2*m[1][1]*m[1][5]-w^2*m[1][3]^2-w^2*m[1][5]+w*m[1][1]*m[2][4]-w*m[2][4]],
        [w*y^2*m[1][1]+w*y^2-m[1][1]+1]];
comment2:="p^8 automorphisms, 1-w*y^2 not square mod p";
#Determinant -1

#Automorphism group has order 2*p^8

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#--------------------------------------

#6 7.492
if lname = "7.492" then

auto1:=[
        [1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns1:=[
        [w*x*m[2][3]^2-2*w*x*m[2][4]+w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5],
         w^2*x*m[1][3]*m[2][3]-w^2*x*m[1][4]+w^2*x*m[2][5]+w^2*m[1][3]^2+2*w^2*m[1][5]]];
comment1:="p^8 automorphisms, 1-w*x^2 not square mod p";
#Determinant 1

auto2:=[
        [-1,0,m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,-1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns2:=[
        [-w*x*m[2][3]^2-2*w*x*m[2][4]-w*m[1][3]*m[2][3]-w*m[1][4]+w*m[2][5],
         -w^2*x*m[1][3]*m[2][3]-w^2*x*m[1][4]+w^2*x*m[2][5]-w^2*m[1][3]^2+2*w^2*m[1][5]]];
comment2:="p^8 automorphisms, 1-w*x^2 not square mod p";
#Determinant 1

auto3:=[
        [w*x*m[1][2]-1,m[1][2],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-w*m[1][2],-w*x*m[1][2]+1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns3:=[
        [w^2*x*m[1][2]*m[1][4]-w^2*x*m[1][2]*m[2][5]-w^2*m[1][2]*m[1][5]+w*x*m[2][3]^2+2*w*x*m[2][4]+w*m[1][2]*m[2][4]+w*m[1][3]*m[2][3]-w*m[1][4]-w*m[2][5],
         w^3*x*m[1][2]*m[1][5]-w^2*x*m[1][2]*m[2][4]+w^2*x*m[1][3]*m[2][3]+w^2*x*m[1][4]+w^2*x*m[2][5]-w^2*m[1][2]*m[1][4]+w^2*m[1][2]*m[2][5]+w^2*m[1][3]^2-2*w^2*m[1][5]],
        [w*x^2*m[1][2]-2*x-m[1][2]]];
comment3:="p^8 automorphisms, 1-w*x^2 not square mod p";
#Determinant -1

auto4:=[
        [w*x*m[1][2]+1,m[1][2],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-w*m[1][2],-w*x*m[1][2]-1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]]];
eqns4:=[
        [-w^2*x*m[1][2]*m[1][4]+w^2*x*m[1][2]*m[2][5]+w^2*m[1][2]*m[1][5]-w*x*m[2][3]^2+2*w*x*m[2][4]-w*m[1][2]*m[2][4]-w*m[1][3]*m[2][3]-w*m[1][4]-w*m[2][5],
         -w^3*x*m[1][2]*m[1][5]+w^2*x*m[1][2]*m[2][4]-w^2*x*m[1][3]*m[2][3]+w^2*x*m[1][4]+w^2*x*m[2][5]+w^2*m[1][2]*m[1][4]-w^2*m[1][2]*m[2][5]-w^2*m[1][3]^2-2*w^2*m[1][5]],
        [w*x^2*m[1][2]+2*x-m[1][2]]];
comment4:="p^8 automorphisms, 1-w*x^2 not square mod p";
#Determinant -1

#Automorphism group has order 4*p^8

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2),
        rec(auto:=auto3,eqns:=eqns3,comment:=comment3),
        rec(auto:=auto4,eqns:=eqns4,comment:=comment4)];

fi;
#--------------------------------------

#7 7.1385
if lname = "7.1385" then

auto1:=[
        [1,0,0,m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,1,0,m[2][4],m[2][5],m[1][5]+2*m[1][6],m[2][7]],
        [0,0,1,m[3][4],m[2][4]+m[3][4],-m[1][4]-m[3][4],m[3][7]]];
eqns1:=[];
comment1:="p^9 automorphisms";
#Determinant 1

auto2:=[
        [m[2][3]-1,-m[2][3],m[2][3]-1,m[1][4],m[1][5],m[1][6],m[1][7]],
        [m[2][3]-1,-m[2][3],m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [-1,1,0,m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns2:=[
        [2*m[2][3]-1,
         -m[1][4]*m[2][3]-m[1][4]+m[1][5]*m[2][3]+m[1][6]*m[2][3]-m[1][6]+m[2][3]*m[2][4]-m[2][3]*m[2][5]-m[2][3]*m[2][6]+m[2][6],
         m[2][3]*m[3][4]-m[2][3]*m[3][5]-m[2][3]*m[3][6]-m[2][5]-m[2][6]-m[3][6],
         m[1][4]*m[2][3]+m[1][4]-m[1][5]*m[2][3]-m[1][5]-m[1][6]*m[2][3]-m[2][3]*m[2][4]+m[2][3]*m[2][5]+m[2][3]*m[2][6]+m[2][3]*m[3][4]
           -m[2][3]*m[3][5]-m[2][3]*m[3][6]-m[2][6]-m[3][4]-m[3][6]]];
comment2:="p^9 automorphisms";
#Determinant 1

#Automorphsim group has order 2*p^9

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#-------------------------------------

#8 7.1386
if lname = "7.1386" then

auto1:=[
        [1+m[2][3],-m[2][3],m[2][3],m[1][4],m[1][5],-m[1][5]+x*m[2][3]+m[2][5]+m[2][6],m[1][7]],
        [m[2][3],1-m[2][3],m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [0,0,1,x*m[2][3],m[3][5],m[3][6],m[3][7]]];
eqns1:=[];
comment1:="p^11 automorphisms";
#Determinant 1

auto2:=[
        [m[2][3]-1,-m[2][3],m[2][3]-1,x*m[2][3]-x+m[2][4],m[1][5],m[1][6],m[1][7]],
        [m[2][3]-1,-m[2][3],m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [-1,1,0,m[3][4],m[3][5],x*m[2][3]-m[3][5],m[3][7]]];
eqns2:=[];
comment2:="p^11 automorphisms";
#Determinant 1

#Automorphism group has order 2*p^11

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#-------------------------------------------

#9 7.1388
if lname = "7.1388" then

auto1:=[
       [m[1][2]-2*m[1][3]+1,m[1][2],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
       [m[1][3],-2*m[1][2]+m[1][3]+1,-m[1][2],m[2][4],m[2][5],m[2][6],m[2][7]],
       [m[1][2]-3*m[1][3],3*m[1][2]-m[1][3],m[1][2]+m[1][3]+1,m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns1:=[
       [-m[1][2]*m[1][4]+m[1][2]*m[1][5]+m[1][2]*m[1][6]-m[1][2]*m[2][4]+m[1][2]*m[2][5]+m[1][2]*m[2][6]+m[1][3]*m[1][4]-m[1][3]*m[1][5]-m[1][3]*m[1][6]
          +m[1][3]*m[2][4]-m[1][3]*m[2][5]-m[1][3]*m[2][6]-m[1][5]+m[1][6]-m[2][5]+m[2][6],
        2*m[1][2]*m[2][4]-2*m[1][2]*m[2][5]-2*m[1][2]*m[2][6]+m[1][2]*m[3][4]-m[1][2]*m[3][5]-m[1][2]*m[3][6]-2*m[1][3]*m[2][4]+2*m[1][3]*m[2][5]
          +2*m[1][3]*m[2][6]-m[1][3]*m[3][4]+m[1][3]*m[3][5]+m[1][3]*m[3][6]-2*m[2][4]+4*m[2][5]-m[3][4]+2*m[3][5],
        -m[1][2]*m[1][4]+m[1][2]*m[1][5]+m[1][2]*m[1][6]-3*m[1][2]*m[2][4]+3*m[1][2]*m[2][5]+3*m[1][2]*m[2][6]-m[1][2]*m[3][4]+m[1][2]*m[3][5]
          +m[1][2]*m[3][6]+m[1][3]*m[1][4]-m[1][3]*m[1][5]-m[1][3]*m[1][6]+3*m[1][3]*m[2][4]-3*m[1][3]*m[2][5]-3*m[1][3]*m[2][6]+m[1][3]*m[3][4]
          -m[1][3]*m[3][5]-m[1][3]*m[3][6]+2*m[1][4]-3*m[1][5]-m[1][6]-3*m[2][5]+3*m[2][6]-m[3][4]+2*m[3][6]],
       [-m[1][2]^2+2*m[1][2]*m[1][3]-m[1][2]-m[1][3]^2-m[1][3]]];
comment1:="p^10 automorphisms";    
#The last equation has p solutions over GF(p) for all primes p
#Any solution to this equation mod p can be lifted to a solution mod p^2
#Determinant 1

#Automorphism group has order p^10

return rec(auto:=auto1,eqns:=eqns1,comment:=comment1);

fi;
#-------------------------------------

#10 7.1389
if lname = "7.1389" then

auto1:=[
        [x*m[2][3]-2*m[2][3]+1,-m[2][3],-m[2][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-m[2][3],x*m[2][3]-2*m[2][3]+1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [x*m[2][3]-m[2][3],-x*m[2][3]+m[2][3],-2*m[2][3]+1,m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns1:=[
        [x*m[1][4]*m[2][3]-x*m[1][5]*m[2][3]-x*m[1][6]*m[2][3]+x*m[1][6]+y*m[1][4]*m[2][3]-y*m[1][5]*m[2][3]-y*m[1][6]*m[2][3]
           +y*m[1][6]+y*m[2][3]*m[2][4]-y*m[2][3]*m[2][5]-y*m[2][3]*m[2][6]+y*m[2][5]-2*m[1][4]*m[2][3]+2*m[1][5]*m[2][3]
           -m[1][5]+2*m[1][6]*m[2][3]-m[1][6]+m[2][3]*m[2][4]-m[2][3]*m[2][5]-m[2][3]*m[2][6]+m[2][6],
         -x*y*m[2][5]-x*m[2][3]*m[2][4]+x*m[2][3]*m[2][5]+x*m[2][3]*m[2][6]-x*m[2][3]*m[3][4]+x*m[2][3]*m[3][5]+x*m[2][3]*m[3][6]
           -2*y*m[2][3]*m[2][4]+2*y*m[2][3]*m[2][5]+2*y*m[2][3]*m[2][6]-y*m[2][3]*m[3][4]+y*m[2][3]*m[3][5]+y*m[2][3]*m[3][6]
           +y*m[2][4]-y*m[2][5]-y*m[3][5]+m[2][3]*m[2][4]-m[2][3]*m[2][5]-m[2][3]*m[2][6]+2*m[2][3]*m[3][4]-2*m[2][3]*m[3][5]
           -2*m[2][3]*m[3][6]-m[2][4]-m[3][4]+m[3][5],
         x^2*m[1][4]*m[2][3]-x^2*m[1][5]*m[2][3]-x^2*m[1][6]*m[2][3]+x*y*m[1][4]*m[2][3]-x*y*m[1][5]*m[2][3]-x*y*m[1][6]*m[2][3]
           +x*y*m[2][3]*m[2][4]-x*y*m[2][3]*m[2][5]-x*y*m[2][3]*m[2][6]+x*y*m[2][5]-3*x*m[1][4]*m[2][3]+x*m[1][4]+3*x*m[1][5]*m[2][3]
           -x*m[1][5]+3*x*m[1][6]*m[2][3]-x*m[1][6]+x*m[2][3]*m[2][4]-x*m[2][3]*m[2][5]-x*m[2][3]*m[2][6]+x*m[2][6]+x*m[3][6]
           -2*y*m[1][4]*m[2][3]+y*m[1][4]+2*y*m[1][5]*m[2][3]+2*y*m[1][6]*m[2][3]-y*m[1][6]+y*m[2][3]*m[3][4]-y*m[2][3]*m[3][5]
           -y*m[2][3]*m[3][6]+y*m[3][6]+m[1][4]*m[2][3]-m[1][5]*m[2][3]-m[1][6]*m[2][3]+m[1][6]+m[2][3]*m[3][4]-m[2][3]*m[3][5]
           -m[2][3]*m[3][6]-m[3][4]],
        [x*m[2][3]^2-3*m[2][3]^2+2*m[2][3]]];
comment1:="x <> -1,3 mod p. 2*p^9 automorphisms when x*y+y^2-y+1 <> 0 mod p";
#Determinant 1

#The matrix above is a special case of the matrix below (taking x*m[2][3]-m[2][3]+m[3][2]=0),
#but the matrix above works when x*y+y^2-y+1 ne 0.

auto2:=[
        [x^2*m[2][3]-x*m[2][3]+x*m[3][2]-m[2][3]-m[3][2]+1,-m[2][3],-x*m[2][3]-m[3][2],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-x*m[2][3]-m[3][2],-m[2][3]-m[3][2]+1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [x^2*m[2][3]+x*m[3][2]-m[2][3],m[3][2],-x*m[2][3]-m[2][3]-m[3][2]+1,m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns2:=[
        [x*y*m[2][3]*m[2][4]-x*y*m[2][3]*m[2][5]-x*y*m[2][3]*m[2][6]+x*m[1][6]+y*m[1][4]*m[2][3]-y*m[1][5]*m[2][3]
           -y*m[1][6]*m[2][3]+y*m[1][6]+y*m[2][4]*m[3][2]-y*m[2][5]*m[3][2]+y*m[2][5]-y*m[2][6]*m[3][2]-m[1][4]*m[2][3]
           -m[1][4]*m[3][2]+m[1][5]*m[2][3]+m[1][5]*m[3][2]-m[1][5]+m[1][6]*m[2][3]+m[1][6]*m[3][2]-m[1][6]+m[2][3]*m[2][4]
           -m[2][3]*m[2][5]-m[2][3]*m[2][6]+m[2][6],
         x*y*m[2][3]*m[2][4]-x*y*m[2][3]*m[2][5]-x*y*m[2][3]*m[2][6]+x*y*m[2][5]+y*m[2][3]*m[2][4]-y*m[2][3]*m[2][5]
           -y*m[2][3]*m[2][6]+y*m[2][3]*m[3][4]-y*m[2][3]*m[3][5]-y*m[2][3]*m[3][6]+y*m[2][4]*m[3][2]-y*m[2][4]
           -y*m[2][5]*m[3][2]+y*m[2][5]-y*m[2][6]*m[3][2]+y*m[3][5]-m[2][3]*m[3][4]+m[2][3]*m[3][5]+m[2][3]*m[3][6]
           -m[2][4]*m[3][2]+m[2][4]+m[2][5]*m[3][2]+m[2][6]*m[3][2]-m[3][2]*m[3][4]+m[3][2]*m[3][5]+m[3][2]*m[3][6]+m[3][4]-m[3][5],
         x^2*y*m[2][3]*m[2][4]-x^2*y*m[2][3]*m[2][5]-x^2*y*m[2][3]*m[2][6]+x*y*m[2][3]*m[3][4]-x*y*m[2][3]*m[3][5]
           -x*y*m[2][3]*m[3][6]+x*y*m[2][4]*m[3][2]-x*y*m[2][5]*m[3][2]+x*y*m[2][5]-x*y*m[2][6]*m[3][2]
           -x*m[1][4]*m[2][3]-x*m[1][4]*m[3][2]+x*m[1][4]+x*m[1][5]*m[2][3]+x*m[1][5]*m[3][2]-x*m[1][5]+x*m[1][6]*m[2][3]
           +x*m[1][6]*m[3][2]-x*m[1][6]+x*m[2][3]*m[2][4]-x*m[2][3]*m[2][5]-x*m[2][3]*m[2][6]+x*m[2][6]+x*m[3][6]
           -y*m[1][4]*m[2][3]-y*m[1][4]*m[3][2]+y*m[1][4]+y*m[1][5]*m[2][3]+y*m[1][5]*m[3][2]+y*m[1][6]*m[2][3]+y*m[1][6]*m[3][2]
           -y*m[1][6]+y*m[3][2]*m[3][4]-y*m[3][2]*m[3][5]-y*m[3][2]*m[3][6]+y*m[3][6]+m[1][4]*m[3][2]-m[1][5]*m[3][2]
           -m[1][6]*m[3][2]+m[1][6]+m[2][3]*m[3][4]-m[2][3]*m[3][5]-m[2][3]*m[3][6]-m[3][4]],
        [x*m[2][3]^2+x*m[2][3]*m[3][2]-x*m[2][3]+m[2][3]^2+m[2][3]*m[3][2]-m[2][3]+m[3][2]^2-m[3][2]]];
comment2:="x <> -1,3 mod p. (p-1)*p^9 automorphisms when x*y+y^2-y+1 = 0 mod p";
#Determinant 1

#Note that since x ne -1,3, if x*y+y^2-y+1=0 then y ne 0,1,-1.
#Now let x*m[2][3]^2+x*m[2][3]*m[3][2]-x*m[2][3]+m[2][3]^2+m[2][3]*m[3][2]-m[2][3]+m[3][2]^2-m[3][2]=0
#and x*y+y^2-y+1=0.
#c:=y*(x*m[2][3]^2+x*m[2][3]*m[3][2]-x*m[2][3]+m[2][3]^2+m[2][3]*m[3][2]-m[2][3]+m[3][2]^2-m[3][2]);
#d:=x*y+y^2-y+1;
#u,v:=Quotrem(c,d);
#a:=-(y-1)^2/2*m[2][3]+y*m[3][2]-y/2;
#b:=(y-1)*(y+1)/2*m[2][3]-y*(y-1)/(2*(y+1));
#v*y eq a^2-b^2-y^3/(y+1)^2;
#For any fixed value of y ne 0,1,-1 the pair (a,b) takes all possible values
#in GF(p)^2 as (m[2][3],m[3][2]) take all possible values.
#So v=0 mod p for exactly p-1 choices of (m[2][3],m[3][2]) in GF(p)^2.
#It is fairly easy to see that you can lift any solution (m[2][3],m[3][2]) of
#x*m[2][3]^2+x*m[2][3]*m[3][2]-x*m[2][3]+m[2][3]^2+m[2][3]*m[3][2]-m[2][3]+m[3][2]^2-m[3][2]=0 mod p
#to a solution (m[2][3]',m[3][2]') mod p^2 (where m[2][3]'=m[2][3] mod p, m[3][2]'= m[3][2] mod p).
#The matrix above gives an actionon L/R(L) with determinant 1.

#Automorphism group has order 2*p^9 or (p-1)*p^9 when x*y+y^2-y+1=0 mod p. (Note x ne-1,3.)

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#--------------------------------------------

#11 7.1391
if lname = "7.1391" then

auto1:=[
        [m[3][3],-m[3][2],-m[3][3]-1,m[1][4],m[1][5],m[1][6],m[1][7]],
        [-m[3][2],-2*m[3][3]-1,m[3][2],m[2][4],m[2][5],m[2][6],m[2][7]],
        [-m[3][3]-1,m[3][2],m[3][3],m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns1:=[
        [-2*m[1][4]*m[3][3]-2*m[1][4]+2*m[1][5]*m[3][2]+2*m[1][6]*m[3][3]+m[2][4]*m[3][2]+2*m[2][5]*m[3][3]+m[2][5]-m[2][6]*m[3][2],
         m[1][4]*m[3][2]+2*m[1][5]*m[3][3]+m[1][5]-m[1][6]*m[3][2]-m[2][4]-m[2][6]+m[3][2]*m[3][4]-m[3][2]*m[3][6]+2*m[3][3]*m[3][5]+m[3][5],
         m[2][4]*m[3][2]+2*m[2][5]*m[3][3]+m[2][5]-m[2][6]*m[3][2]-2*m[3][2]*m[3][5]+2*m[3][3]*m[3][4]-2*m[3][3]*m[3][6]-2*m[3][6]],
        [m[3][2]^2+2*m[3][3]^2+2*m[3][3]]];
comment1:="(p-1)*p^9 automorphism if p=1,3 mod 8, (p+1)*p^9 automorphisms if p=5,7 mod 8";
#Determinant 1

auto2:=[
        [m[3][3],-m[3][2],-m[3][3]+1,m[1][4],m[1][5],m[1][6],m[1][7]],
        [m[3][2],2*m[3][3]-1,-m[3][2],m[2][4],m[2][5],m[2][6],m[2][7]],
        [-m[3][3]+1,m[3][2],m[3][3],m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns2:=[
        [2*m[1][4]*m[3][3]-2*m[1][4]-2*m[1][5]*m[3][2]-2*m[1][6]*m[3][3]+m[2][4]*m[3][2]+2*m[2][5]*m[3][3]-m[2][5]-m[2][6]*m[3][2],
         m[1][4]*m[3][2]+2*m[1][5]*m[3][3]-m[1][5]-m[1][6]*m[3][2]-m[2][4]-m[2][6]+m[3][2]*m[3][4]-m[3][2]*m[3][6]+2*m[3][3]*m[3][5]-m[3][5],
         m[2][4]*m[3][2]+2*m[2][5]*m[3][3]-m[2][5]-m[2][6]*m[3][2]+2*m[3][2]*m[3][5]-2*m[3][3]*m[3][4]+2*m[3][3]*m[3][6]-2*m[3][6]],
        [m[3][2]^2+2*m[3][3]^2-2*m[3][3]]];
comment2:="(p-1)*p^9 automorphism if p=1,3 mod 8, (p+1)*p^9 automorphisms if p=5,7 mod 8";
#Determinant 1

#Automorphism group has order 2*(p-1)*p^9 if p=1 mod 8 or 3 mod 8, 2*(p+1)*p^9 otherwise

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#--------------------------------------------

#12 7.1392
if lname = "7.1392" then

auto1:=[
        [-w*m[1][3]+1,m[2][3],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-w*m[2][3],-2*w*m[1][3]+1,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [w^2*m[1][3],-w*m[2][3],-w*m[1][3]+1,m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns1:=[
        [2*w^2*m[1][3]*m[1][6]-2*w*m[1][3]*m[1][4]-2*w*m[1][3]*m[2][5]+2*w*m[1][5]*m[2][3]-2*w*m[1][6]+w*m[2][3]*m[2][6]-m[2][3]*m[2][4]+m[2][5],
         -2*w^2*m[1][3]*m[1][5]+w^2*m[1][6]*m[2][3]-2*w*m[1][3]*m[3][5]-w*m[1][4]*m[2][3]+w*m[1][5]+w*m[2][3]*m[3][6]-w*m[2][6]-m[2][3]*m[3][4]-m[2][4]+m[3][5],
         -2*w^2*m[1][3]*m[2][5]-2*w^2*m[1][3]*m[3][6]+w^2*m[2][3]*m[2][6]+2*w*m[1][3]*m[3][4]-w*m[2][3]*m[2][4]-2*w*m[2][3]*m[3][5]+w*m[2][5]-2*m[3][4]],
        [2*w*m[1][3]^2-2*m[1][3]+m[2][3]^2]];
comment1:="(p+1)*p^9 automorphism if p=1,3 mod 8, (p-1)*p^9 automorphisms if p=5,7 mod 8";
#Determinant 1

auto2:=[
        [-1-w*m[1][3],-m[2][3],m[1][3],m[1][4],m[1][5],m[1][6],m[1][7]],
        [-w*m[2][3],1+2*w*m[1][3],m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [w^2*m[1][3],w*m[2][3],-1-w*m[1][3],m[3][4],m[3][5],m[3][6],m[3][7]]];
eqns2:=[
        [-2*w^2*m[1][3]*m[1][6]+2*w*m[1][3]*m[1][4]-2*w*m[1][3]*m[2][5]+2*w*m[1][5]*m[2][3]-2*w*m[1][6]-w*m[2][3]*m[2][6]+m[2][3]*m[2][4]-m[2][5],
         -2*w^2*m[1][3]*m[1][5]-w^2*m[1][6]*m[2][3]-2*w*m[1][3]*m[3][5]+w*m[1][4]*m[2][3]-w*m[1][5]-w*m[2][3]*m[3][6]-w*m[2][6]+m[2][3]*m[3][4]-m[2][4]-m[3][5],
         -2*w^2*m[1][3]*m[2][5]+2*w^2*m[1][3]*m[3][6]-w^2*m[2][3]*m[2][6]-2*w*m[1][3]*m[3][4]+w*m[2][3]*m[2][4]-2*w*m[2][3]*m[3][5]-w*m[2][5]-2*m[3][4]],
        [2*w*m[1][3]^2+2*m[1][3]+m[2][3]^2]];
comment2:="(p+1)*p^9 automorphism if p=1,3 mod 8, (p-1)*p^9 automorphisms if p=5,7 mod 8";
#Determinant 1

#Automorphism group has order 2*(p+1)*p^9 if p=1 mod 8 or 3 mod 8, 2*(p-1)*p^9 otherwise

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#--------------------------------------------

#13 7.2943
if lname = "7.2943" then

# determine units
U := RingInvariants(L).units; Add(U, D); 

auto1:=[
        [1,m[1][2],-m[1][2],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,m[2][2],0,m[2][4],m[2][5],m[2][6],m[2][7]],
        [0,0,m[2][2],-m[2][5],m[3][5],-m[2][6],m[3][7]]];
eqns1:=[[D*m[2][2]^2-1,2*m[1][2]*m[2][6]-m[1][6]*m[2][2]-2*m[2][5]]];
comment1:="(p-1)*p^10 automorphisms";
#Determinant m[2][2]^2

auto2:=[
        [-1,m[1][2],-m[1][2],m[1][4],m[1][5],m[1][6],m[1][7]],
        [0,0,m[2][3],m[2][4],m[2][5],m[2][6],m[2][7]],
        [0,m[2][3],0,m[3][4],-m[2][4],-m[2][6],m[3][7]]];
eqns2:=[[D*m[2][3]^2-1,2*m[1][2]*m[2][6]+m[1][6]*m[2][3]+2*m[2][4]]];
comment2:="(p-1)*p^10 automorphisms";
#Determinant m[2][3]^2

#Automorphism group has order 2*(p-1)*p^10

return [rec(auto:=auto1,eqns:=eqns1,comment:=comment1),
        rec(auto:=auto2,eqns:=eqns2,comment:=comment2)];

fi;
#---------------------------------------------
#End of difficult LiePrings

    # extend A from mingenset to basis
    img := ExtendAuto(L, A) * b;

    # evaluate relations
    rel := List([1..k], x -> []);
    for i in [1..d] do
        e := Exponents(p*b[i]);
        r := Exponents(p*img[i]) - Exponents(e*img);
        AddByExpo(k, rel, c, r);
        for j in [1..i-1] do
            e := Exponents(b[i]*b[j]);
            r := Exponents(img[i]*img[j]) - Exponents(e*img);
            AddByExpo(k, rel, c, r);
        od;
    od;

    # add det
    Add(rel[1], Determinant(A{[1..a]}{[1..a]}*D^0)*D - 1);

    # record
    old := Unique(Flat(rel));

    # set up
    tod := rec( var := var, char := k,
                rel := List([1..k], x -> []), 
                eli := []);

    # determine units
    U := RingInvariants(L).units; Add(U, D); 

    # reduce system
    for i in Reversed([1..k]) do
        if Length(rel[i]) > 0 then 
            new := Concatenation(rel[i], Flat(tod.rel));
            new := List( new, x -> EliminateVars(x, tod.eli));
            don := false;
            while not don do

                try := StructuralCopy(new);

                # step 1: eliminate
                new := CheckForElimination(tod, new, i);

                # step 2: use units
                # mrvl change
                #if i = 1 then new := List(new, x -> ReduceByUnits(U, x)); fi;
                new := List(new, x -> ReduceByUnits(U, x));

                # step 3: groebner
                new := CallGroebner(new);

                don := (new=try);
            od;
            tod.rel[i] := new;
        fi;
    od;

# mrvl addition.  This is the key addition to ensure that the substitutions
# in the matrix are also carried out in the relations
    for i in [1..k] do
        if Length(tod.rel[i]) > 0 then 
            new := tod.rel[i];
            new := List( new, x -> EliminateVars(x, tod.eli));
            new := Filtered(new, x -> x <> 0*x);
            tod.rel[i] := new;
        fi;
    od;
# end mrvl addition

    # apply result
    M := 0*A; 
    for e in tod.eli do
        for i in [1..a] do
            for j in [1..d] do
                if A[i][j] = e[1] then 
                    A[i][j] := e[2]; 
                    M[i][j] := e[3];
                fi;
            od;
        od;
    od;

    # a special check
    for i in [2..Length(tod.rel)] do
        for j in [1..Length(tod.rel[i])] do
            r := tod.rel[i][j];
            if r=m[1][1]*m[2][2]-m[1][1] and A[2][2]=A[2][2]^0 then r:=0; fi;
            if r=m[1][1]*m[2][2]-m[2][2] and A[1][1]=A[1][1]^0 then r:=0; fi;
            if r=m[1][1]*m[2][2]-m[1][1] and A[2][2]=A[2][2]^0 then r:=0; fi;
            tod.rel[i][j] := r;
        od;
        tod.rel[i] := Unique(tod.rel[i]);
        tod.rel[i] := Filtered(tod.rel[i], x -> x <> 0*x);
    od;

    wd := List(old, x -> FindWord(tod, x));
    wd := Filtered(wd, x -> x <> false);
    if ForAny(wd, x -> x = fail) then Error("check this"); fi;
    wd := Unique(Flat(wd));
    wd := Filtered(wd, x -> x <> 0*x);
        
    # mrvl change
    #return rec( auto := A, eqns := tod.rel, mods := M, word := wd);
    return rec( auto := A, eqns := tod.rel);
end );

BindGlobal( "SizeByBlocks", function(b, p)
    local e, d, c, i;
    e := 1;
    d := Sum(b);
    c := 0;
    for i in [1..Length(b)] do
        c := c + b[i];
        e := e * SizeGL(b[i], p);
        e := e * p^(b[i]*(d-c));
    od;
    return e;
end );

BindGlobal( "BlockForm", function( M, A )
    local d, b, i, j;

    # set up
    M := MakeInt(M);
    d := Length(M);

    # check first row
    if M[1][1] = 0 then return false; fi;
    b := [1];

    # check other rows
    for i in [2..d] do
        b[i] := First([1..d], x -> M[i][x] <> 0);
        for j in [b[i]..d] do
            if M[i][j] <> A[i][j] then return false; fi;
        od;
        if b[i] < b[i-1] then return false; fi;
    od;

    # determine block sizes
    return List(Collected(b), x -> x[2]);
end );

CHECKCASE := false;
BindGlobal( "SizeOfAutoOnFF", function( L )
    local R, d, n, x, p, z, A, i, j, D, a, e, b, w, g3, g4, h3, g5, g7, g8, g9;

    # determine autos
    R := AutGroupDescription(L);
    d := MinimalGeneratorNumberOfLiePRing(L);
    n := DimensionOfLiePRing(L);
    R.auto := R.auto{[1..d]}{[1..d]};

    # set up
    p := IndeterminateByName("p");
    x := IndeterminateByName("x");
    w := IndeterminateByName("w");
    g3 := IndeterminateByName("(p-1,3)");
    h3 := IndeterminateByName("(p+1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");

    # adjust
    z := p^0;
    R.auto := R.auto*z;

    # generic autos
    A := NullMat(d,n);
    for i in [1..d] do
        for j in [1..n] do
            A[i][j] := Indeterminate(Rationals, Concatenation("A",
                        String(i), String(j)));
        od;
    od;

    # inverse of det
    D := Indeterminate(Rationals, "D");
    a := Determinant(R.auto);

    # check eqns
    if R.eqns = [D*a-1] or R.eqns = [1-D*a] then 
        if CHECKCASE then return 1; fi;

        # the GL case
        if R.auto = A{[1..d]}{[1..d]} then 
            return Product(List([0..d-1], x -> (p^d-p^x))); 
        fi;

        # the block case
        b := BlockForm(R.auto, A);
        if b <> false then return b; SizeByBlocks(b,p); fi;
  
        # dim 2
        if R.auto = [ [ A[1][1], A[1][2] ], [ 0, 1 ] ]*z then 
            return (p-1)*p;

        elif R.auto = [ [ 1, A[1][2] ], [ 0, A[2][2] ] ]*z then
            return (p-1)*p;

        elif R.auto = [ [ D*A[2][2], 0 ], [0, A[2][2]] ]*z then 
            return (p-1)^2;

        elif R.auto = [ [ D^2*A[2][2]^2, 0 ], [0, A[2][2]] ]*z then 
            return (p-1)^2/2;

        elif R.auto = [ [ A[1][1], 0 ], [0, A[1][1]*D] ]*z then 
            return (p-1)^2;

        elif R.auto = [ [ A[1][1], A[1][2] ], [0, A[1][1]*D] ]*z then 
            return (p-1)^2*p;

        elif R.auto = [ [ -1, 0 ], [0, A[2][2]] ]*z then 
            return (p-1);

        elif d = 2 then 
            return fail;
        fi;

        # dim 3
        if R.auto = [ [ A[1][1], A[1][2], A[1][3] ], 
                        [ 0, 1, A[2][3] ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p-1)^2*p^3;

        elif R.auto = [ [ A[1][1], A[1][2], 0 ], 
                        [ 0, 1, 0 ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p-1)^2*p;

        elif R.auto = [ [ A[1][1], A[1][2], A[1][3] ], 
                        [ 0, 1, 0 ], 
                        [ 0, A[3][2], A[3][3] ] ]*z then 
            return (p-1)^2*p^3; 

        elif R.auto = [ [ A[1][1], A[1][2], 0 ], 
                        [ A[2][1], A[2][2], 0 ], 
                        [ 0, 0, A[1][1]*A[2][2]-A[1][2]*A[2][1] ] ]*z then 
            return (p^2-1)*(p^2-p);

        elif R.auto = [ [ 1, 0, A[1][3] ], 
                        [ 0, A[2][2], 0 ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p-1)^2*p;

        elif R.auto = [ [ 1, 0, 0 ], 
                        [ 0, A[2][2], A[2][3] ], 
                        [ 0, 0, A[2][2] ] ]*z then 
            return (p-1)*p;

        elif R.auto = [ [ D*A[3][3]^2, A[1][2], 0 ], 
                        [ 0, A[2][2], 0 ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p-1)^3*p;

        elif R.auto = [ [ D*A[2][2], 0, A[1][3] ], 
                        [ 0, A[2][2], A[2][3] ], 
                        [ 0, 0, 1 ] ]*z then 
            return (p-1)^2*p^2;

        elif R.auto = [ [ D*A[3][3]^2, -D*A[3][3]*A[3][4], 0 ], 
                        [ 0, 1, 0 ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p-1)^2*p;

        elif R.auto = [ [ 1, 0, A[1][3] ], 
                        [ 0, A[2][2], 0 ], 
                        [ 0, 0, 1 ] ]*z then 
            return (p-1)*p;

        elif R.auto = [ [ A[1][1], 0, A[1][3] ], 
                        [ 0, A[1][1]*A[3][3], A[2][3] ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p-1)^2*p^2;

        elif R.auto = [ [ A[1][1], A[1][2], 0 ], 
                        [ A[2][1], A[2][2], 0 ], 
                        [ 0, 0, A[3][3] ] ]*z then 
            return (p^2-1)*(p^2-p)*(p-1);

        elif R.auto = [ [ A[1][1], A[1][2], A[1][3] ], 
                        [ 0, A[2][2], 0 ], 
                        [ 0, 0, A[1][1]*A[2][2] ] ]*z then 
            return (p-1)^2*p^2;

        # 5.23 (compare 5.46)
        elif R.auto = [ [ 1, 0, 0 ], 
                        [ 0, A[2][2], A[3][2]*x ], 
                        [ 0, A[3][2], A[2][2]+A[3][2] ] ]*z then 
            return (p^2-1);

        elif d = 3 then 
            return fail;
        fi;

        # dim 4
        if R.auto = [ [ A[1][1], A[1][2], A[1][3], A[1][4] ], 
                      [ 0, 1, A[2][3], A[2][4] ], 
                      [ 0, 0, A[3][3], A[3][4] ], 
                      [ 0, 0, A[4][3], A[4][4] ] ]*z then 
            return (p^2-1)*(p^2-1)*(p-1)*p^5;
 
        elif R.auto = [ [ A[1][1], A[1][2], 0, A[1][4] ], 
                        [ A[2][1], A[2][2], 0, A[2][4] ], 
                        [ 0, 0, A[1][1]*A[2][2]-A[1][2]*A[2][1], A[3][4] ],
                        [ 0, 0, 0, A[4][4] ] ]*z then 
            return (p-1)*(p^2-1)*(p^2-p)*p^3;

        elif R.auto = [ [ A[3][3]*A[4][4]-A[3][4]*A[4][3], A[1][2], 
                          A[3][2]*A[4][3]-A[3][3]*A[4][2], 
                          A[3][2]*A[4][4]-A[3][4]*A[4][2]], 
                        [ 0, 1, 0, 0 ], 
                        [ 0, A[3][2], A[3][3], A[3][4]],
                        [ 0, A[4][2], A[4][3], A[4][4]]]*z then 
            return (p^3-p)*(p^3-p^2)*p;

        elif d = 4 then 
            return fail;
        fi;

    fi;

    if R.eqns = [a-1] or R.eqns = [1-a] then 
        if CHECKCASE then return 2; fi;

        if R.auto = [ [ A[1][1], A[1][2] ], [ 0, A[2][2] ] ]*z then 
            return (p-1)*p;

        elif R.auto = [ [ A[1][1], 0 ], [ 0, A[2][2] ] ]*z then 
            return (p-1);

        elif R.auto = [ [ A[1][1], A[1][2] ], [ 0, A[1][1] ] ]*z then 
            return 2*p;

        elif R.auto = [ [ A[2][2]^2, 0 ], [ 0, A[2][2] ] ]*z then 
            return g3;

        # 5.46
        elif R.auto = [ [ A[1][1], A[2][1]*x ], 
                        [ A[2][1], A[1][1]+A[2][1] ] ]*z then 

            # A11^2 + A11 A21 - x A21^2 = 1: substituting A11 = a-b
            # and A12 = 2b turns this into a^2 - (1+4x)b^2 = 1 and
            # this has p+1 solutions
            return (p+1);
 
        elif d = 2 then 
            return fail; 
        fi;

        if R.auto = [ [ A[1][1], A[1][2], A[1][3]], 
                      [ 0, 1, A[1][2]*A[3][3] ],
                      [ 0, 0, A[3][3] ]]*z then 
            return (p-1)*p^2;

        elif d = 3 then 
            return fail;
        fi;

    fi;

    if CHECKCASE then return 3; fi;

    if R.eqns = [A[1][1]^2*A[2][2]-1] and 
       R.auto = [[A[1][1], A[1][2]], 
                 [0,       A[2][2]]]*z then 
        return (p-1)*p;

    elif R.eqns = [A[1][1]^4-1] and 
         R.auto = [[A[1][1], A[1][2]], 
                   [0,       A[1][1]^2]]*z then 
        return g4*p;

    # case 5.45
    elif R.eqns = [ -A[2][2]^2*w+A[2][1]^2+w, 
                     D^2-1, 
                     A[1][2]*w-D*A[2][1], 
                    -D*A[2][2]^2+A[1][2]*A[2][1]+D ] and 
         R.auto = [[D*A[2][2], A[1][2]], 
                   [A[2][1], A[2][2]]]*z then 

         # A21 = w/D A12, D = +/- 1, A22^2 - w A21^2 = 1
         # there are (p+1) options for (A22, A21) and two options for D
         return 2*(p+1);

    elif R.eqns = [ (x-1)*(x+1)*A[2][1], 
                    (x-1)*A[2][2]*A[2][1], 
                    (x-1)*(x+1)*A[1][2], 
                    (x-1)*A[2][2]*A[1][2], 
                    (x-1)*A[2][1]*(A[1][2]*A[2][1]-1), 
                    (x-1)*A[1][2]*(A[1][2]*A[2][1]-1), 
                    -A[1][2]*A[2][1]*x + A[1][1]*A[2][2] - 1,
                    (x-1)*A[2][1]*A[1][1], 
                    (x-1)*A[1][2]*A[1][1] ] and 
         R.auto = [[A[1][1], A[1][2]], 
                   [A[2][1], A[2][2]]]*z then 

         # if x = 1 then this is SL(2,p)
         # if x <> 1 then 
              # A22 <> 0 yields A21 = A12 = 0 and A11 A22 = 1, thus (p-1)
              # A22 = 0 yields x A12 A21 = -1 and A11 = 0 and x=-1, thus (p-1)

         return "SL(2,p) if x=1, 2(p-1) if x =-1, (p-1) otherwise";

    elif d = 2 then 
        return fail;
    fi;

    if R.eqns = [D*A[3][3]-1, 
                 A[1][1]*A[2][2]-1] and 
         R.auto = [[A[1][1], A[1][2], A[1][3]],
                   [0,       A[2][2], A[2][3]],
                   [0,       0,       A[3][3]]]*z then 
        return (p-1)^2*p^3;

    elif R.eqns = [A[2][2]^2*A[3][3]-1] and 
         R.auto = [[A[2][2]*A[3][3], -A[2][2]*A[3][4], A[1][3]],
                   [0,                A[2][2],         A[2][3]], 
                   [0,                0,               A[3][3]]] * z then  
        return (p-1)*p^3;

    # 5.20
    elif R.eqns = [ (x-1)*(x+1)*A[3][2], 
                    (x-1)*A[3][3]*A[3][2], 
                    (x-1)*(x+1)*A[2][3],
                    (x-1)*A[3][3]*A[2][3],
                    (x-1)*A[3][2]*A[2][2],
                    (x-1)*A[2][3]*A[2][2],
                    (x-1)*A[3][2]*(D*A[2][3]*A[3][2]-1),
                    (x-1)*A[2][3]*(D*A[2][3]*A[3][2]-1),
                    D*(-A[2][3]*A[3][2]*x+A[2][2]*A[3][3])-1 ] and
         R.auto = [[ D*A[2][3]*A[3][2]*(x-1)+1, 0, 0], 
                   [ 0, A[2][2], A[2][3] ], 
                   [ 0, A[3][2], A[3][3] ] ] * z then 

         # if x = 1, then this is GL(2,p)
         # if x <> 1 then 
            # A33 <> 0 yields A32=A23=0 and D A22 A33 = 1, thus (p-1)^2
            # A33 = 0 yields D A23 A32 = 1 and A22 = 0 and x=-1, thus (p-1)^2 
    
         return "GL(2,p) if x=1, 2(p-1)^2 if x=-1, (p-1)^2 otherwise";

    # 5.22 
    elif R.eqns = [ -A[3][2]^2*w^2+A[2][3]^2, 
                     A[2][2]*A[3][2]*w-A[2][3]*A[3][3], 
                     -A[3][2]*A[3][3]*w+A[2][2]*A[2][3], 
                     A[2][2]^2-A[3][3]^2, 
                     D*(A[3][2]^2*w-A[3][3]^2)+1, 
                     D*A[3][3]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])-A[2][2] ] and
         R.auto = [[ D*(A[2][2]*A[3][3]-A[2][3]*A[3][2]), 0, 0], 
                 [ 0, A[2][2], A[2][3] ],
                 [ 0, A[3][2], A[3][3] ]] * z then 

         # A23 = e w A32 and A33 = e A22 with e = +/-1 and w A32^2 + D = A33^2
         # yields p^2-1 options for (A32, A33) and these determine D
 
         return 2*(p^2-1);

    elif d = 3 then 
        return fail;
   
    fi;

    # 5.6 is a subgroup of a symplectic group with A11 <> 0 and second 
    # row 0,1,0,0
    if LibraryName(L) = "5.6" then 
        return p^3*(p-1)*(p^2-1)*p;
    fi;

    return fail;

end );

BindGlobal( "CheckAGDescription", function(arg)
    local L, p, b, d, a, q, c, k, D, A, i, j, img, rel, e, r, U, 
          new, tod;

    # catch args
    L := arg[1];
    A:=arg[2];

    # set up
    p := PrimeOfLiePRing(L);
    b := BasisOfLiePRing(L);
    d := DimensionOfLiePRing(L);
    a := MinimalGeneratorNumberOfLiePRing(L);
    q := ParametersOfLiePRing(L);
    c := List(b, x -> Length(Factors(LiePOrder(x))));
    k := Maximum(c);

    # variables
    D := Indeterminate(Rationals, "D");

    # extend A from mingenset to basis
    img := ExtendAuto(L, A) * b;

    # evaluate relations
    rel := List([1..k], x -> []);
    for i in [1..d] do
        e := Exponents(p*b[i]);
        r := Exponents(p*img[i]) - Exponents(e*img);
        AddByExpo(k, rel, c, r);
        for j in [1..i-1] do
            e := Exponents(b[i]*b[j]);
            r := Exponents(img[i]*img[j]) - Exponents(e*img);
            AddByExpo(k, rel, c, r);
        od;
    od;

    # set up
#    tod := rec( char := k,
#                rel := List([1..k], x -> []));

    # determine units
#    U := RingInvariants(L).units; Add(U, D); 

    # reduce system
#    for i in Reversed([1..k]) do
#        if Length(rel[i]) > 0 then 
#            new := Concatenation(rel[i], Flat(tod.rel));

                # step 2: use units
#                new := List(new, x -> ReduceByUnits(U, x));

                # step 3: groebner
#                new := CallGroebner(new);

#            tod.rel[i] := new;
#        fi;
#    od;

    return rel;
end );

