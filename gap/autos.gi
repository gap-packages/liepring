ExtendAuto := function(L, mat)
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
end;

AutoOnMult := function( C, mat )
    local dim, sml, new;
    dim := DimensionOfLiePRing(C);
    sml := dim - DimensionOfLiePRing(C!.mult);
    new := ExtendAuto(C,mat);
    return new{[sml+1..dim]}{[sml+1..dim]};
end;

IndetsOfPoly := function( pol )
    local e, ind, i, j, x;
    e := ExtRepPolynomialRatFun(pol);
    ind := [];
    for i in [1,3..Length(e)-1] do
        for j in [1,3..Length(e[i])-1] do
            x := Indeterminate(Rationals, e[i][j]);
            AddSet(ind,x);
        od;
    od;
    return ind;
end;

MyVals := function( pol, vars, vals )
    if IsInt(pol) then return pol; fi;
    if IsBool(pol) then return pol; fi;
    return Value(pol, vars, vals);
end;
   
TryEliminate := function( pol )
    local v, x, q, r, w;
    v := IndetsOfPoly(pol);
    w := IndeterminateByName("w");
    for x in v do
        q := PolynomialReduction( pol, [x], ORDER);
        r := q[1];
        q := MakeInt(q[2][1]);
        if q in [1,-1] then 
            return [x,-r/q];
        fi;
    od;
    return fail;
end;

CHECKCASE := false;
 
FindAutos := function(L)
    local p, b, d, a, s, x, f, A, D, i, j, h, k, img, rel, e, r, w, new, don,
          uni, u, vars, y, z;

    # set up
    p := PrimeOfLiePRing(L);
    b := BasisOfLiePRing(L);
    d := DimensionOfLiePRing(L);
    a := MinimalGeneratorNumberOfLiePRing(L);
    w := IndeterminateByName("w");
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");

    # generic autos
    A := MutableNullMat(a,d);
    for i in [1..a] do
        for j in [1..d] do
            A[i][j] := Indeterminate(Rationals, Concatenation("A",
                        String(i), String(j)));
        od;
    od;
    vars := Flat(A);

    # inverse of det
    D := Indeterminate(Rationals, "D");
    Add(vars, D);

    # create generic auto
    img := ExtendAuto(L, A) * b;

    # evaluate relations
    rel := [];
    for i in [1..d] do
        e := Exponents(p*b[i]);
        r := Exponents(p*img[i]) - Exponents(e*img);
        Add(rel, r);
        for j in [1..i-1] do
            e := Exponents(b[i]*b[j]);
            r := Exponents(img[i]*img[j]) - Exponents(e*img);
            Add(rel, r);
        od;
    od;

    # add det
    Add(rel, Determinant(A{[1..a]}{[1..a]})*D - 1);

    # prune
    rel := Flat(rel);
    rel := Filtered(rel, x-> x <> 0*x);
    rel := List(rel, x -> NumeratorOfRationalFunction(x));

    # determine units
    uni := [w, D];
    if L!.LibraryName in ["5.42", "5.20", "5.46", "5.23"] then 
        Add(uni, x); 
    elif L!.LibraryName in ["7.1798"] then 
        Add(uni, y); Add(uni, z);
    elif L!.LibraryName in ["7.1799"] then 
        Add(uni, z);
    fi;

    # reduce system
    repeat
        new := StructuralCopy(rel);
        for i in [1..Length(new)] do
            for u in uni do
                if IsPolynomial(new[i]/u) then new[i] := new[i]/u; fi;
            od;
        od;
        Print("  -- calling groebner \n");
        new := GroebnerBasis(new, ORDER);
        Print("  ---- done \n");
        don := (new=rel);
        rel := new;
    until don;

    # create result
    for i in [1..Length(rel)] do
        if not IsBool(rel[i]) and not IsInt(rel[i]) then 
            e := TryEliminate(rel[i]);
            if not IsBool(e) then 
                for h in [1..a] do
                    for k in [1..d] do
                        A[h][k] := MyVals(A[h][k], [e[1]], [e[2]]);
                    od;
                od;
                for j in [1..Length(rel)] do
                    rel[j] := MyVals(rel[j], [e[1]], [e[2]]);
                od;
                rel[i] := false;
            fi;
        fi;
    od;
    
    # prune once more
    rel := Filtered(rel, x -> x <> false);
    rel := Filtered(rel, x -> x <> 0*x);
    if Length(rel) > 0 then 
        Print("  -- calling groebner \n");
        rel := GroebnerBasis(rel, ORDER);
        Print("  ---- done \n");
    fi;
            
    return rec( auto := A, eqns := rel, vars := vars);
end;

SizeOfGL := function(n, p)
    return Product(List([0..n-1], x -> (p^n - p^x)));
end;

SizeByBlocks := function(b, p)
    local e, d, c, i;
    e := 1;
    d := Sum(b);
    c := 0;
    for i in [1..Length(b)] do
        c := c + b[i];
        e := e * SizeOfGL(b[i], p);
        e := e * p^(b[i]*(d-c));
    od;
    return e;
end;

BlockForm := function( M, A )
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
end;

SizeOfAutoOnFF := function( L )
    local R, d, n, x, p, z, A, i, j, D, a, e, b, w, g3, g4, h3, g5, g7, g8, g9;

    # determine autos
    R := FindAutos(L);
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
    A := MutableNullMat(d,n);
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

end;


