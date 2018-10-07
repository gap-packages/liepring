
EvaluatePorcPoly := function( f, P )
    local inds, vals, p, g3, g4, g5, g6, g7, g8, g9, h3, s16;

    if not IsPrimeInt(P) or P = 2 then return false; fi;

    # check
    if IsInt(f) then return f; fi;

    # get prime
    p := IndeterminateByName("p");

    # get indets
    g3 := IndeterminateByName("(p-1,3)");
    h3 := IndeterminateByName("(p+1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");
    s16 := IndeterminateByName("(p^2-1,16)");
    inds := [p, g3, h3, g4, g5, g7, g8, g9, s16];

    # get values
    vals := [P, Gcd(P-1,3), Gcd(P+1,3), Gcd(P-1,4), Gcd(P-1,5), 
             Gcd(P-1,7), Gcd(P-1,8), Gcd(P-1,9), Gcd(P^2-1, 16)];

    # return
    return Value( f, inds, vals );
end;

DegreePorcPoly := function( f )
    local p, a, b;
    if IsInt(f) then return 0; fi;
    p := IndeterminateByName("p");
    a := DegreeIndeterminate(NumeratorOfRationalFunction(f),p);
    b := DegreeIndeterminate(DenominatorOfRationalFunction(f),p);
    if b <> 0 then Error("something wrong with degrees "); fi;
    return a;
end;

EliminateDenominator := function( f )
    local g3, g4, g5, g7, g8, a, b;
   
    # check the easy case
    if IsInt(f) or IsPolynomial(f) then return f; fi;

    # translate
    g3 := IndeterminateByName("(p-1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");

    b := DenominatorOfRationalFunction(f);
    a := NumeratorOfRationalFunction(f);
    if b = g3 then return a*(4-b)/3; fi;
    if b = g4 then return a*(6-b)/8; fi;
    if b = g5 then return a*(6-b)/5; fi;
    if b = g7 then return a*(8-b)/7; fi;
    if b = g8 then return a*(b^2/64 - 7*b/32 + 7/8); fi;
 
    Error("unknown denominator");
end;

SimplifyPorcPoly := function( f )
    local g, g3, g4, g5, g7, o, g8, g9, p, h3, s16;

    if IsInt(f) then return f; fi;

    # get the indeterminates
    p := IndeterminateByName("p");
    g3 := IndeterminateByName("(p-1,3)");
    h3 := IndeterminateByName("(p+1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");
    s16 := IndeterminateByName("(p^2-1,16)");

    o := MonomialLexOrdering([g8,g4,g9,g3,g5,g7,h3,s16,p]);
    if DegreeIndeterminate(f, g9) > 1 then 
        f := PolynomialReducedRemainder( f, [g9^2-12*g9+8*g3+3], o );
        if DegreeIndeterminate(f,g9)>1 then Error("wrong reduction g9"); fi;
    fi;
    if DegreeIndeterminate(f, g3) > 1 then 
        f := PolynomialReducedRemainder( f, [g3^2-4*g3+3], o );
        if DegreeIndeterminate(f,g3)>1 then Error("wrong reduction g3"); fi;
    fi;
    if DegreeIndeterminate(f, g5) > 1 then 
        f := PolynomialReducedRemainder( f, [g5^2-6*g5+5], o );
        if DegreeIndeterminate(f,g5)>1 then Error("wrong reduction g5"); fi;
    fi;
    if DegreeIndeterminate(f, g7) > 1 then 
        f := PolynomialReducedRemainder( f, [g7^2-8*g7+7], o );
        if DegreeIndeterminate(f,g7)>1 then Error("wrong reduction g7"); fi;
    fi;
    if DegreeIndeterminate(f, g8) > 1 then 
        f := PolynomialReducedRemainder( f, [g8^2-12*g8+6*g4+8], o );
        if DegreeIndeterminate(f,g8)>1 then Error("wrong reduction g8"); fi;
    fi;
    if DegreeIndeterminate(f, g4) > 1 then 
        f := PolynomialReducedRemainder( f, [g4^2-6*g4+8], o );
        if DegreeIndeterminate(f,g4)>1 then Error("wrong reduction g4"); fi;
    fi;
    if DegreeIndeterminate(f, g8) > 0 and DegreeIndeterminate(f, g4) > 0 then 
        f := PolynomialReducedRemainder( f, [g4*g8-4*g8-2*g4+8], o );
    fi;

    return f;
end;
    
NumberOfLiePRingsInFamilyNC := function( L )
    local S, d, p, w, x, y, z, t, j, k, m, n, r, s, u, v, para, c, 
          g3, h3, g4, g5, g7, g8, g9, s16, a;

    # get some data
    S := SCTable(Zero(L));
    d := DimensionOfLiePRing(L);

    # the trivial case
    if not IsBound(S.param) then return 1; fi;

    # get number of parameters of L
    w := IndeterminateByName("w");
    para := Difference( S.param, [w] );

    # get the indeterminates
    p := IndeterminateByName("p");
    w := IndeterminateByName("w");
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    t := IndeterminateByName("t");
    j := IndeterminateByName("j");
    k := IndeterminateByName("k");
    m := IndeterminateByName("m");
    n := IndeterminateByName("n");
    r := IndeterminateByName("r");
    s := IndeterminateByName("s");
    u := IndeterminateByName("u");
    v := IndeterminateByName("v");

    # get gcd's
    g3 := IndeterminateByName("(p-1,3)");
    h3 := IndeterminateByName("(p+1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");
    s16 := IndeterminateByName("(p^2-1,16)");

    # special case
    a := IndeterminateByName("a");
    ## a equals the number of non-zero a in GF(p) such that 
    ## 2*w^2*a^3+5*w*a^2+2*a is a square

    # get comment
    c := LibraryConditions(L)[1];

    # 0 parameters
    if Length(para) = 0 then return 1; fi;

    # 1 parameter
    if Length(para) = 1 then 

        # para = [x]
        if (Length(c) = 0) or (c = "all x") then 
            return p;
        elif c = "x ne -1" then 
            return p-1;
        elif c = "x ne -1,-1/2" then 
            return p-2;
        elif c = "x ne -2" then 
            return p-1;
        elif c = "x ne 0" then 
            return p-1;
        elif c = "x ne 0,-1" then 
            return p-2;
        elif c = "x ne 0,-1,2,1/2" then 
            return p-4;
        elif c = "x ne 0,-1/4" then 
            return p-2;
        elif c = "x ne 0,-w,2w,w/2" then 
            return p-4;
        elif c = "x ne 0,1" then 
            return p-2;
        elif c = "x ne 1" then 
            return p-1;
        elif c = "x=0,1,w" then 
            return 3;
        elif c = "x=1/2" then 
            return 1;
        elif c = "x=1/4" then 
            return 1;
        elif c = "x=0,1,w,w^2,w^3" then 
            return 5;
        elif c = "x=w,w^2,w^3,w^4" or c = "x=w^i, i=1,2,3,4" then 
            return 4;
        elif c = "x=w,w^2,w^3,w^4,w^5,w^6" or c = "x=w^i, i=1,2,3,4,5,6" then 
            return 6;
        elif c = "x=w^2,w^3,w^4,w^5" or c = "x=w^i, i=2,3,4,5" then 
            return 4;
        elif c = "x=w^4,w^5,w^6,w^7" then 
            return 4;
        elif c = "x=w/2" then 
            return 1;
        elif c = "x=w^i, 0 le i lt 2gcd[p-1,3]" then 
            return 2*g3;
        elif c = "x=w^i, 0 le i lt gcd[p-1,8]" then 
            return g8;
        elif c = "all x, x~-x" or c = "all x,x~-x" or c = "x~-x" then 
            return (p+1)/2;
        elif c = "x ne 0, x~-x" then 
            return (p-1)/2;
        elif c = "x ne 0, x~x^-1" or c = "x ne 0,x~x^-1" then 
            return (p+1)/2;
        elif c = "x ne 0,1,-1, x~x^-1" then 
            return (p-3)/2;
        elif c = "all x, x~1-x" then 
            return (p+1)/2;
        elif c = "x ne 0,-2, x~-x-2" then 
            return (p-1)/2;
        elif c = "all x, x~w-x" then 
            return (p+1)/2;
        elif c = "x ne 0,-w, x~-w-x" then 
            return (p-1)/2;
        elif c = "1+4x not a square mod p" or
             c = "all x, 1+4x not a square mod p" then 
            return (p-1)/2;
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p" then 
            return (p-1)/g5;
        elif c = "x ne 0,  x~x' if x^4=x'^4 mod p" or 
             c = "x ne 0, x~x' if x^4=x'^4 mod p" then 
            return (p-1)/g4;
        elif c = "x ne 0, x~x' if x^3=x'^3  mod  p" or
             c = "x ne 0, x~x' if x^3=x'^3 mod p" or
             c = "x ne 0,x~x' if x^3=x'^3 mod p" then 
            return (p-1)/g3;
        elif c = "x~x' if x^4=x'^4 mod p" or
             c = "all x, x~x' if x^4=x'^4 mod p" then 
            return 1+(p-1)/g4;
        elif c = "all x, x~x' if x^3=x'^3 mod p" then 
            return 1+(p-1)/g3;
        elif c = "x ne 0, x~x' if x^6=x'^6 mod p" then 
            return (p-1)/(2*g3);
        elif c = "x ne 0, x~x' if x^7=x'^7 mod p" then 
            return (p-1)/g7;
        elif c = "x FIXED with x^2-w not a square mod p" then 
            return 1;
        elif c = "x ne 0, equivalence classes {x,-x,1/x,-1/x}" then 
            #return (p-1+g4)/2;
            return (p-1+g4)/4;
    
        # para = [y]
        elif c = "all y" then 
            return p;
        elif c = "y=0,1" then 
            return 2;
        elif c = "y=0,1,w" then 
            return 3;
        elif c = "y=w,w^2" then 
            return 2;
        elif c = "y=w^3,w^4,w^5,w^6,w^7,w^8" then 
            return 6;
        elif c = "y=w^i, i=1,2,3,4" then 
            return 4;
        elif c = "y=w^i, i=2,3,4,5" then
            return 4;
        elif c = "y=w,w^2" then
            return 2;
        elif c = "y=w^3,w^4,w^5,w^6,w^7,w^8" then 
            return 6;
        elif c = "y FIXED such that 1-wy^2 is not a square mod p" then
            return 1;
        elif c = "all y, y~-y" then 
            return (p+1)/2;
        elif c = "all y, y~-y-2" then 
            return (p+1)/2;
        elif c = "all y, y~y' if y^4=y'^4 mod p" then 
            return 1+(p-1)/g4;
        elif c = "y ne 0, y~y' if y^3=y'^3 mod p" then 
            #return (p-1)/3;
            return (p-1)/g3;
        elif c = "y ne 0, y~y' if y^4=y'^4 mod p" then 
            return (p-1)/g4;

        # para = [z]
        elif c = "all z" then 
            return p;
        elif c = "z ne 0" then 
            return p-1;
        elif c = "z ne 0, z~-z" then 
            return (p-1)/2; 
        elif c = "z=0,1" then 
            return 2;
        elif c = "z=1/2" then 
            return 1;
        elif c = "z=w^i, i=1,2,3,4" then 
            return 4;

        # para = [t]
        elif c = "all t" then 
            return p;
        elif c = "t=0,1,w" then 
            return 3;
        fi;
  
    # 2 parameter
    elif Length(para) = 2 then 

        # sort parameters correctly
        if x in para then # [x,y], [x,z], [x,t]
            para := [x, Difference(para,[x])[1]]; 
        elif y in para then # [y,z], [y,t]
            para := [y, Difference(para,[y])[1]]; 
        elif z in para then # [z, t]
            para := [z, Difference(para,[z])[1]]; 
        fi;

        # para = [x,y];
        if Length(c) = 0 or c = "all x,y" then 
            return p^2;
        elif c = "all x, y ne 0" then 
            return p^2-p;
        elif c = "x ne 0, all y" then 
            return p^2-p;
        elif c = "x ne 0, all y, y~-y" then 
            return (p^2-1)/2;
        elif c = "all x,y, x~-x" or c = "all x, x~-x, all y" then 
            return p*(p+1)/2;
        elif c = "all x, x~-x, y ne 0" then 
            return (p^2-1)/2;
        elif c = "x ne 0, all y, x~-x" or c = "x ne 0, x~-x, all y" then 
            return p*(p-1)/2;
        elif c = "x ne 0, all y, x~-x, y~-y" or 
             c = "x ne 0, x~-x, all y, y~-y" then 
            return (p^2-1)/4;
        elif c = "all x,y, y~-y" or c = "all x, all y, y~-y" then 
            return p*(p+1)/2;
        elif c = "all x, y ne 0, y~-y" then 
            return p*(p-1)/2;
        elif c = "y=1,w, all x, x~-x" then 
            return p+1;
        elif c = "x ne 0, x~-x, y=1,w" then 
            return p-1;
        elif c = "all x, x~-x, y=1/2" then 
            return (p+1)/2;
        elif c = "x ne 0, x~-x, y=w^-1" then 
            return (p-1)/2;
        elif c = "x ne 0, y ne 0,1,-1, y~y^-1" then 
            return (p-1)*(p-3)/2;
        elif c = "x ne 0, y ne 0 ,y~y^-1" then 
            return (p-1)*(p+1)/2;
        elif c = "all x, y=w^2,w^3,w^4,w^5" then 
            return 4*p;
        elif c = "all x, y=0,1" then 
            return 2*p;
        elif c = "all x, y ne 1/2, y~1-y" then 
            return p*(p-1)/2;
        elif c = "x ne -1, y=[1+x]^-1 mod p" then 
            return p-1;
        elif c = "x ne -2, all y" then 
            return p*(p-1);
        elif c = "x ne -2w, all y" then 
            return p*(p-1);
        elif c = "y=1,w,w^2, all x, x1~x2 if x1^3=x2^3 mod p" or
             c = "y=1,w,w^2, all x, x~x' if x^3=x'^3 mod p" then 
            return 3*(1+(p-1)/g3);
        elif c = "x ne 0, all y, x~x' if x^3=x'^3 mod p" then 
            return p*(p-1)/g3;
        elif c = "all x,y, x~x' if x^3=x'^3 mod p" then 
            return p*(1+(p-1)/g3);
        elif c = "y=1,w, all x, x~x' if x^8=x'^8 mod p" then 
            return 2*(1+(p-1)/g8);
        elif c = "y=1,w,w^2, all x, x1~x2 if x1^3=x2^3 mod p" then 
            return 3*(1+(p-1)/g3);
        elif c = "y=1,w,w^2,w^3,w^4,w^5, all x, x~x' if x^6=x'^6 mod p" then 
            return 6*(1+(p-1)/(2*g3));
        elif c = "y=w^2,w^3, all x, x~x' if x^8=x'^8 mod p" then 
            return 2*(1+(p-1)/g8);
        elif c = "y=w^4,w^5,w^6,w^7, all x, x~x' if x^8=x'^8 mod p" then 
            return 4*(1+(p-1)/g8);
        elif c = "x ne 0, all y, x~x' if x^3=x'^3 mod p" or 
             c = "x ne 0, x~x' if x^3=x'^3 mod p, all y" then 
            return p*(p-1)/g3;
        elif c = "x ne 0, x~x' if x^6=x'^6 mod p, all y" then 
            return p*(p-1)/(2*g3);
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p, all y" then 
            return p*(p-1)/g5;
        elif c = "y=w^i, i=1,2,3,4, x ne 0, x~x' if x^5=x'^5 mod p" then 
            return 4*(p-1)/g5;
        elif c = "y=w^i, i=1,2,3,4,5,6, x ne 0, x~x' if x^7=x'^7 mod p" then 
            return 6*(p-1)/g7;
        elif c = "y=w^i, i=2,3,4,5, x ne 0, x~x' if x^3=x'^3 mod p" then
            return 4*(p-1)/g3;
        elif c = "x ne 1, all y" then 
            return p*(p-1);
        elif c = "x=0,-1, y=0,1" then 
            return 4;
        elif c = "x=0,1, all y" then 
            return 2*p;
        elif c = "x=0,1, y=0,1" then 
            return 4;
        elif c = "x ne 1-w, x~-x+2[1-w], all y, y~-y" then 
            return (p^2-1)/4;
        elif c = "x ne 1-w^2, x~-x+2[1-w^2], all y, y~-y" then 
            return (p^2-1)/4;
        elif c = "x ne 1-w^3, x~-x+2[1-w^3], all y, y~-y" then 
            return (p^2-1)/4;
        elif c = "all x, y ne 0, y~y' if y^4=y'^4 mod p" then 
            return p*(p-1)/g4;
        elif c = "all x, y ne 0, y~y' if y^5=y'^5 mod p" or
             c = "all x, y ne 0, y~y' if y^5=y'^5 mod p" then 
            return p*(p-1)/g5;
        elif c = "all x,y, x~-x, y~y' if y^3 eq y'^3  mod  p" or
             c = "all x,y, x~-x, y~y' if y^3=y'^3 mod p" then 
            return (p+1)*(1+(p-1)/g3)/2;
        elif c = "all x, y ne 0, y~y' if y^3=y'^3 mod p" or 
             c = "y ne 0, y~y' if y^3=y'^3 mod p, all x" then 
            return p*(p-1)/g3;
        elif c = "all x, y ne 0, y~y' if y^4=y'^4 mod p" then 
            return p*(p-1)/g4;
        elif c = "all x,y, y~y' if y^4=y'^4 mod p" then 
            return p*(1+(p-1)/g4);
        elif c = "x=w^i, 0 le i lt 2gcd[p-1,3], y ne 0, y~y' if y^6=y'^6 mod p"
        then 
            return p-1;
        elif c = "x=w^i, 0 le i lt gcd[p-1,8], y ne 0, y~y' if y^8=y'^8 mod p"
        then 
            return p-1;
        elif c = "all x,y, [x,y]~[y,x]" then 
            return p*(p+1)/2;
        elif c = "x,y ne 0, x ne y, [x,y]~[y,x]" then 
            return (p-1)*(p-2)/2;
        elif c = "all x,y such that xy ne 1, [x,y]~[y,x]" then 
            return (p^2-1)/2;
        elif c = "all x,y, [x,y]~[x',y'] if y^2-wx^2=y'^2-wx'^2 mod p" then 
            return p;
        elif c ="x ne 0, x~-x, y FIXED such that 1-wy^2 is not a square mod p"
        then 
            return (p-1)/2;
        elif c = "x ne 0, all y, 4x+y^2 not a square mod p" then 
            return p*(p-1)/2;
        elif c = "If 2 is a square mod p and y^2=2, x ne 0, x~-x, y~-y" then 
            return (s16-8)*(p-1)/16;
        elif c = "If 2 is a square mod p and y^2=2w^2, x ne 0, x~-x, y~-y" then
            return (s16-8)*(p-1)/16;
        elif c = "If 2 is not a square mod p and wy^2=2, x ne 0, x~-x, y~-y" 
        then
            return (16-s16)*(p-1)/16;
        elif c = "If 2 is not a square mod p and y^2=2w^3, x ne 0, x~-x, y~-y"
        then 
            return (16-s16)*(p-1)/16;
        elif d = 6 and c = "See note6.62" then 
            return p;
        elif d = 7 and c = "See Notes6.178" then 
            return (p+1)/2;
        elif c = "See note1dec5.1" then 
            return p;

        # para = [x,z]
        elif c = "all x,z" then 
            return p^2;
        elif c = "x ne 0, all z" then 
            return p*(p-1);
        elif c = "all x, z=1,w" then
            return 2*p;
        elif c = "x,z ne 0, [x,z]~[z,x]" then 
            return p*(p-1)/2;
        elif c = "x ne 0, x~-x, z=w^-1" then 
            return (p-1)/2;
        elif c = "x ne 0, x~x' if x^4=x'^4 mod p, z=w^-1" then 
            return (p-1)/g4;
        elif c = "x ne 0, x~x' if x^6=x'^6 mod p, all z" then 
            return p*(p-1)/(2*g3);
        elif c = "z=w^i, i=1,2,3,4, x ne 0, x~x' if x^5=x'^5 mod p" then 
            return 4*(p-1)/g5;
        elif c = "x ne -1,3, See Notes6.114" then 
            return 3*p-7;

        # para = [x,t]
        elif c = "all x, t=0,1" then 
            return 2*p;
        elif c = "t=1,-1, x ne 0" then 
            return 2*(p-1);
        elif c = "t=1,-1, x ne 0, x~-x" then 
            return p-1;
        elif c = "all x,t, [x,t]~[t,x]" then 
            return p*(p+1)/2;

        # para = [y,z]
        elif c = "all y, z ne 0, y~-y, z~-z" then 
            return (p^2-1)/4;
        elif c = "all y, z=0,1,w" then 
            return 3*p;
        elif c = "all y,z" then 
            return p^2;
        elif c = "y,z ne 0, y~-y" then 
            return (1/2)*(p-1)^2;
        elif c = "y,z ne 0, y~y' if y^4=y'^4 mod p" then 
            return (p-1)*(p-1)/g4;
        elif c = "y=0,1, all z" then 
            return 2*p;
        elif c = "y=0,1, z=0,1,w" then 
            return 6;
        elif c = "y=1,w, all z" then 
            return 2*p;

        # para = [y,t]
        elif c = "all t, y=0,1" then 
            return 2*p;
        elif c = "all t,y" then 
            return p^2;
        elif c = "t=0,1,w, y=0,1" then 
            return 6;
        elif c = "t=1,w, all y" then
            return 2*p;

        # para = [z,t]
        elif c = "z=0,1,w, t ne 0,1,-1, t~t^-1" then 
            return 3*(p-3)/2;
        elif c = "t=1,-1, z=0,1,w" then 
            return 6;
        elif c = "all z,t" then 
            return p^2;
        fi;

    elif Length(para) = 3 then 
    
        # sort para
        if x in para and y in para and z in para then 
            para := [x,y,z];
        elif x in para and y in para and t in para then 
            para := [x,y,t];
        elif x in para and z in para and t in para then 
            para := [x,z,t];
        elif y in para and z in para and t in para then 
            para := [y,z,t];
        fi;

        # para = [x,y,z]
        if Length(c) = 0 or c = "all x,y,z" then 
            return p^3;
        elif c="all x,y, [x,y]~[y,x], z FIXED so that z^2-4 not a square mod p"
        then 
            return p*(p+1)/2;
        elif c = "all x,y, z=1,w,w^2,w^3,w^4, y~y' if y^5=y'^5 mod p" then 
            return 5*p*(1+(p-1)/g5);
        elif c = "x ne 0,all y,z, x~-x" then 
            return p^2*(p-1)/2;
        elif c = "all x,y,z, z~-z" then 
            return p^2*(p+1)/2;
        elif c = "all x,z, y ne 0, y~-y" then 
            return p^2*(p-1)/2;
        elif c = "x ne -2, all y,z, z~-z" then 
            return p*(p^2-1)/2;
        elif c = "x ne -2w, all y,z, z~-z" then 
            return p*(p^2-1)/2;
        elif c = "x ne 0, all y,z" then 
            return p^2*(p-1);
        elif c = "x ne 0, all y,z, x~-x" then 
            return p^2*(p-1)/2;
        elif c = "x,z ne 0, all y, y~-y, z~-z" then 
            return (p-1)*(p^2-1)/4;
        elif c = "x=0,1, y=0,1, z ne 0,-1" then 
            return 4*(p-2);
        elif c = "all x,y,z, [x,y,z]~[-x,-y,-z]" then 
            return (p^3+1)/2;
        elif c = "x ne 0, all y,z, x~x' if x^4=x'^4 mod p" then 
            #return p^2*(1+(p-1)/g4);
            return p^2*(p-1)/g4;
        elif c = "x ne 0, x~x' if x^3=x'^3 mod p, all y,z" then 
            return p^2*(p-1)/g3;
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p, all y,z" then 
            return p^2*(p-1)/g5;
        elif c = "x=w^i, 0 le i lt 2gcd[p-1,3], y ne 0, y~y' if y^3=y'^3 mod p, all z, z~-z" then 
            return (p^2-1);
        elif c = "y=w^i, i=2,3,4,5, x ne 0, x~x' if x^6=x'^6 mod p, all z" 
        then 
            return 2*p*(p-1)/g3;
        elif c = "z=w^i, i=1,2,3,4, x ne 0, x~x' if x^5=x'^5 mod p, all y" 
        then
            return 4*p*(p-1)/g5;
        elif c = "See Notes6.150, Case 3" then 
            return (p+1+(p+3)*g3+g4)/2;
        elif c = "See Notes6.150, Case 4" then 
            return 3+g3*(p+3+g4)/4;
        elif c = "See Notes6.150, Case 5" then 
            return 2+g3*(p+7-g4)/4;

        # para = [x,y,t]
        elif c = "all t, x=0,1, y=0,1" then 
            return 4*p;
        elif c = "t=0,1,w, x=0,1, y=0,1" then
            return 12;

        # para = [x,z,t]
        elif c = "t ne 0, all x,z, [x,z]~[tz,x/t]" then 
            return (p-1)*(p^2+p)/2;
        elif c="t ne 0,1, all x,z such that [x+t][1+z]=1 mod p, [x,z]~[tz,x/t]" 
        then 
            return (p^2-2*p-1)/2;

        # para = [y,z,t]
        elif c = "all y,z,t" then 
            return p^3;
        elif c = "all y,z,t, t~-t" then 
            return p^2*(p+1)/2;
        elif c = "y=0,1, all z,t" then 
            return 2*p^2;
        elif c = "y=0,1, z=0,1, t=0,1" then 
            return 8;
        elif c = "y=1,w, all z, t ne 0,1,-1, t~t^-1" then
            return p*(p-3);
        elif c = "all y,z,t, (z,t)~(t,z)" then 
            return p^2*(p+1)/2;
        fi;

    elif Length(para) = 4 then 
 
        # sort para
        para := [x,y,z,t];

        if Length(c) = 0 or c = "all x,y,z,t" then 
            return p^4;
        elif c = "all x,y,z,t, y~-y, if y=0 then t~-t" then 
            return (p^4+p^2)/2;
        elif c = "all y,z,t, x ne 0, x~-x" or
             c = "x ne 0, x~-x, all y,z,t" then 
            return p^3*(p-1)/2;
        elif c = "x ne 0, x~x' if x^5=x'^5 mod p, all y,z,t" then
            return p^3*(p-1)/g5;
        elif c = "all x,y,z,t, [x,y,z,t]~[t+1,z+1,y-1,x-1]" then 
            return (p^4+p^2)/2;
        elif d = 6 and c = "See note6.178" then 
            return p^2+(p+1-g4)/2;
        elif d = 7 and c = "See note5.38" and not w in S.param then 
            return (g3*(p^2+3*p+11)+1)/2;
        elif d = 7 and c = "See note5.38" and w in S.param then 
            return (g3*(p^2+p+1)+5)/2;
        elif d = 7 and c = "See Notes5.3, Case 4" then 
            return p^2+(p+1-g4)/2;
        elif c = "See Notes5.12" and w in S.param then 
            return p^2+(3*p+5)/2;
        elif c = "See Notes5.12" then 
            return p^2+(7*p+15)/2;
        elif c = "See Notes6.150, Case 8" then 
            return (5*p-7+(p^2-5)*g3-g4)/2;
        elif c = "See note2dec5.1" then 
            #return (p^2+2-g3)/2;
            return (p^2-2+h3)/2;
        elif c = "See Notes4.1 Case 4" then 
            return (p+1)^2/2;
        fi;

    elif Length(para) = 5 then 

        # sort para
        if x in para then
            para := [x,y,z,t,u];
        else
            para := [y,z,t,u,v];
        fi;

        # para = [x,y,z,t,u];
        if Length(c) = 0 then 
            return p^5;
        elif c = "See Notes6.163a" then  
            return 2*p^2-(5*p-1)/2;
        elif c = "See Notes6.163b" then  
            return (p^3-5*p+p*g4)/2;

        # para = [y,z,t,u,v];
        elif c = "See Notes6.173" then 
            return 3*p+3+(p^2+2*p+3)*g3/2;
        elif d = 7 and c = "See Notes6.178a" then 
            return (p-3)*(p+1)/2+a*(p-1)/2; 
        elif d = 7 and c = "See Notes6.178b" then 
            return p*(p-1)/2+(p-1)*(p+1-a)/2; 
        fi;

    elif Length(para) = 6 then 

        # sort para
        para := [x,y,z,t,u,v];

        # para = [x,y,z,t,u,v];
        if Length(c) = 0 then 
            return p^6;
        elif c = "See Notes5.14, Case 1" then 
            return 3*p+22;
        elif c = "See Notes5.14, Case 2" then 
            return 5*p+13+g3+g4;
        elif c = "See Notes5.14, Case 3" then 
            return 2*p+8+g4;
        elif c = "See Notes5.14, Case 4" then 
            return 6*p+8+2*g3+g4+g5;
        elif c = "See Notes5.14, Case 5" then 
            return 5*p+13+2*g3+g4;
        elif c = "See Notes5.14, Case 6" then 
            return p^2+3*p-3+(p+2)*g3+(p+1)*g4+(p+1)*g5;
        elif c = "See Notes5.14, Case 7" then 
            return 2*p^2+11*p+43+g4;
        elif c = "See Notes5.14, Case 8" then 
            return p^3+4*p^2+6*p+(p+5)*g3+3*g4+g5;
        elif c = "See Notes5.14, Case 9" then
            return p^3+5*p^2/2+7*p+19/2+(p+4)*g4/2;
        elif c = "See Notes5.14, Case 10" then 
            return p^3+5*p^2/2+7*p+19/2+(p+4)*g4/2;
        elif c = "See Notes5.14, Case 11" then
            return (p^4+p^3+4*p^2+p-1+(p^2+2*p+3)*g3+(p+2)*g4)/2;
        elif c = "See Notes5.14, Case 12" then 
            return (p^4+p^3+4*p^2+p-1+(p^2+2*p+3)*g3+(p+2)*g4)/2;
        elif c = "See Notes5.14, Case 13" then 
            return 2*p^2+11*p+27+g4;
        elif c = "See Notes5.14, Case 14" then 
            return p^3+2*p^2+6*p+10+(p+4)*g3;
        elif c = "See Notes5.14, Case 15" then 
            return p^3+(7*p^2+17*p+59+5*g3+(p+1)*g4)/2;
        elif c = "See Notes5.14, Case 16" then 
            return 2*p^4+4*p^3+8*p^2+14*p+11+4*g3+3*g4;
        elif c = "See Notes5.14, Case 17" then  
            return (p^4+2*p^3+3*p^2+4*p+2)*(p-1)/g3+3*p+4+(p^2+p+1)*g4/2;
        elif c = "See Notes5.14, Case 18" then 
            return (2*p^5+2*p^4+2*p^3+2*p^2+14*p+17)/3;
        elif c = "See Notes5.14, Case 19" then 
            return 2*p^2+11*p+27+g4;
        elif c = "See Notes5.14, Case 20" then 
            return 2*p^4+4*p^3+6*p^2+11*p+11+2*g3+(p+1)*g4;
        elif c = "See Notes5.14, Case 21" then 
            return 2*p^3+6*p^2+7*p+7+(p+1)*g4;
        elif c = "See Notes5.14, Case 22" then 
            return (2*p^3+3*p^2+3*p+13-g3+(p+1)*g4)/2;
        elif c = "See Notes5.14, Case 23" then 
            return (p^5+p^4+p^3+p^2)*g3/3+p+2+(p^2+p+1)*g4/2;
        fi;

    elif Length(para) = 7 then 

        # sort para
        para := [x,y,z,t,u,v,s];

        # para = [x,y,z,t,u,v,s];
        if c = "See Notes5.14, Case 24" then 
            return 2*(p^5+p^4+p^3+p^2)/3+2*p+3;
        fi;

    elif Length(para) = 8 then 

        # sort para
        para := [x,y,z,t,r,s,u,v];

        # para = [x,y,z,t,r,s,u,v];
        if Length(c) = 0 then 
            return p^8;
        elif c = "See Notes5.3, Case 6" then 
            return p^3/2+(2+g4)*p^2/2+(7*g4-9)*p/2+(-7-5*g3+7*g4+2*g3*g4)/2;
        elif c = "See Notes5.3, Case 7" then
            return p^3/2+(8-g4)*p^2/2+(33-7*g4)*p/2+(27-5*g4+7*g3-2*g3*g4)/2;
        fi;

    elif Length(para) = 12 then 

        # sort para
        para := [x,y,z,t,j,k,m,n,r,s,u,v];

        # para = [x,y,z,t,j,k,m,n,r,s,u,v];
        if Length(c) = 0 then 
            return p^12;
        elif c = "See Notes4.1, Case 5" then 
            return p^5+p^4+4*p^3+6*p^2+15*p+16+(p+1)*g3+(g3+h3-4)/2;
        elif c = "See Notes4.1, Case 6" then
            return (p^4+p^3+6*p^2+9*p+13)/2;
        fi;
    fi;

    return fail;
end;

NumberOfLiePRingsInFamily := function( L )
    local f, c, g, g3, g4, g5, g7, g8, g9, s16, h3;

    # set up
    f := NumberOfLiePRingsInFamilyNC(L);
    c := LibraryConditions(L)[2];
    g := 1;

    # get gcd's
    g3 := IndeterminateByName("(p-1,3)");
    h3 := IndeterminateByName("(p+1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");
    s16 := IndeterminateByName("(p^2-1,16)");

    # multiply
    if c = "p=1 mod 3" then
        g := (g3-1)/2;
    elif c = "p=2 mod 3" then
        g := -(g3-3)/2;
    elif c = "p=1 mod 4" then
        g := (g4-2)/2;
    elif c = "p=3 mod 4" then
        g := -(g4-4)/2;
    elif c = "p=1 mod 5" then
        g := (g5-1)/4;
    elif c = "p ne 1 mod 5" then
        g := -(g5-5)/4;
    elif c = "p=1 mod 7" then
        g := (g7-1)/6;
    elif c = "p=1 mod 8" then
        g := (g8-2)*(g8-4)/24;
    elif c = "p=1 mod 9" then
        g := (g9-1)*(g9-3)/48;
    elif c = "p=1 mod 4, p=1 mod 3" then
        g := (g3-1)*(g4-2)/4;
    elif c = "p=3 mod 4, p=1 mod 3" then
        g := -(g3-1)*(g4-4)/4;
    fi;

    f := EliminateDenominator( f*g );
    return SimplifyPorcPoly( f );

end;

