
BindGlobal( "EvaluatePorcPoly", function( f, P )
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
end );

BindGlobal( "DegreePorcPoly", function( f )
    local p, a, b;
    if IsInt(f) then return 0; fi;
    p := IndeterminateByName("p");
    a := DegreeIndeterminate(NumeratorOfRationalFunction(f),p);
    b := DegreeIndeterminate(DenominatorOfRationalFunction(f),p);
    if b <> 0 then Error("something wrong with degrees "); fi;
    return a;
end );

BindGlobal( "EliminateDenominator", function( f )
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
end );

BindGlobal( "SimplifyPorcPoly", function( f )
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
end );
    
BindGlobal( "CharFunctionByIsomRestrictions", function( L )
    local d, c, p, l, para, w, x, y, z, t, j, k, m, n, r, s, u, v, 
          g3, h3, g4, g5, g7, g8, g9, s16, a;

    # get some data
    d := DimensionOfLiePRing(L);
    c := LibraryConditions(L);
    p := PrimeOfLiePRing(L);
    l := LibraryName(L);
    para := ParametersOfLiePRing(L);

    # the trivial cases
    if Length(para) = 0 then return 1; fi;
    if Length(c[1]) = 0 then return p^Length(para); fi;

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

    ## a = number of non-zero elms in GF(p) with 2*w^2*a^3+5*w*a^2+2*a square
    a := IndeterminateByName("a");

    if l = "6.178" then # See note6.178 
        return p^2+(p+1-g4)/2;
    elif l = "6.62" then # See note6.62 
        return p;
    elif l = "7.62" then # See note5.38 
        return (g3*(p^2+3*p+11)+1)/2;
    elif l = "7.63" then # See note5.38 
        return (g3*(p^2+p+1)+5)/2;
    elif l = "7.729" then # See Notes5.12 
        return p^2+(7*p+15)/2;
    elif l = "7.730" then # See Notes5.12 
        return p^2+(3*p+5)/2;
    elif l = "7.757" then # See Notes5.14, Case 1 
        return 3*p+22;
    elif l = "7.758" then # See Notes5.14, Case 2 
        return 5*p+13+g3+g4;
    elif l = "7.759" then # See Notes5.14, Case 3 
        return 2*p+8+g4;
    elif l = "7.760" then # See Notes5.14, Case 4 
        return 6*p+8+2*g3+g4+g5;
    elif l = "7.761" then # See Notes5.14, Case 5 
        return 5*p+13+2*g3+g4;
    elif l = "7.762" then # See Notes5.14, Case 6 
        return p^2+3*p-3+(p+2)*g3+(p+1)*g4+(p+1)*g5;
    elif l = "7.763" then # See Notes5.14, Case 7 
        return 2*p^2+11*p+43+g4;
    elif l = "7.764" then # See Notes5.14, Case 8 
        return p^3+4*p^2+6*p+(p+5)*g3+3*g4+g5;
    elif l = "7.765" then # See Notes5.14, Case 9 
        return p^3+5*p^2/2+7*p+19/2+(p+4)*g4/2;
    elif l = "7.766" then # See Notes5.14, Case 10 
        return p^3+5*p^2/2+7*p+19/2+(p+4)*g4/2;
    elif l = "7.767" then # See Notes5.14, Case 11 
        return (p^4+p^3+4*p^2+p-1+(p^2+2*p+3)*g3+(p+2)*g4)/2;
    elif l = "7.768" then # See Notes5.14, Case 12 
        return (p^4+p^3+4*p^2+p-1+(p^2+2*p+3)*g3+(p+2)*g4)/2;
    elif l = "7.769" then # See Notes5.14, Case 13 
        return 2*p^2+11*p+27+g4;
    elif l = "7.770" then # See Notes5.14, Case 14 
        return p^3+2*p^2+6*p+10+(p+4)*g3;
    elif l = "7.771" then # See Notes5.14, Case 15 
        return p^3+(7*p^2+17*p+59+5*g3+(p+1)*g4)/2;
    elif l = "7.772" then # See Notes5.14, Case 16 
        return 2*p^4+4*p^3+8*p^2+14*p+11+4*g3+3*g4;
    elif l = "7.773" then # See Notes5.14, Case 17 
        return (p^4+2*p^3+3*p^2+4*p+2)*(p-1)/g3+3*p+4+(p^2+p+1)*g4/2;
    elif l = "7.774" then # See Notes5.14, Case 18 
        return (2*p^5+2*p^4+2*p^3+2*p^2+14*p+17)/3;
    elif l = "7.775" then # See Notes5.14, Case 19 
        return 2*p^2+11*p+27+g4;
    elif l = "7.776" then # See Notes5.14, Case 20 
        return 2*p^4+4*p^3+6*p^2+11*p+11+2*g3+(p+1)*g4;
    elif l = "7.777" then # See Notes5.14, Case 21 
        return 2*p^3+6*p^2+7*p+7+(p+1)*g4;
    elif l = "7.778" then # See Notes5.14, Case 22 
        return (2*p^3+3*p^2+3*p+13-g3+(p+1)*g4)/2;
    elif l = "7.779" then # See Notes5.14, Case 23 
        return (p^5+p^4+p^3+p^2)*g3/3+p+2+(p^2+p+1)*g4/2;
    elif l = "7.780" then # See Notes5.14, Case 24 
        return 2*(p^5+p^4+p^3+p^2)/3+2*p+3;
    elif l = "7.1691" then # See Notes6.150, Case 3 
       return (p+1+(p+3)*g3+g4)/2;
    elif l = "7.1692" then # See Notes6.150, Case 4 
        return 3+g3*(p+3+g4)/4;
    elif l = "7.1693" then # See Notes6.150, Case 5 
        return 2+g3*(p+7-g4)/4;
    elif l = "7.1709" then # See Notes6.150, Case 8 
        return (5*p-7+(p^2-5)*g3-g4)/2;
    elif l = "7.1763" then # See Notes6.163a 
        return 2*p^2-(5*p-1)/2;
    elif l = "7.1764" then # See Notes6.163b 
        return (p^3-5*p+p*g4)/2;
    elif l = "7.1777" then # See Notes6.173 
        return 3*p+3+(p^2+2*p+3)*g3/2;
    elif l = "7.1797" then # See Notes6.178 
        return (p+1)/2;
    elif l = "7.1798" then # See Notes6.178a 
        return (p-3)*(p+1)/2+a*(p-1)/2;
    elif l = "7.1799" then # See Notes6.178b 
        return p*(p-1)/2+(p-1)*(p+1-a)/2;
    elif l = "7.3068" then # See Notes4.1 Case 4 
        return (p+1)^2/2;
    elif l = "7.3285" then # See Notes4.1, Case 5 
       return p^5+p^4+4*p^3+6*p^2+15*p+16+(p+1)*g3+(g3+h3-4)/2;
    elif l = "7.3286" then # See Notes4.1, Case 6 
        return (p^4+p^3+6*p^2+9*p+13)/2;
    elif l = "7.3387" then # See Notes5.3, Case 4 
        return p^2+(p+1-g4)/2;
    elif l = "7.3437" then # See Notes5.3, Case 6 
        return p^3/2+(2+g4)*p^2/2+(7*g4-9)*p/2+(-7-5*g3+7*g4+2*g3*g4)/2;
    elif l = "7.3438" then # See Notes5.3, Case 7 
        return p^3/2+(8-g4)*p^2/2+(33-7*g4)*p/2+(27-5*g4+7*g3-2*g3*g4)/2;
    elif l = "7.4670" then # See note1dec5.1 
        return p;
    elif l = "7.4723" then # See note2dec5.1 
        return (p^2-2+h3)/2;
    elif l = "7.1389" then # See Notes6.114
        return 3*p-7;
    fi;

    c := c[1];

    # 1 parameter
    if Length(para) = 1 then 

        if c = "1+4x not a square" then
            return (p-1)/2;
        elif c = "unique x so that 1-wx^2 is not a square" then
            return 1;
        elif c = "unique x so that x^2-w is not a square" then
            return 1;
        elif c = "x ne -1" then
            return p-1;
        elif c = "x ne -1,-1/2" then
            return p-2;
        elif c = "x ne -2" then
            return p-1;
        elif c = "x ne 0" then
            return p-1;
        elif c = "x ne 0, equivalence classes {x,-x,1/x,-1/x}" then
            return (p-1+g4)/4;
        elif c = "x ne 0, x~-x" then
            return (p-1)/2;
        elif c = "x ne 0, x~x^-1" then
            return (p+1)/2;
        elif c = "x ne 0,-2, x~-x-2" then
            return (p-1)/2;
        elif c = "x ne 0, x~ax if a^3=1" then
            return (p-1)/g3;
        elif c = "x ne 0, x~ax if a^4=1" then
            return (p-1)/g4;
        elif c = "x ne 0, x~ax if a^5=1" then
            return (p-1)/g5;
        elif c = "x ne 0, x~ax if a^6=1" then
            return (p-1)/(2*g3);
        elif c = "x ne 0, x~ax if a^7=1" then
            return (p-1)/g7;
        elif c = "x ne 0,-1" then
            return p-2;
        elif c = "x ne 0,-1, x~-1-x" then
            return (p-1)/2;
        elif c = "x ne 0,-1,2,1/2" then
            return p-4;
        elif c = "x ne 0,-1/4" then
            return p-2;
        elif c = "x ne 0,-w, x~-w-x" then
            return (p-1)/2;
        elif c = "x ne 0,-w,2w,w/2" then
            return p-4;
        elif c = "x ne 0,1" then
            return p-2;
        elif c = "x ne 1" then
            return p-1;
        elif c = "x~-1-x" then
            return (p+1)/2;
        elif c = "x~-x" then
            return (p+1)/2;
        elif c = "x~-x-2" then
            return (p+1)/2;
        elif c = "x~1-x" then
            return (p+1)/2;
        elif c = "x~w-x" then
            return (p+1)/2;
        elif c = "x~ax if a^3=1" then
            return 1+(p-1)/g3;
        elif c = "x~ax if a^4=1" then
            return 1+(p-1)/g4;
        elif c = "x=0,1" then
            return 2;
        elif c = "x=0,1,w" then
            return 3;
        elif c = "x=0,1,w,w^2,w^3" then
            return 5;
        elif c = "x=1,w,...,w^(2gcd(p-1,3)-1)" then
            return 2*g3;
        elif c = "x=1,w,...,w^(gcd(p-1,8)-1)" then
            return g8;
        elif c = "x=w,w^2" then
            return 2;
        elif c = "x=w,w^2,...,w^((p-3)/2)" then
            return (p-3)/2;
        elif c = "x=w,w^2,...,w^6" then
            return 6;
        elif c = "x=w,w^2,w^3,w^4" then
            return 4;
        elif c = "x=w^2,w^3,w^4,w^5" then
            return 4;
        elif c = "x=w^3,w^4,...,w^8" then
            return 6;
        elif c = "x=w^4,w^5,w^6,w^7" then
            return 4;
        fi;
    fi;
  
    # 2 parameter - todo
    if Length(para) = 2 then 

       if c = "[x,y]~[-x,y]" then
            return p*(p+1)/2; 
        elif c = "[x,y]~[x,-y]" then
            return p*(p+1)/2; 
        elif c = "[x,y]~[y,x]" then
            return p*(p+1)/2;
        elif c = "[x,y]~[x',y'] if y^2-wx^2=y'^2-wx'^2" then
            return p;
        elif c = "[x,y]~[ax,y] if a^3=1" then
            return p*(1+(p-1)/g3);
        elif c = "[x,y]~[x,ay] if a^4=1" then
            return p*(1+(p-1)/g4);
        elif c = "[x,y]~[±x,ay] if a^3=1" then
            return (p+1)/2*(1+(p-1)/g3);
        elif c = "x ne -1, (1+x)y=1" then
            return p-1;
        elif c = "x ne -2" then
            return p*(p-1);
        elif c = "x ne -2w" then
            return p*(p-1);
        elif c = "x ne 0" then
            return p*(p-1);
        elif c = "x ne 0, 4x+y^2 not a square" then
            return p*(p-1)/2;
        elif c = "x ne 0, [x,y]~[-x,-y-2]" then
            return p*(p-1)/2;
        elif c = "x ne 0, [x,y]~[-x,-y]" then
            return p*(p-1)/2;
        elif c = "x ne 0, [x,y]~[-x,y]~[x,-y]" then
            return (p-1)*(p+1)/4;
        elif c = "x ne 0, [x,y]~[a^4x,ay] if a^5=1" then
            return p*(p-1)/g5;
        elif c = "x ne 0, [x,y]~[ax,a^2y] if a^3=1" then
            return p*(p-1)/g3;
        elif c = "x ne 0, [x,y]~[ax,a^2y] if a^6=1" then
            return p*(p-1)/(2*g3);
        elif c = "x ne 0, [x,y]~[ax,a^3y] if a^5=1" then
            return p*(p-1)/g5;
        elif c = "x ne 0, [x,y]~[ax,a^3y] if a^6=1" then
            return p*(p-1)/(2*g3);
        elif c = "x ne 0, [x,y]~[ax,ay] if a^3=1" then
            return p*(p-1)/g3;
        elif c = "x ne 0, [x,y]~[x,-y]" then
            return (p-1)*(p+1)/2;
        elif c = "x ne 0, [x,y]~[x,-y]~[-x,iy] if i^2=-1" then
            return (p-1)*(p+1)/4;
        elif c = "x ne 0, unique y so that 1-wy^2 is not a square, [x,y]~[-x,y]" then
            return (p-1)/2;
        elif c = "x ne 0, unique y so that wy^2=2, [x,y]~[-x,y]" then
            return (p-1)/2;
        elif c = "x ne 0, unique y so that y^2=2, [x,y]~[-x,y]" then
            return (p-1)/2;
        elif c = "x ne 0, unique y so that y^2=2w^2, [x,y]~[-x,y]" then
            return (p-1)/2;
        elif c = "x ne 0, unique y so that y^2=2w^3, [x,y]~[-x,y]" then
            return (p-1)/2;
        elif c = "x ne 0, y=1,-1" then
            return 2*(p-1);
        elif c = "x ne 0, y=1,-1, [x,y]~[-x,y]" then
            return (p-1);
        elif c = "x ne 0, y=1,w, [x,y]~[-x,y]" then
            return (p-1);
        elif c = "x ne 0, y=w,w^2,...,w^((p-3)/2)" then
            return (p-1)*(p-3)/2;
        elif c = "x ne 0, y=w,w^2,...,w^6, [x,y]~[ax,y] if a^7=1" then
            return 6*(p-1)/g7;
        elif c = "x ne 0, y=w,w^2,w^3,w^4, [x,y]~[ax,y] if a^5=1" then
            return 4*(p-1)/g5;
        elif c = "x ne 0, y=w^2,w^3,w^4,w^5, [x,y]~[ax,y] if a^3=1" then
            return 4*(p-1)/g3;
        elif c = "x ne 1" then
            return p*(p-1);
        elif c = "x ne 1-w, [x,y]~[x,-y]~[-x+2(1-w),iy] if i^2=-1" then
            return (p+1)*(p-1)/4;
        elif c = "x ne 1-w^2, [x,y]~[x,-y]~[-x+2(1-w^2),iy] if i^2=-1" then
            return (p+1)*(p-1)/4;
        elif c = "x ne 1-w^3, [x,y]~[x,-y]~[-x+2(1-w^3),iy] if i^2=-1" then
            return (p+1)*(p-1)/4;
        elif c = "x,y ne 0, [x,y]~[-x,-y]" then
            return (p-1)^2/2;
        elif c = "x,y ne 0, [x,y]~[ax,a^3y] if a^4=1" then
            return (p-1)^2/g4;
        elif c = "x,y ne 0, [x,y]~[xy^-2,y^-1]" then
            return (p-1)*(p+1)/2;
        elif c = "x,y ne 0, [x,y]~[y,x]" then
            return p*(p-1)/2;
        elif c = "x,y ne 0, x ne y, [x,y]~[y,x]" then
            return (p-1)*(p-2)/2;
        elif c = "x=0,-1, y=0,1" then
            return 4;
        elif c = "x=0,1" then
            return 2*p;
        elif c = "x=0,1, y=0,1" then
            return 4;
        elif c = "x=0,1,w" then
            return 3*p;
        elif c = "x=0,1,w, y=0,1" then
            return 6;
        elif c = "x=0,1,w, y=1,-1" then
            return 6;
        elif c = "x=0,1,w, y=w,w^2,...,w^((p-3)/2)" then
            return 3*(p-3)/2;
        elif c = "x=1,w" then
            return 2*p;
        elif c = "x=1,w,...,w^(2gcd(p-1,3)-1), y ne 0, [x,y]~[x,ay] if a^6=1" then
            return (p-1);
        elif c = "x=1,w,...,w^(gcd(p-1,8)-1), y ne 0, [x,y]~[x,ay] if a^8=1" then
            return (p-1);
        elif c = "xy ne 1, [x,y]~[y,x]" then
            return (p^2-1)/2;
        elif c = "y ne 0, [x,y]~[-x,y]" then
            return (p-1)*(p+1)/2;
        elif c = "y ne 0, [x,y]~[a^2x,ay] if a^4=1" then
            return p*((p-1)/g4);
        elif c = "y ne 0, [x,y]~[x,-y]" then
            return p*(p-1)/2;
        elif c = "y ne 0, [x,y]~[x,ay] if a^3=1" then
            return p*((p-1)/g3);
        elif c = "y ne 0, [x,y]~[x,ay] if a^4=1" then
            return p*((p-1)/g4);
        elif c = "y ne 0, [x,y]~[x,ay] if a^5=1" then
            return p*((p-1)/g5);
        elif c = "y ne 1/2, [x,y]~[-x,1-y]" then
            return p*(p-1)/2;
        elif c = "y=1,w, [x,y]~[-x,y]" then
            return (p+1);
        elif c = "y=1,w, [x,y]~[ax,y] if a^8=1" then
            return 2*(1+(p-1)/g8);
        elif c = "y=1,w,...,w^5, [x,y]~[ax,y] if a^6=1" then
            return 6*(1+(p-1)/(2*g3));
        elif c = "y=1,w,w^2, [x,y]~[ax,y] if a^3=1" then
            return 3*(1+(p-1)/g3);
        elif c = "y=w^2,w^3, [x,y]~[ax,y] if a^8=1" then
            return 2*(1+(p-1)/g8);
        elif c = "y=w^2,w^3,w^4,w^5" then
            return 4*p;
        elif c = "y=w^4,w^5,w^6,w^7, [x,y]~[ax,y] if a^8=1" then
            return 4*(1+(p-1)/g8);
        fi;
    fi;

    # 3 parameter
    if Length(para) = 3 then 
    
        if c = "[x,y,z]~[-x,-y,-z]" then
            return (p^3+1)/2;
        elif c = "[x,y,z]~[-x,y,z]" then
            return p^2*(p+1)/2;
        elif c = "[x,y,z]~[x,y,-z]" then
            return p^2*(p+1)/2;
        elif c = "[x,y,z]~[z,y,x]" then
            return p^2*(p+1)/2;
        elif c = "unique z so that z^2-4 is not a square, [x,y,z]~[y,x,z]" then
            return p*(p+1)/2;
        elif c = "x ne -2, [x,y,z]~[x,y,-z]" then
            return p*(p^2-1)/2;
        elif c = "x ne -2w, [x,y,z]~[x,y,-z]" then
            return p*(p^2-1)/2;
        elif c = "x ne 0" then
            return p^2*(p-1);
        elif c = "x ne 0, [x,y,z]~[-x,-y,z]" then
            return p^2*(p-1)/2;
        elif c = "x ne 0, [x,y,z]~[a^3x,a^4y,az] if a^5=1" then
            return p^2*(p-1)/g5;
        elif c = "x ne 0, [x,y,z]~[ax,a^2y,az] if a^4=1" then
            return p^2*(p-1)/g4;
        elif c = "x ne 0, [x,y,z]~[ax,y,a^2z] if a^3=1" then
            return p^2*(p-1)/g3;
        elif c = "x ne 0, y=w^2,w^3,w^4,w^5, [x,y,z]~[ax,y,a^2z] if a^6=1" then
            return 2*p*(p-1)/g3;
        elif c = "x ne 0, z=w,w^2,w^3,w^4, [x,y,z]~[ax,a^3y,z] if a^5=1" then
            return 4*p*(p-1)/g5;
        elif c = "x,z ne 0, [x,y,z]~[-x,y,-z]~[x,-y,z]" then
            return (p-1)*(p^2-1)/4;
        elif c = "x=0,1" then
            return 2*p^2;
        elif c = "x=0,1, y=0,1" then
            return 4*p;
        elif c = "x=0,1, y=0,1, z ne 0,-1" then
            return 4*(p-2);
        elif c = "x=0,1, y=0,1, z=0,1" then
            return 8;
        elif c = "x=0,1, y=0,1, z=0,1,w" then
            return 12;
        elif c = "x=1,w,...,w^(2gcd(p-1,3)-1), y ne 0, [x,y,z]~[x,ay,±a^2z] if a^3=1" then
            return (p-1)*(p+1);
        elif c = "x=w,w^2,...,w^((p-3)/2), y=1,w" then
            return p*(p-3);
        elif c = "y ne 0, [x,y,z]~[x,-y,-z]" then
            return p^2*(p-1)/2;
        elif c = "y ne 0, [x,y,z]~[zy,y,x/y]" then
            return p*(p-1)*(p+1)/2;
        elif c = "y ne 0,1, (x+y)(1+z)=1, [x,y,z]~[zy,y,x/y]" then
            return (p^2-2*p-1)/2;
        elif c = "z=1,w,w^2,w^3,w^4, [x,y,z]~[x,ay,z] if a^5=1" then
            return 5*p*(1+(p-1)/g5);
        fi;
    fi;

    # 4 parameters
    if Length(para) = 4 then 
        if c = "[x,y,z,t]~[t+1,z+1,y-1,x-1]" then
            return (p^4+p^2)/2;
        elif c = "[x,y,z,t]~[x,-y,z,-t]" then
            return (p^4+p^2)/2;
        elif c = "x ne 0, [x,y,z,t]~[-x,y,t,z]" then
            return p^3*(p-1)/2;
        elif c = "x ne 0, [x,y,z,t]~[-x,y,z,-t]" then
            return p^3*(p-1)/2;
        elif c = "x ne 0, [x,y,z,t]~[ax,a^3y,a^4z,at] if a^5=1" then
            return p^3*(p-1)/g5;
        elif c = "x ne 0, [x,y,z,t]~[ax,a^3y,a^4z,t] if a^5=1" then
            return p^3*(p-1)/g5;
        fi;
    fi;
end );

BindGlobal( "CharFunctionByPrimeRestrictions", function( L )

    local g3, g4, g5, g7, g8, g9, s16, h3, c;

    # get gcd's
    g3 := IndeterminateByName("(p-1,3)");
    h3 := IndeterminateByName("(p+1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");
    s16 := IndeterminateByName("(p^2-1,16)");

    c := LibraryConditions(L)[2];

    # check cases
    if c = "p=1 mod 3" then
        return (g3-1)/2;
    elif c = "p=2 mod 3" then
        return -(g3-3)/2;
    elif c = "p=1 mod 4" then
        return (g4-2)/2;
    elif c = "p=3 mod 4" then
        return -(g4-4)/2;
    elif c = "p=1 mod 5" then
        return (g5-1)/4;
    elif c = "p ne 1 mod 5" then
        return -(g5-5)/4;
    elif c = "p=1 mod 7" then
        return (g7-1)/6;
    elif c = "p=1 mod 8" then
        return (g8-2)*(g8-4)/24;
    elif c = "p=1 mod 9" then
        return (g9-1)*(g9-3)/48;
    elif c = "p=1 mod 4, p=1 mod 3" then
        return (g3-1)*(g4-2)/4;
    elif c = "p=3 mod 4, p=1 mod 3" then
        return -(g3-1)*(g4-4)/4;
    elif c = "p=±1 mod 8" then
        return (s16-8)/8;
    elif c = "p=±3 mod 8" then
        return (16-s16)/8;
    fi;
    return 1;
end );

BindGlobal( "NumberOfLiePRingsInFamily", function( L )
    local a, b, f;
    a := CharFunctionByIsomRestrictions(L);
    b := CharFunctionByPrimeRestrictions(L);
    f := EliminateDenominator(a*b);
    return SimplifyPorcPoly(f);
end );
