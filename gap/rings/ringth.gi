##
## l = LibraryConditions
## d = Dimension
## pp = Parameters
##

BindGlobal( "RingInvariantsByData", function( l, d, pp )
    local p, w, x, y, z, t, u, v, s, g3, g4, g8, units, zeros;

    p := IndeterminateByName("p");
    w := IndeterminateByName("w");
    x := IndeterminateByName("x");
    y := IndeterminateByName("y");
    z := IndeterminateByName("z");
    t := IndeterminateByName("t");
    u := IndeterminateByName("u");
    v := IndeterminateByName("v");
    s := IndeterminateByName("s");
    g3 := IndeterminateByName("(p-1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g8 := IndeterminateByName("(p-1,8)");

    units := [];
    zeros := [];

    if Length(pp) = 0 then 
        return rec( units := units, zeros := CallGroebner(zeros));
    elif Length(l)=0 then 
        return rec( units := units, zeros := CallGroebner(zeros));
    fi;

    #  special cases 
    if d=7 and l = "See Notes4.1 Case 4" then 
        Append(units, [x*t-y*z] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes4.1, Case 5" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes4.1, Case 6" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.12" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 1" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 10" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 11" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 12" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 13" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 14" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 15" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 16" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 17" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 18" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 19" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 2" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 20" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 21" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 22" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 23" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 24" then 
        Append(units, [s] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 3" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 4" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 5" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 6" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 7" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 8" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.14, Case 9" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.3, Case 4" then 
        Append(units, [x*t-y*z] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.3, Case 6" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes5.3, Case 7" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.150, Case 3" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.150, Case 4" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.150, Case 5" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.150, Case 8" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.163a" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.163b" then 
        Append(units, [x, u] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.173" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.178" then 
        Append(units, [y] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.178a" then 
        Append(units, [y, z] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See Notes6.178b" then 
        Append(units, [z] ); 
        Append(zeros, [t] ); 
    elif d=7 and l = "See note1dec5.1" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See note2dec5.1" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    elif d=7 and l = "See note5.38" then 
        Append(units, [] ); 
        Append(zeros, [] ); 
    fi; 

    # 1 parameters 
    if Length(pp) = 1 then 
        if l = "1+4x not a square" then 
            Append(units, [x, x-2] ); 
            Append(zeros, [] ); 
        elif l = "unique x so that 1-wx^2 is not a square" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "unique x so that x^2-w is not a square" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne -1" then 
            Append(units, [x+1] ); 
            Append(zeros, [] ); 
        elif l = "x ne -1,-1/2" then 
            Append(units, [x+1, x+1/2] ); 
            Append(zeros, [] ); 
        elif l = "x ne -2" then 
            Append(units, [x+2] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, equivalence classes {x,-x,1/x,-1/x}" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~-x" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~ax if a^3=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~ax if a^4=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~ax if a^5=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~ax if a^6=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~ax if a^7=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, x~x^-1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-1" then 
            Append(units, [x, x+1] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-1, x~-1-x" then 
            Append(units, [x, x+1] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-1,2,1/2" then 
            Append(units, [x, x+1, x-2, x-1/2] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-1/4" then 
            Append(units, [x, x+1/4] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-2, x~-x-2" then 
            Append(units, [x, x+2] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-w, x~-w-x" then 
            Append(units, [x, x+w] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,-w,2w,w/2" then 
            Append(units, [x, x+w, x-2*w, x-w/2] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0,1" then 
            Append(units, [x, x-1] ); 
            Append(zeros, [] ); 
        elif l = "x ne 1" then 
            Append(units, [x-1] ); 
            Append(zeros, [] ); 
        elif l = "x=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)] ); 
        elif l = "x=0,1,w" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)*(x-w)] ); 
        elif l = "x=0,1,w,w^2,w^3" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)*(x-w)*(x-w^2)*(x-w^3)] ); 
        elif l = "x=1,w,...,w^(2gcd(p-1,3)-1)" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x=1,w,...,w^(gcd(p-1,8)-1)" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x=w,w^2" then 
            Append(units, [x, x-1, x+1] ); 
            Append(zeros, [(x-w)*(x-w^2)] ); 
        elif l = "x=w,w^2,...,w^((p-3)/2)" then 
            Append(units, [x, x-1, x+1] ); 
            Append(zeros, [] ); 
        elif l = "x=w,w^2,...,w^6" then 
            Append(units, [x] ); 
            Append(zeros, [Product(List([1..6], a -> (x-w^a)))] ); 
        elif l = "x=w,w^2,w^3,w^4" then 
            Append(units, [x] ); 
            Append(zeros, [Product(List([1..4], a -> (x-w^a)))] ); 
        elif l = "x=w^2,w^3,w^4,w^5" then 
            Append(units, [x] ); 
            Append(zeros, [Product(List([2..5], a -> (x-w^a)))] ); 
        elif l = "x=w^3,w^4,...,w^8" then 
            Append(units, [x] ); 
            Append(zeros, [Product(List([3..8], a -> (x-w^a)))] ); 
        elif l = "x=w^4,w^5,w^6,w^7" then 
            Append(units, [x] ); 
            Append(zeros, [Product(List([4..7], a -> (x-w^a)))] ); 
        elif l = "x~-1-x" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x~-x" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x~-x-2" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x~1-x" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x~ax if a^3=1" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x~ax if a^4=1" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x~w-x" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        fi;
    fi; 

    # 2 parameters 
    if Length(pp) = 2 then 
        if l = "[x,y]~[-x,y]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y]~[ax,y] if a^3=1" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y]~[x',y'] if y^2-wx^2=y'^2-wx'^2" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y]~[x,-y]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y]~[x,ay] if a^4=1" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y]~[y,x]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y]~[±x,ay] if a^3=1" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x ne -1, (1+x)y=1" then 
            Append(units, [y, x+1] ); 
            Append(zeros, [(1+x)*y - 1] ); 
        elif l = "x ne -1,3, See Notes6.114" then 
            Append(units, [x+1, x-3] ); 
            Append(zeros, [] ); 
        elif l = "x ne -2" then 
            Append(units, [x+2] ); 
            Append(zeros, [] ); 
        elif l = "x ne -2w" then 
            Append(units, [x+2*w] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, 4x+y^2 not a square" then 
            Append(units, [x, (x-1)+y-(x-1)*y, (x-2)+(y-1)-(x-2)*(y-1)] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[-x,-y-2]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[-x,-y]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[-x,y]~[x,-y]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[a^4x,ay] if a^5=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[ax,a^2y] if a^3=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[ax,a^2y] if a^6=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[ax,a^3y] if a^5=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[ax,a^3y] if a^6=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[ax,ay] if a^3=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[x,-y]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y]~[x,-y]~[-x,iy] if i^2=-1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, unique y so that 1-wy^2 is not a square, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, unique y so that wy^2=2, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [w*y^2-2] ); 
        elif l = "x ne 0, unique y so that y^2=2, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [y^2-2] ); 
        elif l = "x ne 0, unique y so that y^2=2w^2, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [y^2-2*w^2] ); 
        elif l = "x ne 0, unique y so that y^2=2w^3, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [y^2-2*w^3] ); 
        elif l = "x ne 0, y=1,-1" then 
            Append(units, [x, y] ); 
            Append(zeros, [(y-1)*(y+1)] ); 
        elif l = "x ne 0, y=1,-1, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [(y-1)*(y+1)] ); 
        elif l = "x ne 0, y=1,w, [x,y]~[-x,y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [(y-1)*(y-w)] ); 
        elif l = "x ne 0, y=w,w^2,...,w^((p-3)/2)" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, y=w,w^2,...,w^6, [x,y]~[ax,y] if a^7=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [Product(List([1..6], a -> (y-w^a)))] ); 
        elif l = "x ne 0, y=w,w^2,w^3,w^4, [x,y]~[ax,y] if a^5=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [Product(List([1..4], a -> (y-w^a)))] ); 
        elif l = "x ne 0, y=w^2,w^3,w^4,w^5, [x,y]~[ax,y] if a^3=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [Product(List([2..5], a -> (y-w^a)))] ); 
        elif l = "x ne 1" then 
            Append(units, [x-1] ); 
            Append(zeros, [] ); 
        elif l = "x ne 1-w, [x,y]~[x,-y]~[-x+2(1-w),iy] if i^2=-1" then 
            Append(units, [x-(1-w)] ); 
            Append(zeros, [] ); 
        elif l = "x ne 1-w^2, [x,y]~[x,-y]~[-x+2(1-w^2),iy] if i^2=-1" then 
            Append(units, [x-(1-w^2)] ); 
            Append(zeros, [] ); 
        elif l = "x ne 1-w^3, [x,y]~[x,-y]~[-x+2(1-w^3),iy] if i^2=-1" then 
            Append(units, [x-(1-w^3)] ); 
            Append(zeros, [] ); 
        elif l = "x,y ne 0, [x,y]~[-x,-y]" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x,y ne 0, [x,y]~[ax,a^3y] if a^4=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x,y ne 0, [x,y]~[xy^-2,y^-1]" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x,y ne 0, [x,y]~[y,x]" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x,y ne 0, x ne y, [x,y]~[y,x]" then 
            Append(units, [x, y, x-y] ); 
            Append(zeros, [] ); 
        elif l = "x=0,-1, y=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x+1), y*(y-1)] ); 
        elif l = "x=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)] ); 
        elif l = "x=0,1, y=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1), y*(y-1)] ); 
        elif l = "x=0,1,w" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)*(x-w)] ); 
        elif l = "x=0,1,w, y=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)*(x-w), y*(y-1)] ); 
        elif l = "x=0,1,w, y=1,-1" then 
            Append(units, [y] ); 
            Append(zeros, [x*(x-1)*(x-w), (y+1)*(y-1)] ); 
        elif l = "x=0,1,w, y=w,w^2,...,w^((p-3)/2)" then 
            Append(units, [y] ); 
            Append(zeros, [x*(x-1)*(x-w)] ); 
        elif l = "x=1,w" then 
            Append(units, [x] ); 
            Append(zeros, [(x-1)*(x-w)] ); 
        elif l = "x=1,w,...,w^(2gcd(p-1,3)-1), y ne 0, [x,y]~[x,ay] if a^6=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x=1,w,...,w^(gcd(p-1,8)-1), y ne 0, [x,y]~[x,ay] if a^8=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "xy ne 1, [x,y]~[y,x]" then 
            Append(units, [x*y-1] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[-x,y]" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[a^2x,ay] if a^4=1" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[x',y'] if b(wyy'+xx'-w)=a(x-x') and b(xy'+yx')=a(y-y') with a^2-wb^2 ne 0" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[x,-y]" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[x,ay] if a^3=1" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[x,ay] if a^4=1" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y]~[x,ay] if a^5=1" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 1/2, [x,y]~[-x,1-y]" then 
            Append(units, [y-1/2] ); 
            Append(zeros, [] ); 
        elif l = "y=1,w, [x,y]~[-x,y]" then 
            Append(units, [y] ); 
            Append(zeros, [(y-1)*(y-w)] ); 
        elif l = "y=1,w, [x,y]~[ax,y] if a^8=1" then 
            Append(units, [y] ); 
            Append(zeros, [(y-1)*(y-w)] ); 
        elif l = "y=1,w,...,w^5, [x,y]~[ax,y] if a^6=1" then 
            Append(units, [y] ); 
            Append(zeros, [Product(List([0..5], a -> (y-w^a)))] ); 
        elif l = "y=1,w,w^2, [x,y]~[ax,y] if a^3=1" then 
            Append(units, [y] ); 
            Append(zeros, [(y-1)*(y-w)*(y-w^2)] ); 
        elif l = "y=w^2,w^3, [x,y]~[ax,y] if a^8=1" then 
            Append(units, [y] ); 
            Append(zeros, [(y-w^2)*(y-w^3)] ); 
        elif l = "y=w^2,w^3,w^4,w^5" then 
            Append(units, [y] ); 
            Append(zeros, [Product(List([2..5], a -> (y-w^a)))] ); 
        elif l = "y=w^4,w^5,w^6,w^7, [x,y]~[ax,y] if a^8=1" then 
            Append(units, [y] ); 
            Append(zeros, [Product(List([4..7], a -> (y-w^a)))] ); 
        fi;
    fi; 

    # 3 parameters 
    if Length(pp) = 3 then 
        if l = "[x,y,z]~[-x,-y,-z]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y,z]~[-x,y,z]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y,z]~[x,y,-z]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y,z]~[z,y,x]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "unique z so that z^2-4 is not a square, [x,y,z]~[y,x,z]" then 
            Append(units, [z] ); 
            Append(zeros, [] ); 
        elif l = "x ne -2, [x,y,z]~[x,y,-z]" then 
            Append(units, [x+2] ); 
            Append(zeros, [] ); 
        elif l = "x ne -2w, [x,y,z]~[x,y,-z]" then 
            Append(units, [x+2*w] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z]~[-x,-y,z]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z]~[a^3x,a^4y,az] if a^5=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z]~[ax,a^2y,az] if a^4=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z]~[ax,y,a^2z] if a^3=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, y=w^2,w^3,w^4,w^5, [x,y,z]~[ax,y,a^2z] if a^6=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [Product(List([2..5], a -> (y-w^a)))] ); 
        elif l = "x ne 0, z=w,w^2,w^3,w^4, [x,y,z]~[ax,a^3y,z] if a^5=1" then 
            Append(units, [x, z] ); 
            Append(zeros, [Product(List([1..4], a -> (z-w^a)))] ); 
        elif l = "x,z ne 0, [x,y,z]~[-x,y,-z]~[x,-y,z]" then 
            Append(units, [x,z] ); 
            Append(zeros, [] ); 
        elif l = "x=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1)] ); 
        elif l = "x=0,1, y=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1), y*(y-1)] ); 
        elif l = "x=0,1, y=0,1, z ne 0,-1" then 
            Append(units, [z, z+1] ); 
            Append(zeros, [x*(x-1), y*(y-1)] ); 
        elif l = "x=0,1, y=0,1, z=0,1" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1), y*(y-1), z*(z-1)] ); 
        elif l = "x=0,1, y=0,1, z=0,1,w" then 
            Append(units, [] ); 
            Append(zeros, [x*(x-1), y*(y-1), z*(z-1)*(z-w)] ); 
        elif l = "x=1,w,...,w^(2gcd(p-1,3)-1), y ne 0, [x,y,z]~[x,ay,±a^2z] if a^3=1" then 
            Append(units, [x, y] ); 
            Append(zeros, [] ); 
        elif l = "x=w,w^2,...,w^((p-3)/2), y=1,w" then 
            Append(units, [x, y] ); 
            Append(zeros, [(y-1)*(y-w)] ); 
        elif l = "y ne 0, [x,y,z]~[x,-y,-z]" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0, [x,y,z]~[zy,y,x/y]" then 
            Append(units, [y] ); 
            Append(zeros, [] ); 
        elif l = "y ne 0,1, (x+y)(1+z)=1, [x,y,z]~[zy,y,x/y]" then 
            Append(units, [y, y-1, x+y, z+1] ); 
            Append(zeros, [(x+y)*(1+z)-1] ); 
        elif l = "z=1,w,w^2,w^3,w^4, [x,y,z]~[x,ay,z] if a^5=1" then 
            Append(units, [z] ); 
            Append(zeros, [Product(List([0..4], a -> (z-w^a)))] ); 
        fi;
    fi; 

    # 4 parameters 
    if Length(pp) = 4 then 
        if l = "A=[[t,x],[y,z]] non-singular, A~A' if det(P)*A*P=P*A' with P=[[a,b],[ewb,ea]] non-singular, e=±1" then 
            Append(units, [t*z-x*y] ); 
            Append(zeros, [] ); 
        elif l = "[x,y,z,t]~[t+1,z+1,y-1,x-1]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "[x,y,z,t]~[x,-y,z,-t]" then 
            Append(units, [] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z,t]~[-x,y,t,z]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z,t]~[-x,y,z,-t]" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z,t]~[ax,a^3y,a^4z,at] if a^5=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        elif l = "x ne 0, [x,y,z,t]~[ax,a^3y,a^4z,t] if a^5=1" then 
            Append(units, [x] ); 
            Append(zeros, [] ); 
        fi;
    fi; 

    return rec( units := units, zeros := CallGroebner(zeros));
end );



