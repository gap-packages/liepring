
LIE_TABLE := [];

CallFuncList( function()
    local p, g3, g4, g5, g6, g7, g8, g9;
    p := IndeterminateByName("p");
    g3 := IndeterminateByName("(p-1,3)");
    g4 := IndeterminateByName("(p-1,4)");
    g5 := IndeterminateByName("(p-1,5)");
    g7 := IndeterminateByName("(p-1,7)");
    g8 := IndeterminateByName("(p-1,8)");
    g9 := IndeterminateByName("(p-1,9)");

LIE_TABLE := [

["gapdec6.0", 1],
["2gen/gapdec4.7", 4 ],
["2gen/gapdec5.37", p^2+8*p+25 ],
["2gen/gapdec5.38", p+6+(p^2+3*p+10)*g3 ],
["2gen/gapdec5.60", 2*p^2+p+3+2*(p+1)*g3 + (2*p+4)*g4 + g8 ],
["2gen/gapdec5.65", p^3+p^2+p-2+2*g3+g4+(p+1)*g5 ],
["2gen/gapdec6.366", 2 ],
["2gen/gapdec6.367", 2 ],
["2gen/gapdec6.368", 9+6*g3 ],
["2gen/gapdec6.369", 5+g3 ],
["2gen/gapdec6.370", 2*p+2+(p+1)/2+g3+g4/2 ],
["2gen/gapdec6.371", 2*p+2+(p+1)/2+g3+g4/2 ],
["2gen/gapdec6.372", (p+1)/2+g3+g4/2 ],
["2gen/gapdec6.373", p ],
["2gen/gapdec6.374", p ],
["2gen/gapdec6.375", (p+1)/2+g3+g4/2 ],
["2gen/gapdec6.376", p ],
["2gen/gapdec6.377", p ],
["2gen/gapdec6.378", p+1+g3 ],
["2gen/gapdec6.379", 5+2*g3+g4 ],
["2gen/gapdec6.380", p ],
["2gen/gapdec6.381", 6+2*g3+g4 ],
["2gen/gapdec6.382", p+5+(p-1)/2 ],
["2gen/gapdec6.383", p+5+(p-1)/2 ],
["2gen/gapdec6.384", 2 ],
["2gen/gapdec6.386", p+5 ],
["2gen/gapdec6.388", 4*p+5+2*g3+g4 ],
["2gen/gapdec6.394", p^2+3*p+10+2*g3+3*g4+2*g5+g8+g9 ],
["2gen/gapdec6.399", (1/2)*p^3+(3/2)*p^2+4*p+7+3*g3+(p+3)*g4/2+g5+g8/2 ],
["2gen/gapdec6.404", (p^3+p^2+2*p+2+(p+1)*g4+g8)/2 ],
["2gen/gapdec6.408", (1/2)*(p+1)*(p-1+g3) ],
["2gen/gapdec6.411", 2+2*g3+g4 ],
["2gen/gapdec6.412", p-1+(1/2)*(p+1)*g4 ],
["2gen/gapdec6.414", (1/2)*(p+1)*(p-1+g3) ],
["2gen/gapdec6.417", 2+2*g3+g4 ],
["2gen/gapdec6.418", p-1+(1/2)*(p+1)*g4 ],
["2gen/gapdec6.420", 5*p+1+g4 ],
["2gen/gapdec6.421", 2*p^2-p-1+(p+1)*g3 ],
["2gen/gapdec6.424", p^2 ],
["2gen/gapdec6.425", 3 ],
["2gen/gapdec6.426", p ],
["2gen/gapdec6.427", p+1 ],
["2gen/gapdec6.428", 2 ],
["2gen/gapdec6.429", 1 ],
["2gen/gapdec6.431", 4 ],
["2gen/gapdec6.435", 7+2*g3+2*g4 ],
["2gen/gapdec6.436", 2+2*g3+g4+2*g5 ],
["2gen/gapdec6.442", 2*p-2+(p+2)*g3 ],
["2gen/gapdec6.445", 2*p-2+g3 ],
["2gen/gapdec6.448", (p+1)/2 ],
["2gen/gapdec6.451", (p+1)/2 ],
["2gen/gapdec6.454", 8+2*g3+4*g4 ],
["2gen/gapdec6.455", 4 ],
["2gen/gapdec6.459", 4*p+2*g3+g4+2*g5 ],
["2gen/gapdec6.460", 3*p-3+g4 ],
["2gen/gapdec6.467", 2 ],
["2gen/gapdec6.469", 4+2*g3+3*g5+g8 ],
["2gen/gapdec6.475", 4*p-1+g5+g7 ],
["2gen/gapdec6.507", 2*p^2+p+2*p*g3+p*g5 ],
["2gen/gapdec6.518", 2 ],
["3gen/gapdec3.1", p+14 ],
["3gen/gapdec5.10", 2*p+7 ],
["3gen/gapdec5.12", 3*p^2+17*p+53+g3+g4 ],
["3gen/gapdec5.14", 2*p^5+7*p^4+19*p^3+49*p^2+128*p+256+(p^2+7*p+29)*g3+(p^2+7*p+24)*g4+(p+3)*g5 ],
["3gen/gapdec5.15", 3*p^2 + 12*p + 14 + (p+2)*g4 ],
["3gen/gapdec5.16", p^4 + 2*p^3 + 5*p^2 + 14*p ],
["3gen/gapdec5.18", 3*p^3 + 6*p^2 + 6*p + 11 + (p+7)*g3 + (p+1)*g4 + g5 ],
["3gen/gapdec5.8", p+8 ],
["3gen/gapdec5.9", 4*p^2+26*p+107+5*g3+(p+4)*g4 ],
["3gen/gapdec6.100", 3 ],
["3gen/gapdec6.101", 3*(p-1)/2 ],
["3gen/gapdec6.102", p+3 ],
["3gen/gapdec6.103", p+3 ],
["3gen/gapdec6.104", 5*p + 24 + 3*g3 ],
["3gen/gapdec6.105", 11 + 5*g3 + g5 ],
["3gen/gapdec6.106", p^2 + 10*p + 32 + g3 + g4 ],
["3gen/gapdec6.108", (p-1)*(12 + g3)/2 + g4/2 ],
["3gen/gapdec6.109", 7+2*g3 ],
["3gen/gapdec6.110", 2+g3+g4/2 ],
["3gen/gapdec6.111", (p-1)*(4 + g3)/2 ],
["3gen/gapdec6.112", p^2 + 4*p + 3 + 2*g3 ],
["3gen/gapdec6.113", 5*p+4 ],
["3gen/gapdec6.114", 4*p-4 ],
["3gen/gapdec6.115", 3 ],
["3gen/gapdec6.116", 3 ],
["3gen/gapdec6.117", p+1 ],
["3gen/gapdec6.118", 11+4*g3 ],
["3gen/gapdec6.119", (p-1)/2 + 3 + 2*g3 + g4 ],
["3gen/gapdec6.120", (p-1)/2 + 3 + 2*g3 + g4 ],
["3gen/gapdec6.121", 2*p+4+2*g3 ],
["3gen/gapdec6.122", p+2 ],
["3gen/gapdec6.125", p+1 ],
["3gen/gapdec6.127", 3*p + 4 + 6*g3 + g4 ],
["3gen/gapdec6.131", 15 + (p+10)*g3 + g4 + g7],
["3gen/gapdec6.132", (p + 1 + 3*(p+1)*g3 + g4)/2 ],
["3gen/gapdec6.133", (p + 1 + 3*(p+1)*g3 + g4)/2 ],
["3gen/gapdec6.134", 3*p-1+g3 ],
["3gen/gapdec6.135", p^2 + 2*p + 3 + g3 + g4 + g5 ],
["3gen/gapdec6.138", p*(p+1)/2 ],
["3gen/gapdec6.139", (p+3)/2 + g3 ],
["3gen/gapdec6.140", (3*p+1)/2 ],
["3gen/gapdec6.142", p*(p+1)/2 ],
["3gen/gapdec6.143", (p+3)/2+g3 ],
["3gen/gapdec6.144", (3*p+1)/2 ],
["3gen/gapdec6.146", 2*p^2 + 3*p ],
["3gen/gapdec6.148", p^3+p^2+p+(p+2)*g3+g5 ],
["3gen/gapdec6.150", 4*p+14+(p^2/2+2*p+13/2)*g3+g4 ],
["3gen/gapdec6.151", p^2+p+2+(p+1)*g3 ],
["3gen/gapdec6.152", (p^2 + 2*p + 3 + 2*g3 + (p+1)*g4)/2 ],
["3gen/gapdec6.153", p*(p+1)/2 ],
["3gen/gapdec6.154", (p^2 + 2*p + 3 + 2*g3 + (p+1)*g4)/2 ],
["3gen/gapdec6.155", p*(p+1)/2 ],
["3gen/gapdec6.156", p ],
["3gen/gapdec6.157", 2*p - 1 ],
["3gen/gapdec6.158", p ],
["3gen/gapdec6.159", (p^3 + p^2)/2 ],
["3gen/gapdec6.160", (p^3 + p^2)/2 ],
["3gen/gapdec6.160a", (p+3)/2 ],
["3gen/gapdec6.161", 2*p - 1 ],
["3gen/gapdec6.162", 2*p-1 ],
["3gen/gapdec6.163", (p^3 + 4*p^2 - 10*p + 1 + p*g4)/2 ],
["3gen/gapdec6.168", p^3 + 2*p^2 + 2*p + 2 + g3 ],
["3gen/gapdec6.172", (p^4 + p^2)/2 ],
["3gen/gapdec6.173", 3*p + 3 + (p^2+2*p+3)*g3/2 ],
["3gen/gapdec6.174", (p^3 + p^2 + p + 1)/4 ],
["3gen/gapdec6.175", (p^3 + p^2 + p + 1)/4 ],
["3gen/gapdec6.176", p*(p-1) + p*g4/2 ],
["3gen/gapdec6.178", (3*p^2-1)/2 ],
["3gen/gapdec6.179", (p^4+p^2)/2 ],
["3gen/gapdec6.182", 4 ],
["3gen/gapdec6.183", 2 ],
["3gen/gapdec6.184", 9 ],
["3gen/gapdec6.187", 11+4*g3+2*g4 ],
["3gen/gapdec6.188", 3 + 2*g3 ],
["3gen/gapdec6.189", (5*p + 7)/2 ],
["3gen/gapdec6.190", 2 ],
["3gen/gapdec6.191", (5*p + 7)/2 ],
["3gen/gapdec6.192", 2 ],
["3gen/gapdec6.197", 11 + 4*g3 + 2*g4 ],
["3gen/gapdec6.198", 4 + g3 ],
["3gen/gapdec6.207", 5 ],
["3gen/gapdec6.212", 4 ],
["3gen/gapdec6.215", 3 ],
["3gen/gapdec6.216", 4 ],
["3gen/gapdec6.218", 11 + 4*g3 + 2*g4 ],
["3gen/gapdec6.222", 2*p + 5 + 3*g3 + g4 ],
["3gen/gapdec6.228", p+1 ],
["3gen/gapdec6.231", p^2+5*p+14+(p+17)*g3+2*g4+g7 + g8 ],
["3gen/gapdec6.256", 20 + (p+11)*g3 + 4*g4 ],
["3gen/gapdec6.261", 4*p + 2 + (p^2+3*p+1)*g3 + (p+1)*g4 ],
["3gen/gapdec6.267", p + 3 + 3*g3 + g4 + g5 ],
["3gen/gapdec6.269", (p + 5)/2 ],
["3gen/gapdec6.271", (p+5)/2 ],
["3gen/gapdec6.273", 1 + g3 + g4 ],
["3gen/gapdec6.274", p + 2 + (2*p+3)*g3 + g5 ],
["3gen/gapdec6.275", 1 ],
["3gen/gapdec6.276", 1 ],
["3gen/gapdec6.277", p ],
["3gen/gapdec6.278", (5*p + 3)/2 ],
["3gen/gapdec6.279", p ],
["3gen/gapdec6.280", (5*p + 3)/2 ],
["3gen/gapdec6.281", 2*p^2 + 4*p + 4 + 2*g3 ],
["3gen/gapdec6.282", p^3 + p^2 + 2*p + 2 + g3 ],
["3gen/gapdec6.289", 6*p ],
["3gen/gapdec6.290", 2*p^2+p ],
["3gen/gapdec6.294", p + 4 + 5*g3 + g4 ],
["3gen/gapdec6.295", (p + 5)/2 ],
["3gen/gapdec6.296", (p+5)/2 ],
["3gen/gapdec6.297", 2*p + p*g3 ],
["3gen/gapdec6.298", (3*p - 1)/2 ],
["3gen/gapdec6.299", (3*p-1)/2 ],
["3gen/gapdec6.303", p^2 + p + (p+1)*g3 + g4 ],
["3gen/gapdec6.304", p*(p+1)/2 ],
["3gen/gapdec6.305", p*(p+1)/2 ],
["3gen/gapdec6.312", 4*p + 2 + 2*g3 + 4*g4 ],
["3gen/gapdec6.313", p + g4 + g5 ],
["3gen/gapdec6.322", 3 ],
["3gen/gapdec6.325", 19 + 5*g3 + 6*g4 ],
["3gen/gapdec6.326", 4*p + 5 + 4*g3 + 2*g4 + 4*g5 ],
["3gen/gapdec6.327", 3*p + 12 + 2*g3 + 3*g4 + g7 ],
["3gen/gapdec6.328", p^2 + 3*p - 2 + (p+3)*g3 + 2*g4 + 2*g5 ],
["3gen/gapdec6.362", p^2 + 7*p + 3 + 2*g3 + 3*g4 + g5 ],
["3gen/gapdec6.85", 3 ],
["3gen/gapdec6.86", p+17 ],
["3gen/gapdec6.87", 5 ],
["3gen/gapdec6.88", 26 ],
["3gen/gapdec6.89", p+6 ],
["3gen/gapdec6.90", 30 ],
["3gen/gapdec6.91", 2*p + 88 + g4 ],
["3gen/gapdec6.92", 2*p + 13 ],
["3gen/gapdec6.93", p + 15 + 2*g3 + g4 ],
["3gen/gapdec6.94", 2*p + 15 + 3*g3 + g4 ],
["3gen/gapdec6.95", 5*p + 10 ],
["3gen/gapdec6.96", 2*p + 26 ],
["3gen/gapdec6.97", 3*p + 18 + g3 ],
["3gen/gapdec6.98", (5*p + 1)/2 ],
["3gen/gapdec6.99", p+4 ],
["4gen/gapdec4.1", p^5+2*p^4+7*p^3+25*p^2+88*p+270+(p+4)*g3+g4],
["4gen/gapdec5.3", p^4+5*p^3+19*p^2+64*p+140+(p+6)*g3+(p+7)*g4+g5 ],
["4gen/gapdec6.10", 8 ],
["4gen/gapdec6.11", 39 ],
["4gen/gapdec6.12", 8 ],
["4gen/gapdec6.13", 11 ],
["4gen/gapdec6.14", 39 ],
["4gen/gapdec6.15", 2*p+6 ],
["4gen/gapdec6.16", 6 ],
["4gen/gapdec6.17", p+7 ],
["4gen/gapdec6.18", 3*p+5 ],
["4gen/gapdec6.19", 97 ],
["4gen/gapdec6.20", 6*p+35 ],
["4gen/gapdec6.21", 13*p+27 ],
["4gen/gapdec6.23", 3*p+41+8*g3+2*g4 ],
["4gen/gapdec6.24", 5 ],
["4gen/gapdec6.29", 10*p+69+g3+g4 ],
["4gen/gapdec6.33", 39+g3 ],
["4gen/gapdec6.34", 5*p+38+2*g4 ],
["4gen/gapdec6.35", 10*p+26 ],
["4gen/gapdec6.36", p^3+9*p^2+20*p+18+g3+g4 ],
["4gen/gapdec6.48", 27+3*g3+2*g4 ],
["4gen/gapdec6.51", 8*p+15+(2*p+5)*g3+(p+2)*g4+g5 ],
["4gen/gapdec6.52", 4*p^2+15*p+15+(p+1)*g3+g4 ],
["4gen/gapdec6.60", 7+g3 ],
["4gen/gapdec6.63", 4 ],
["4gen/gapdec6.67", 18+8*g3+3*g4 ],
["4gen/gapdec6.72", 17+(p+7)*g3+g4 ],
["4gen/gapdec6.9", 7 ],
["5gen/gapdec5.1", p^2+15*p+125 ],
["5gen/gapdec6.2", 5],
["5gen/gapdec6.3", 25],
["6gen/gapdec6.1", 9],
["gap7.1", 1] ];
end, []);

InstallGlobalFunction( LiePRingsDim7ByFile, function( arg )
    local rem, L, spe, v;

    # back up
    rem := StructuralCopy(LIE_DATA[7]);

    # read and construct the desired lie p-rings
    LIE_DATA[7] := [];
    LiePRing_ReadPackage(Concatenation("lib/dim7/",LIE_TABLE[arg[1]][1]));
    L := List(LIE_DATA[7], x -> LiePRingByData(7, x));

    # restore
    LIE_DATA[7] := rem;

    # this is it if there is no prime given
    if Length(arg) = 1 then return L; fi;

    # add primes
    spe := List( L, X -> LiePRingsInFamily( X, arg[2] ) );
    spe := Flat(Filtered(spe, l -> l <> fail ));

    # check
    v := EvaluatePorcPoly( LIE_TABLE[arg[1]][2], arg[2] );
    if Length(spe) <> v then
        Print("WARNING: ",Length(spe)-v," too many \n");
        return fail;
    fi;
    return spe;
end );

BindGlobal( "GetDim7FileByNumber", function( P, nr )
    local i, j, v;
    i := 0;
    j := 0;
    while j < nr do
        i := i+1;
        v := EvaluatePorcPoly( LIE_TABLE[i][2], P );
        j := j + v;
    od;
    return [i, nr-j+v];
end );

