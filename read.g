#############################################################################
##
##
##

OFFSET_VARS := 1000;
BindGlobal( "IndeterminateByName", function( name )
    local l, i, v;
    l := ["p", "w",
          "x", "y", "z", "t", "j", "k", "m", "n", "r", "s", "u", "v",
          "(p-1,3)", "(p-1,4)", "(p-1,5)", "(p-1,7)", "(p-1,8)", "(p-1,9)",
          "(p+1,3)", "(p^2-1,16)", "a"];
    i := Position( l, name );
    if i = fail then return fail; fi;
    v := Indeterminate(Integers, OFFSET_VARS+i);
    if not HasName(v) then SetName(v,name); fi;
    return v;
end );


ReadPackage( "liepring", "gap/singular.gi");


BindGlobal( "LiePRing_ReadPackage", function(relpath)
    local preamble, code, c, filename, func;

    preamble := "local p, w, x, y, z, t, j, k, m, n, r, s, u, v;\n";
    for c in "pwxyztjkmnrsuv" do
        # c := IndeterminateByName("c");
        Add(preamble, c);
        Append(preamble, " := IndeterminateByName(\"");
        Add(preamble, c);
        Append(preamble, "\");\n");
    od;

    filename:= Filename( DirectoriesPackageLibrary( "liepring", "" ), relpath );
    code := StringFile( filename );
    code := Concatenation(preamble, code);

    # Parse this into a function
    func := ReadAsFunction( InputTextString( code ) );

    # Execute the function to populate LIE_DATA
    func();
end );

#############################################################################
##
## read data
##

LIE_DATA := [];
LiePRing_ReadPackage( "lib/data1.gi"); ## dimension 1 - 4
LiePRing_ReadPackage( "lib/data2.gi"); ## dimension 5

LIE_DATA[6] := [];
LiePRing_ReadPackage( "lib/dim6/gapdec5.0"); # 1
LiePRing_ReadPackage( "lib/dim6/gapdec4.6"); # 2
LiePRing_ReadPackage( "lib/dim6/gapdec4.7"); # p+15
LiePRing_ReadPackage( "lib/dim6/gapdec4.8"); # 1
LiePRing_ReadPackage( "lib/dim6/gapdec5.37"); # p+8
LiePRing_ReadPackage( "lib/dim6/gapdec5.38"); # 5+3(p-1,3)
LiePRing_ReadPackage( "lib/dim6/gapdec5.39"); # 1+(p-1,3)+(1/2)(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.40"); # 1+(p-1,3)+(1/2)(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.41"); # p+1+(p-1,3)
LiePRing_ReadPackage( "lib/dim6/gapdec5.42"); # p+1
LiePRing_ReadPackage( "lib/dim6/gapdec5.45"); # p
LiePRing_ReadPackage( "lib/dim6/gapdec5.47"); # 2
LiePRing_ReadPackage( "lib/dim6/gapdec5.48"); # 1
LiePRing_ReadPackage( "lib/dim6/gapdec5.49"); # 4
LiePRing_ReadPackage( "lib/dim6/gapdec5.50"); # 3+2(p-1,3)+(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.51"); # 3(p+1)/2
LiePRing_ReadPackage( "lib/dim6/gapdec5.52"); # 3(p+1)/2
LiePRing_ReadPackage( "lib/dim6/gapdec5.54"); # 3+2(p-1,3)+(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.58"); # 2
LiePRing_ReadPackage( "lib/dim6/gapdec5.60"); # 7+2(p-1,3)+3(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.65"); # 2p+2(p-1,3)+(p-1,4)+2(p-1,5)
LiePRing_ReadPackage( "lib/dim6/gapdec5.73"); # 2
LiePRing_ReadPackage( "lib/dim6/gapdec3.1"); # 3p+27 
LiePRing_ReadPackage( "lib/dim6/gapdec4.3"); # 3p^2+13p+37+(p-1,3)+(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.8"); # 4
LiePRing_ReadPackage( "lib/dim6/gapdec5.9"); # 23
LiePRing_ReadPackage( "lib/dim6/gapdec5.10"); # 5
LiePRing_ReadPackage( "lib/dim6/gapdec5.11"); # 4
LiePRing_ReadPackage( "lib/dim6/gapdec5.12"); # 12
LiePRing_ReadPackage( "lib/dim6/gapdec5.13"); # p+1
LiePRing_ReadPackage( "lib/dim6/gapdec5.14"); # 35
LiePRing_ReadPackage( "lib/dim6/gapdec5.15"); # 2p+13
LiePRing_ReadPackage( "lib/dim6/gapdec5.16"); # 4p+8
LiePRing_ReadPackage( "lib/dim6/gapdec5.18"); # 2p+13+3(p-1,3)+(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.19"); # 3
LiePRing_ReadPackage( "lib/dim6/gapdec5.24"); # 3
LiePRing_ReadPackage( "lib/dim6/gapdec5.27"); # 11+4(p-1,3)+2(p-1,4)
LiePRing_ReadPackage( "lib/dim6/gapdec5.32"); # 4+2(p-1,3)
LiePRing_ReadPackage( "lib/dim6/gapdec4.1"); # 4p+48
LiePRing_ReadPackage( "lib/dim6/gapdec5.2"); # 4
LiePRing_ReadPackage( "lib/dim6/gapdec5.3"); # 18
LiePRing_ReadPackage( "lib/dim6/gapdec5.1"); # 7
LiePRing_ReadPackage( "lib/dim6/gap6.1"); # 1 

LIE_DATA[7] := [];
LiePRing_ReadPackage( "lib/dim7/gapdec6.0");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec4.7");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec5.37");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec5.38");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec5.60");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec5.65");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.366");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.367");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.368");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.369");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.370");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.371");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.372");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.373");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.374");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.375");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.376");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.377");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.378");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.379");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.380");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.381");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.382");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.383");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.384");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.386");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.388");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.394");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.399");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.404");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.408");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.411");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.412");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.414");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.417");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.418");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.420");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.421");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.424");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.425");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.426");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.427");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.428");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.429");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.431");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.435");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.436");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.442");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.445");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.448");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.451");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.454");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.455");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.459");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.460");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.467");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.469");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.475");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.507");
LiePRing_ReadPackage( "lib/dim7/2gen/gapdec6.518");

LiePRing_ReadPackage( "lib/dim7/3gen/gapdec3.1");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.10");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.12");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.14");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.15");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.16");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.18");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.8");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec5.9");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.100");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.101");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.102");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.103");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.104");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.105");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.106");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.108");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.109");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.110");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.111");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.112");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.113");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.114");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.115");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.116");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.117");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.118");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.119");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.120");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.121");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.122");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.125");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.127");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.131");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.132");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.133");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.134");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.135");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.138");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.139");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.140");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.142");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.143");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.144");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.146");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.148");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.150");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.151");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.152");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.153");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.154");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.155");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.156");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.157");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.158");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.159");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.160");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.160a");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.161");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.162");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.163");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.168");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.172");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.173");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.174");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.175");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.176");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.178");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.179");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.182");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.183");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.184");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.187");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.188");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.189");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.190");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.191");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.192");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.197");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.198");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.207");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.212");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.215");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.216");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.218");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.222");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.228");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.231");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.256");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.261");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.267");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.269");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.271");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.273");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.274");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.275");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.276");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.277");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.278");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.279");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.280");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.281");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.282");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.289");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.290");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.294");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.295");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.296");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.297");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.298");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.299");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.303");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.304");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.305");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.312");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.313");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.322");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.325");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.326");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.327");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.328");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.362");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.85");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.86");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.87");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.88");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.89");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.90");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.91");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.92");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.93");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.94");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.95");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.96");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.97");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.98");
LiePRing_ReadPackage( "lib/dim7/3gen/gapdec6.99");

LiePRing_ReadPackage( "lib/dim7/4gen/gapdec4.1");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec5.3");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.10");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.11");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.12");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.13");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.14");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.15");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.16");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.17");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.18");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.19");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.20");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.21");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.23");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.24");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.29");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.33");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.34");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.35");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.36");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.48");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.51");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.52");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.60");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.63");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.67");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.72");
LiePRing_ReadPackage( "lib/dim7/4gen/gapdec6.9");

LiePRing_ReadPackage( "lib/dim7/5gen/gapdec5.1");
LiePRing_ReadPackage( "lib/dim7/5gen/gapdec6.2");
LiePRing_ReadPackage( "lib/dim7/5gen/gapdec6.3");

LiePRing_ReadPackage( "lib/dim7/6gen/gapdec6.1");
LiePRing_ReadPackage( "lib/dim7/gap7.1");

#############################################################################
##
## read functions
##
ReadPackage( "liepring", "gap/rings/polys.gi");
ReadPackage( "liepring", "gap/rings/zeruni.gi");
ReadPackage( "liepring", "gap/rings/ringth.gi");
ReadPackage( "liepring", "gap/rings/number.gi");

ReadPackage( "liepring", "gap/basic/general.gi");
ReadPackage( "liepring", "gap/basic/echelon.gi"); # used in evals
ReadPackage( "liepring", "gap/basic/linalg.gi");  # used in cover
ReadPackage( "liepring", "gap/basic/collect.gi");
ReadPackage( "liepring", "gap/basic/lieelms.gi");
ReadPackage( "liepring", "gap/basic/liering.gi");
ReadPackage( "liepring", "gap/basic/subring.gi");
ReadPackage( "liepring", "gap/basic/series.gi");


ReadPackage( "liepring", "gap/evals/valsfun.gi");
ReadPackage( "liepring", "gap/evals/indx27a.gi");
ReadPackage( "liepring", "gap/evals/reps27a.gi");
ReadPackage( "liepring", "gap/evals/vals27a.gi");
ReadPackage( "liepring", "gap/evals/vals27b.gi");
ReadPackage( "liepring", "gap/evals/vals28.gi");
ReadPackage( "liepring", "gap/evals/special.gi");

ReadPackage( "liepring", "gap/class/pgroup.gi");
ReadPackage( "liepring", "gap/class/porcpoly.gi");
ReadPackage( "liepring", "lib/data.gi");
ReadPackage( "liepring", "lib/table.gi");

ReadPackage( "liepring", "gap/advan/elements.gi");
ReadPackage( "liepring", "gap/advan/autos.gi");
ReadPackage( "liepring", "gap/advan/schur.gi");
ReadPackage( "liepring", "gap/advan/cover.gi");


#############################################################################
##
## add-on: the lie-p-rings of dimension 8 and maximal class
##
LIE_DATA_MC8 := [];
LiePRing_ReadPackage( "lib/dim8/maxclass");
