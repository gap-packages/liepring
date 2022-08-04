
##
## This file is mainly due to Willem de Graaf.
##
##
## A is the tree corr to BCH, we evaluate the formula in x, y elements 
## of the Lie ring L
##

BindGlobal( "NewEvalBCH", function( A, tors, x, y )
   local a, k, val, coef, e, T, T0, i, vall, valr;

   a := x + y;
   k:= Maximum( tors ); 
   val := x * y;
   e := A.label;
   if not k = infinity then e := e mod k; fi;
   a := a + e * val;
   T := [ [ A, val ] ];
   while Length( T ) > 0 do
      T0 := [ ];
      for i in [ 1 .. Length( T ) ] do
         if T[ i ][ 1 ].isleaf = false then
            val := T[ i ][ 2 ];

            if not IsZero( val ) then
               vall := x * val;
               if not IsZero(vall) then
                  e := T[ i ][ 1 ].left.label;
                  if not k = infinity then e := e mod k; fi;
                  a := a + e * vall;
               fi;
               Add( T0, [ T[ i ][ 1 ].left, vall ] );
               valr := y * val;
               if not IsZero( valr ) then
                  e := T[ i ][ 1 ].right.label;
                  if not k = infinity then e := e mod k; fi;
                  a := a + e * valr;
               fi;
               Add( T0, [ T[ i ][ 1 ].right, valr ] );
            fi;
         fi;
      od;
      T := T0;
   od;
   return a;
end );

##
## Determine the p-group corresponding to L
##
InstallGlobalFunction( "PGroupByLiePRing", function( L )
   local G_BCH, p, dim, F, rels, i, u, c, l, x0, y0, z0, w, j, G, g, GtoL, 
         LtoG, v1, v2, v3, gens, tors, k, s;
 
   if not IsParentLiePRing(L) then return fail; fi;
   if not IsBound( LRPrivateFunctions.LAZARDTrec.G_BCH ) then
      LRPrivateFunctions.InitialiseLAZARD();
   fi;
   G_BCH:= LRPrivateFunctions.LAZARDTrec.G_BCH;
 
   # set up
   p := PrimeOfLiePRing(L);
   gens := BasisOfLiePRing(L);
   dim := Length( gens );

   # check
   if PClassOfLiePRing(L) > p-1 then return fail; fi;

   # compute torsion
   tors:= [ ];
   for i in [1..Length(gens)] do
       k:= 1; while not IsZero(p^k*gens[i]) do k:= k+1; od;
       Add( tors, p^k );
   od;

   # construct group
   F := FreeGroup( dim );
   rels:= [ ];
   for i in [1..dim-1] do
      u := Exponents( p*gens[i] );
      c := [ ];
      for l in [ 1 .. dim ] do
         Add( c, u[ l ] );
         y0:= u*gens;
         x0 := - u[ l ] * gens[l];
         z0 := NewEvalBCH( G_BCH, tors, x0, y0 );
         u := Exponents( z0 );
      od;

      w := Product( [ 1 .. dim ], s -> F.( s )^( c[ s ] ) );
      Add( rels, F.( i )^p/w );
   od;
   Add( rels, F.( dim )^p );
   for i in [ 1 .. dim ] do
      for j in [ i + 1 .. dim ] do
         v1 := NewEvalBCH( G_BCH, tors, gens[ j ], gens[ i ] );
         v2 := NewEvalBCH( G_BCH, tors, -gens[ j ], -gens[ i ] );
         v3 := NewEvalBCH( G_BCH, tors, v2, v1 );
         u := Exponents(v3);
         c := [ ];

         for l in [ 1 .. dim ] do
            Add( c, u[ l ] );
            y0 := u*gens;
            x0 := - u[ l ] * gens[l];
            z0 := NewEvalBCH( G_BCH, tors, x0, y0 );
            u := Exponents( z0 );
         od;

         w := Product( [ 1 .. dim ], s -> F.( s )^( c[ s ] ) );
         w := F.( j ) * w;
         Add( rels, F.( j )^F.( i )/( w ) );
      od;
   od;

   return PcGroupFpGroupNC( F/rels );
end );

BindGlobal( "GroupToLiePRing", function( L, G, g0 )
    local b, cf, x0, i;
    b := BasisOfLiePRing(L);
    cf:= ExponentsOfPcElement( Pcgs(G), g0 );
    x0:= cf[1]*b[1];
    for i in [2..Length(b)] do
        x0:= NewEvalBCH(LRPrivateFunctions.LAZARDTrec.G_BCH,L,x0,cf[i]*b[i]);
    od;
    return x0;
end );

BindGlobal( "LiePRingToGroup", function( L, G, x0 )
    local b, g, cf, exps, i;
    b := BasisOfLiePRing(L);
    g := Pcgs(G);
    cf := ExponentsLPR(L, x0);
    exps:= [ ];
    for i in [1..Length(b)] do
        Add( exps, cf[i] mod RelativeOrders(g)[i] );
        x0:= NewEvalBCH(LRPrivateFunctions.LAZARDTrec.G_BCH,L,-cf[i]*b[i],x0);
        cf:= ExponentsLPR(L, x0);
    od;
    return PcElementByExponents(g, exps);
end );

BindGlobal( "PGroupByLiePRing_Old", function( L )
    local d, p, l, F, f, r, i, j, K;

    d := DimensionOfLiePRing(L);
    p := PrimeOfLiePRing(L);
    l := GeneratorsOfRing(L);
    if not IsInt(p) then return fail; fi;

    F := FreeLieRing( Integers, d );
    f := GeneratorsOfLeftOperatorRing( F );
    r := [];

    for i in [1..d] do
        Add(r, p*f[i] - LinearCombination(Exponents(p*l[i]), f));
    od;

    for i in [1..d] do
        for j in [i+1..d] do
            Add(r, f[i]*f[j] - LinearCombination(Exponents(l[i]*l[j]), f));
        od;
    od;

    K := FpLieRing( F, r );

    return LieRingToPGroup(K).pgroup;
end );

BindGlobal( "GroupsViaLiePRings", function( arg )
    local d, P, L, G, i;
    d := arg[1];
    P := arg[2];
    L := LiePRingsByLibrary(d);
    G := List(L, x -> true);
    for i in [1..Length(L)] do
        G[i] := LiePRingsInFamily(L[i], P);
        if not IsBool(G[i]) then 
            G[i] := List( G[i], x -> PGroupByLiePRing(x) );
            G[i] := Filtered( G[i], x -> not IsBool(x) );
            if Length(arg) = 3 and arg[3] = true then 
                G[i] := List( G[i], x -> CodePcGroup(x) );
            fi;
        fi;
    od;
    G := Flat(G);
    G := Filtered(G, x -> not IsBool(x));
    return G;
end );

BindGlobal( "PrintGroupsViaLiePRings", function( arg )
    local d, P, L, G, i, j;
    d := arg[1];
    P := arg[2];
    L := LiePRingsByLibrary(d);
    PrintTo(arg[3], "[ \n");
    for i in [1..Length(L)] do
        G := LiePRingsInFamily(L[i], P);
        if not IsBool(G) then 
            for j in [1..Length(G)] do
                G[j] := PGroupByLiePRing(G[j]);
                if not IsBool(G[j]) then 
                    G[j] := CodePcGroup(G[j]);
                    AppendTo(arg[3], G[j],",\n");
                    G[j] := false;
                fi;
            od;
        fi;
        G := false;
    od;
end );

BindGlobal( "Print2GenGroupsViaLiePRings", function( arg )
    local d, P, L, G, i, j;
    d := arg[1];
    P := arg[2];
    L := LiePRingsByLibrary(d);
    PrintTo(arg[3], "[ \n");
    for i in [1..Length(L)] do
        Print("starting ",i," of ",Length(L)," \n");
        if L[i]!.MinimalGeneratorNumberOfLiePRing <= 2 then 
        G := LiePRingsInFamily(L[i], P);
        if not IsBool(G) then
            for j in [1..Length(G)] do
                G[j] := PGroupByLiePRing(G[j]);
                if not IsBool(G[j]) then
                    G[j] := CodePcGroup(G[j]);
                    AppendTo(arg[3], G[j],",\n");
                    G[j] := false;
                fi;
            od;
        fi;
        G := false;
        fi;
    od;
end );

