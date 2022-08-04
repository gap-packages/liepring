BindGlobal( "Linearise27a", function(A)
    return IntVecFFE( [ A[4][1], A[4][2], A[4][3], A[3][1], 
                        A[1][1], A[1][2], A[1][3], A[2][1],
                        A[2][2], A[2][3], A[3][2], A[3][3] ] );   
end );

BindGlobal( "ValsFunction27a", function ( P )
    local  F, W, curoots, reps1, reps2, reps3, reps4, reps5, reps6, reps7, 
    reps8, reps9, reps10, reps11, repstable, Btable, Ctable, x, v, y, z, t, 
    u, G, A, B, C, D, a, b, c, d, n, params, l, ind1, ind, index, i, j, new, 
    translft, transrt, urange, zrange;
    F := GF( P );
    W := PrimitiveRootMod( P );
    curoots := List( [ 1 .. P - 1 ], function ( x )
            return 0;
        end );
    for i  in [ 1 .. P - 1 ]  do
        j := i ^ 3 mod P;
        curoots[j] := i;
    od;
    reps1 
     := 
      [ [ [ 0, 0, 0 ], [ 0, 0, 0 ] ], [ [ 0, 0, 0 ], [ 1, 0, 0 ] ], 
          [ [ 0, 0, 0 ], [ 0, 1, 0 ] ], [ [ 1, 0, 0 ], [ 0, 1, 0 ] ], 
          [ [ 1, 0, 0 ], [ 0, 0, 1 ] ], [ [ 0, 0, 0 ], [ 0, 0, 1 ] ], 
          [ [ 0, 0, 0 ], [ 1, 0, 1 ] ], [ [ 0, 0, 0 ], [ W, 0, 1 ] ], 
          [ [ 0, 1, 0 ], [ 0, 0, 1 ] ], [ [ 0, 1, 0 ], [ 1, 0, 1 ] ], 
          [ [ 0, 1, 0 ], [ W, 0, 1 ] ] ] * One( F );
    repstable := [  ];
    Btable := [  ];
    Ctable := [  ];
    for v  in [ 0, 1 ]  do
        for y  in [ 0, 1 ]  do
            for z  in [ 0 .. P - 1 ]  do
                for t  in [ 0 .. P - 1 ]  do
                    for u  in [ 0 .. P - 1 ]  do
                        A := [ [ v, z, t ], [ 0, y, u ] ] * One( F );
                        if EchelonForm( A )[1] = A  then
                            index 
                             := 1 + P ^ 3 * (2 * v + y) + P ^ 2 * z + P * t 
                              + u;
                            G := GetReps1Slow( F, W, A );
                            repstable[index] := G[1];
                            B := G[2];
                            C := G[3];
                            a := B[1][1];
                            b := B[1][2];
                            c := B[2][1];
                            d := B[2][2];
                            n := C[2][1];
                            x := C[2][2];
                            B 
                             := 
                              [ [ a, b, 0, 0 ], [ c, d, 0, 0 ], 
                                  [ c * n, d * n, d * x, (- c) * x ], 
                                  [ (- a) * n, (- b) * n, (- b) * x, a * x ] ]
                              * One( F );
                            Btable[index] := B;
                            Ctable[index] := C;
                        fi;
                    od;
                od;
            od;
        od;
    od;
    #Print( "repstable done \n" );
    reps2 := GetReps2( P );
    reps3 := GetReps3( P );
    reps4 := GetReps4( P );
    reps5 := GetReps5( P );
    reps6 := GetReps6( P );
    reps7 := GetReps7( P );
    reps8 := GetReps8( P );
    reps9 := GetReps9( P );
    #Print( "reps2-9 done \n" );
    translft := [  ];
    transrt := [  ];
    translft[1] 
     := MyCutVector( [ 0, 0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0 ], 4 ) 
      * One( F );
    transrt[1] := MyCutVector( [ 0, 0, 1, 0, 1, 0, 1, 0, 0 ], 3 ) * One( F );
    for i  in [ 1 .. P - 1 ]  do
        Add( translft, 
         MyCutVector( [ 1, 0, 0, (- i), 0, 1, i, 0, 0, 0, 1, 0, 0, 0, 0, 1 ], 4 
             ) * One( F ) );
        Add( transrt, MyCutVector( [ 1, 2 * i, i ^ 2, 0, 1, i, 0, 0, 1 ], 3 ) 
          * One( F ) );
    od;
    #Print( "transversal done \n" );
    params := [  ];
    for i  in [ 1 .. 11 ]  do
        A := NullMat( 4, 3, F );
        A[3] := reps1[i][1];
        A[4] := reps1[i][2];
        Add( params, Linearise27a(A) );
    od;
    #Print( "reps1 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps2 ) ]  do
        A := reps2[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][1] ) + P ^ 4 * IntFFE( A[3][2] ) 
                + P ^ 3 * IntFFE( A[3][3] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 2  then
                new := false;
                break;
            fi;
            if n > 2  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex2( P, W, D );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps2 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps3 ) ]  do
        A := reps3[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][1] ) + P ^ 4 * IntFFE( A[3][2] ) 
                + P ^ 3 * IntFFE( A[3][3] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 3  then
                new := false;
                break;
            fi;
            if n > 3  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex3( P, W, D );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps3 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps4 ) ]  do
        A := reps4[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][1] ) + P ^ 4 * IntFFE( A[3][2] ) 
                + P ^ 3 * IntFFE( A[3][3] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 4  then
                new := false;
                break;
            fi;
            if n > 4  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex4( P, W, D );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps4 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps5 ) ]  do
        A := reps5[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[4][3] ) + P ^ 4 * IntFFE( A[3][3] ) 
                + P ^ 3 * IntFFE( A[4][1] ) + P ^ 2 * IntFFE( A[3][2] ) 
            + P * IntFFE( A[3][1] ) + IntFFE( A[4][2] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 5  then
                new := false;
                break;
            fi;
            if n > 5  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex5( P, W, D, curoots );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps5 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps6 ) ]  do
        A := reps6[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][1] ) + P ^ 4 * IntFFE( A[3][2] ) 
                + P ^ 3 * IntFFE( A[3][3] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 6  then
                new := false;
                break;
            fi;
            if n > 6  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex6( P, W, D, curoots );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps6 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps7 ) ]  do
        A := reps7[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][1] ) + P ^ 4 * IntFFE( A[3][2] ) 
                + P ^ 3 * IntFFE( A[3][3] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 7  then
                new := false;
                break;
            fi;
            if n > 7  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex7( P, W, D );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps7 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps8 ) ]  do
        A := reps8[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][1] ) + P ^ 4 * IntFFE( A[3][2] ) 
                + P ^ 3 * IntFFE( A[3][3] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 8  then
                new := false;
                break;
            fi;
            if n > 8  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex7( P, W, D );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps8 dealt with: ", Length( params ), "\n" );
    for i  in [ 1 .. Length( reps9 ) ]  do
        A := reps9[i];
        new := true;
        index 
         := 
          P ^ 5 * IntFFE( A[3][3] ) + P ^ 4 * IntFFE( A[4][2] ) 
                + P ^ 3 * IntFFE( A[3][2] ) + P ^ 2 * IntFFE( A[4][1] ) 
            + P * IntFFE( A[3][1] ) + IntFFE( A[4][3] );
        for j  in [ 1 .. P ]  do
            B := translft[j];
            C := transrt[j];
            D := B * A * C ^ -1;
            G := GetReps1( P, D, repstable );
            n := G[1];
            ind := G[2];
            B := G[3];
            if n < 9  then
                new := false;
                break;
            fi;
            if n > 9  then
                continue;
            fi;
            a := B[1][1];
            b := B[1][2];
            c := B[2][1];
            d := B[2][2];
            B 
             := 
              MyCutVector( [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, (- c), 0, 0, 
                    (- b), a ], 4 ) * One( F );
            D := B * D;
            B := Btable[ind];
            C := Ctable[ind];
            D := B * D * C ^ -1;
            ind1 := GetIndex9( P, W, D, curoots );
            if ind1 < index  then
                new := false;
                break;
            fi;
        od;
        if new  then
            Add( params, Linearise27a(A) );
        fi;
    od;
    #Print( "reps9 dealt with: ", Length( params ), "\n" );
    for x  in [ 0 .. (P - 1) / 2 ]  do
        zrange := [ 0 .. P - 1 ];
        if x = 0  then
            zrange := [ 0 .. (P - 1) / 2 ];
        fi;
        for y  in [ 0 .. P - 1 ]  do
            for z  in zrange  do
                urange := [ 0 .. P - 1 ];
                if x + z = 0  then
                    urange := [ 0 .. (P - 1) / 2 ];
                fi;
                for t  in [ 0 .. P - 1 ]  do
                    for u  in urange  do
                        for v  in [ 0 .. P - 1 ]  do
                            A 
                             := 
                              MyCutVector( [ 0, 1, 0, 1, 0, 1, x, y, z, t, u, v 
                                   ], 4 ) * One( F );
                            new := true;
                            index 
                             := 
                              P ^ 5 * IntFFE( A[3][1] ) 
                                      + P ^ 4 * IntFFE( A[3][2] ) 
                                    + P ^ 3 * IntFFE( A[3][3] ) 
                                  + P ^ 2 * IntFFE( A[4][1] ) 
                                + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
                            for j  in [ 1 .. P ]  do
                                B := translft[j];
                                C := transrt[j];
                                D := B * A * C ^ -1;
                                G := GetReps1( P, D, repstable );
                                n := G[1];
                                ind := G[2];
                                B := G[3];
                                if n < 10  then
                                    new := false;
                                    break;
                                fi;
                                if n > 10  then
                                    continue;
                                fi;
                                a := B[1][1];
                                b := B[1][2];
                                c := B[2][1];
                                d := B[2][2];
                                B 
                                 := 
                                  MyCutVector( 
                                     [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, 
                                        (- c), 0, 0, (- b), a ], 4 ) 
                                  * One( F );
                                D := B * D;
                                B := Btable[ind];
                                C := Ctable[ind];
                                D := B * D * C ^ -1;
                                ind1 := GetIndex10( P, W, D );
                                if ind1 < index  then
                                    new := false;
                                    break;
                                fi;
                            od;
                            if new  then
                                Add( params, Linearise27a(A) );
                            fi;
                        od;
                    od;
                od;
            od;
        od;
    od;
    #Print( "reps10 dealt with: ", Length( params ), "\n" );
    for x  in [ 0 .. (P - 1) / 2 ]  do
        zrange := [ 0 .. P - 1 ];
        if x = 0  then
            zrange := [ 0 .. (P - 1) / 2 ];
        fi;
        for y  in [ 0 .. P - 1 ]  do
            for z  in zrange  do
                urange := [ 0 .. P - 1 ];
                if x + z = 0  then
                    urange := [ 0 .. (P - 1) / 2 ];
                fi;
                for t  in [ 0 .. P - 1 ]  do
                    for u  in urange  do
                        for v  in [ 0 .. P - 1 ]  do
                            A 
                             := 
                              MyCutVector( [ 0, 1, 0, W, 0, 1, x, y, z, t, u, v 
                                   ], 4 ) * One( F );
                            new := true;
                            index 
                             := 
                              P ^ 5 * IntFFE( A[3][1] ) 
                                      + P ^ 4 * IntFFE( A[3][2] ) 
                                    + P ^ 3 * IntFFE( A[3][3] ) 
                                  + P ^ 2 * IntFFE( A[4][1] ) 
                                + P * IntFFE( A[4][2] ) + IntFFE( A[4][3] );
                            for j  in [ 1 .. P ]  do
                                B := translft[j];
                                C := transrt[j];
                                D := B * A * C ^ -1;
                                G := GetReps1( P, D, repstable );
                                n := G[1];
                                ind := G[2];
                                B := G[3];
                                if n < 11  then
                                    new := false;
                                    break;
                                fi;
                                if n > 11  then
                                    continue;
                                fi;
                                a := B[1][1];
                                b := B[1][2];
                                c := B[2][1];
                                d := B[2][2];
                                B 
                                 := 
                                  MyCutVector( 
                                     [ a, b, 0, 0, c, d, 0, 0, 0, 0, d, 
                                        (- c), 0, 0, (- b), a ], 4 ) 
                                  * One( F );
                                D := B * D;
                                B := Btable[ind];
                                C := Ctable[ind];
                                D := B * D * C ^ -1;
                                ind1 := GetIndex10( P, W, D );
                                if ind1 < index  then
                                    new := false;
                                    break;
                                fi;
                            od;
                            if new  then
                                Add( params, Linearise27a(A) );
                            fi;
                        od;
                    od;
                od;
            od;
        od;
    od;
    #Print( "reps11 dealt with: ", Length( params ), "\n" );
    if P = 3 then 
        l := 550;
    elif P mod 3 = 1  then
        l := P ^ 5 + P ^ 4 + 4 * P ^ 3 + 6 * P ^ 2 + 18 * P + 19;
    else
        l := P ^ 5 + P ^ 4 + 4 * P ^ 3 + 6 * P ^ 2 + 16 * P + 17;
    fi;
    if Length( params ) <> l  then
        Error( "wrong number of params" );
    fi;
    if Length( params ) <> Length( Set( params ) )  then
        Error( "duplicate params" );
    fi;
    return params;
end );
