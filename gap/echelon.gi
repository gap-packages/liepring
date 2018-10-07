EchelonForm := function(M)
    local A, n, m, F, B, c, i, j, C, r, k;

    # copy, get dims and field and get rid off some trivial cases
    A := List(M, ShallowCopy);
    n := Length(A);
    if n = 0 then return fail; fi;
    m := Length(A[1]);
    if m = 0 then return fail; fi;
    F := Field(A[1][1]);
    B := IdentityMat(n, F);
    if A = 0*A then return [A,B]; fi;

    c := 1;
    for i in [1..m] do
      
        j := First([c..n], x -> not IsZero(A[x][i]));
        if j = fail then
            continue;
        fi;

        # swap rows
        if j <> c then 
            B{[j,c]} := B{[c,j]};
            A{[j,c]} := A{[c,j]};
        fi;

        # norm privot
        if A[c][i] <> One(F) then 
            B[c] := B[c] / A[c][i];
            A[c] := A[c] / A[c][i];
        fi;

        # clear column
        for k in [1..n] do
            if k <> c and not IsZero(A[k][i]) then
                B[k] := B[k] - A[k][i] * B[c];
                A[k] := A[k] - A[k][i] * A[c];
            fi;
        od;

        # go to next column
        c := c+1;
    od;
   
    if B*M <> A then Error("echelon form wrong"); fi;
    return [A,B];
end;
