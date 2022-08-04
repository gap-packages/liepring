
BindGlobal( "EchelonForm", function(M)
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
      
        j := First([c..n], x -> A[x][i] <> Zero(F));
        if not IsBool(j) then

            # swop rows
            if j <> c then 
                C := IdentityMat(n, F);
                r := C[j];
                C[j] := C[c];
                C[c] := r;
                r := A[j];
                A[j] := A[c];
                A[c] := r;
                B := C * B;
            fi;

            # norm privot
            if A[c][i] <> One(F) then 
                C := IdentityMat(n, F);
                C[c] := C[c] / A[c][i];
                A[c] := A[c] / A[c][i];
                B := C * B;
            fi;

            # clear column
            for k in [1..n] do
                if k <> c and A[k][i] <> Zero(F) then 
                    C := IdentityMat(n, F);
                    C[k] := C[k] - A[k][i] * C[c];
                    A[k] := A[k] - A[k][i] * A[c];
                    B := C * B;
                fi;
            od;
 
            # iterate
            c := c+1;
        fi;
    od;
   
    if B*M <> A then Error("echelon form wrong"); fi;
    return [A,B];
end );
                
                     
