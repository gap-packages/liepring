#############################################################################
##
## variables and polynomials are a problem - GAP is very slow in computing
## with polynomials and Singular is fast, but has bugs.
##

IsSingularAvailable_known:= false;

if not IsBound( StartSingular ) then
  # Avoid warnings about unbound global variables.
  StartSingular:= "dummy";
fi;

BindGlobal( "StartSingularIfAvailable", function()
  local IsExec, sing_exec;

  # Once we got a positive answer, we need not check.
  if IsSingularAvailable_known = true then
    return true;
  fi;

  # If the Singular package is not available
  # then we cannot use the Singular executable.
  if not IsPackageMarkedForLoading( "Singular", "" ) then
    return false;
  fi;

  # Use the code from Singular's
  # 'CheckSingularExecutableAndTempDir' in order to determine
  # whether a Singular executable will be found by the Singular package.
  # Note that the fact that the Singular package had been loaded
  # does not imply that a Singular executable will be found.
  IsExec:= path -> IsString( path ) and not IsDirectoryPath( path )
                   and IsExecutableFile( path );

  if IsBoundGlobal( "sing_exec" ) then
    # The Singular package has been loaded.
    sing_exec:= ValueGlobal( "sing_exec" );
    if not IsString( sing_exec ) then
      # This may happen after a failed 'StartSingular()' call.
      sing_exec:= "singular";
    fi;
  else
    # Use the default path.
    sing_exec:= "singular";
  fi;
  if IsDirectoryPath( sing_exec ) then
    sing_exec:= Filename( Directory( sing_exec ), "Singular" );
  elif not IsExecutableFile( sing_exec ) then
    sing_exec:= Filename( DirectoriesSystemPrograms(), sing_exec );
  fi;
  if not IsExec( sing_exec ) then
    sing_exec:= Filename( DirectoriesSystemPrograms(), "Singular" );
  fi;
  if IsExec( sing_exec ) then
    IsSingularAvailable_known:= true;
    StartSingular();
    return true;
  fi;
  return false;
end );

if IsString( StartSingular ) then
  Unbind( StartSingular );
fi;

StartSingularIfAvailable();
ORDER := MonomialGrlexOrdering();

if not IsBound( SINGULARGBASIS ) then
  SINGULARGBASIS:= "dummy";
fi;

BindGlobal( "CallGroebner", function( list )
    local oldvar, new;
    if Length(list) = 0 then return list; fi;
    if StartSingularIfAvailable() then
      # Use the Singular standalone.
      oldvar := GBASIS;
      GBASIS := SINGULARGBASIS;
      new := GroebnerBasis(list, ORDER);
      GBASIS := oldvar;
    else
      # Use GAP's implementation.
      new := GroebnerBasis(list, ORDER);
    fi;
    return ReducedGroebnerBasis(new, ORDER);
end );

if IsString( SINGULARGBASIS ) then
  Unbind( SINGULARGBASIS );
fi;

USE_GAP_FACS := false;

if not IsBound( SingularSetBaseRing ) then
  # Avoid warnings about unbound global variables.
  SingularSetBaseRing:= "dummy";
  FactorsUsingSingularNC:= "dummy";
fi;

BindGlobal( "MyFactors", function(pp, h)
    local R, k, w, q;

    if USE_GAP_FACS then return Factors(h); fi;
    if not StartSingularIfAvailable() then return Factors(h); fi;

    # singular is faster, but it has bugs ...
    w := IndeterminateByName("w");
    q := Concatenation(pp, [w]);
    R := PolynomialRing(Rationals, q);
    SingularSetBaseRing(R);
    k := FactorsUsingSingularNC(h);
    if Product(k) = h then return k{[2..Length(k)]}; fi;

    # fall back
    return Factors(h);
end );

if IsString( SingularSetBaseRing ) then
  Unbind( SingularSetBaseRing );
  Unbind( FactorsUsingSingularNC );
fi;
