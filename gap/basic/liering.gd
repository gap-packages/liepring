
#############################################################################
##
## Introduce the LR (Lie p-ring) elements
##
DeclareCategory( "IsLPRElement", IsRingElement );
DeclareCategoryFamily( "IsLPRElement" );
DeclareCategoryCollections( "IsLPRElement" );
DeclareRepresentation( "IsLPRElementRep",
                        IsComponentObjectRep,
                        ["sctable", "exponents", "name" ] );

DeclareGlobalFunction( "LPRElementConstruction" );
DeclareGlobalFunction( "LPRElementByExponentsNC" );
DeclareGlobalFunction( "LPRElementByExponents" );
DeclareGlobalFunction( "LPRElementByWordList" );

DeclareOperation( "Exponents", [ IsLPRElementRep ] );
DeclareOperation( "SCTable", [ IsLPRElementRep ] );
DeclareOperation( "NameTag", [ IsLPRElementRep ] );

#############################################################################
##
## Introduce LR (Lie p-rings)
##
DeclareGlobalFunction( "CreateLiePRing");
DeclareGlobalFunction( "LiePRingBySCTableNC" );
DeclareGlobalFunction( "LiePRingBySCTable" );
DeclareGlobalFunction( "ViewPCPresentation");
DeclareGlobalFunction( "ViewShortPresentation");
DeclareGlobalFunction( "ParametersOfLiePRing");
DeclareGlobalFunction( "LiePRingsByLibrary");
DeclareGlobalFunction( "LiePRingsDim7ByFile");
DeclareGlobalFunction( "LiePRingsByLibraryMC8");

DeclareProperty( "IsLiePRing", IsRing );
DeclareAttribute( "ClassOfLiePRing", IsLiePRing );
DeclareAttribute( "PClassOfLiePRing", IsLiePRing );
DeclareAttribute( "MinimalGeneratorNumberOfLiePRing", IsLiePRing );
DeclareAttribute( "DimensionOfLiePRing", IsLiePRing );
DeclareAttribute( "PrimeOfLiePRing", IsLiePRing );
DeclareProperty( "IsParentLiePRing", IsLiePRing );
DeclareAttribute( "BasisOfLiePRing", IsLiePRing );

DeclareAttribute( "LibraryConditions", IsLiePRing );
DeclareAttribute( "LibraryName", IsLiePRing );
DeclareAttribute( "LiePValues", IsLiePRing );
DeclareAttribute( "ShortPresentation", IsLiePRing );


