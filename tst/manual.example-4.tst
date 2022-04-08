gap> LiePRingsByLibrary(4);
[ <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p>,
  <LiePRing of dimension 4 over prime p> ]
gap> LiePRingsByLibrary(3, 5);
[ <LiePRing of dimension 3 over prime 5>,
  <LiePRing of dimension 3 over prime 5>,
  <LiePRing of dimension 3 over prime 5>,
  <LiePRing of dimension 3 over prime 5>,
  <LiePRing of dimension 3 over prime 5> ]
gap> LiePRingsByLibrary(5, 2, 4);
[ <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p> ]
gap> LiePRingsByLibrary(5, 7, 2, 4);
[ <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7>,
  <LiePRing of dimension 5 over prime 7> ]
gap> List([1..7], x -> NumberOfLiePRings(x));
[ 1, 2, 5, 15, 75, 542, 4773 ]
gap> L := LiePRingsByLibrary(7)[780];
<LiePRing of dimension 7 over prime p with parameters
[ x, y, z, t, s, u, v ]>
gap> NumberOfLiePRingsInFamily(L);
-1/3*p^5*(p-1,3)+p^5-1/3*p^4*(p-1,3)+p^4-1/3*p^3*(p-1,3)+p^3-1/3*p^2*(p-1,3)
+p^2-p*(p-1,3)+3*p-3/2*(p-1,3)+9/2
gap> L := LiePRingsByLibrary(7)[118];
<LiePRing of dimension 7 over prime p with parameters [ x, y ]>
gap> LibraryConditions(L);
[ "[x,y]~[x,-y]", "p=1 mod 4" ]
gap> LiePRingsInFamily(L, 7);
fail
gap> Length(LiePRingsInFamily(L,13));
91
gap> 13^2;
169
gap> L := LiePRingsByLibrary(5);;
gap> L := Filtered(L, x -> PClassOfLiePRing(x)=4);
[ <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p>,
  <LiePRing of dimension 5 over prime p> ]
gap> K := List(L, x-> LiePRingsInFamily(x, 29));
[ [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ], fail, fail,
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ], fail, fail,
  [ <LiePRing of dimension 5 over prime 29> ],
  [ <LiePRing of dimension 5 over prime 29> ] ]
gap> K := Filtered(Flat(K), x -> x<>fail);
[ <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29>,
  <LiePRing of dimension 5 over prime 29> ]
gap> L := LiePRingsByLibrary(7)[118];
<LiePRing of dimension 7 over prime p with parameters [ x, y ]>
gap> LibraryName(L);
"7.118"
gap> LibraryConditions(L);
[ "[x,y]~[x,-y]", "p=1 mod 4" ]
gap> L := LiePRingsByLibrary(7)[118];
<LiePRing of dimension 7 over prime p with parameters [ x, y ]>
gap> K := SpecialiseLiePRing(L, 13, ParametersOfLiePRing(L), [0,0]);
<LiePRing of dimension 7 over prime 13>
gap> LibraryName(K);
"7.118"
gap> LibraryConditions(K);
[ "[x,y]~[x,-y]", "p=1 mod 4" ]
gap> L := LiePRingsByLibrary(7);;
gap> Filtered(L, x -> LibraryName(x) = "7.1010")[1];
<LiePRing of dimension 7 over prime p>
gap> LIE_TABLE[100];
[ "3gen/gapdec6.139", 1/2*p+(p-1,3)+3/2 ]
gap> LiePRingsDim7ByFile(100);
[ <LiePRing of dimension 7 over prime p>,
  <LiePRing of dimension 7 over prime p>,
  <LiePRing of dimension 7 over prime p>,
  <LiePRing of dimension 7 over prime p>,
  <LiePRing of dimension 7 over prime p with parameters [ x ]> ]
gap> LiePRingsDim7ByFile(100, 7);
[ <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7>,
  <LiePRing of dimension 7 over prime 7> ]
