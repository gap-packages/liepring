gap> SC := rec( dim := 3, prime := 2, tab := [] );;
gap> L := LiePRingBySCTable(SC);
<LiePRing of dimension 3 over prime 2>
gap> l := BasisOfLiePRing(L);
[ l1, l2, l3 ]
gap> l[1]*l[2];
0
gap> 2*l[1];
0
gap> l[1] + l[2];
l1 + l2
gap> SC := rec( dim := 4, prime := 5, tab := [ [], [3, 1], [], [4, 1]]);;
gap> L := LiePRingBySCTableNC(SC);;
gap> ViewPCPresentation(L);
[l2,l1] = l3
[l3,l1] = l4
gap> p := IndeterminateByName("p");;
gap> x := IndeterminateByName("x");;
gap> S := rec( dim := 5,
>              param := [ x ],
>              prime := p,
>              tab := [ [ 4, 1 ], [ 3, 1 ], [ 5, x ], [ 4, 1 ], [ 5, 1 ] ] );;
gap> L := LiePRingBySCTable(S);
<LiePRing of dimension 5 over prime p with parameters [ x ]>
gap> ViewPCPresentation(L);
p*l1 = l4
p*l2 = x*l5
[l2,l1] = l3
[l3,l1] = l4
[l3,l2] = l5
gap> l := BasisOfLiePRing(L);
[ l1, l2, l3, l4, l5 ]
gap> p*l[1];
l4
gap> l[1]+l[2];
l1 + l2
gap> l[1]*l[2];
-1*l3
gap> p := IndeterminateByName("p");;
gap> w := IndeterminateByName("w");;
gap> x := IndeterminateByName("x");;
gap> y := IndeterminateByName("y");;
gap> S := rec( dim := 7,
>              param := [ w, x, y ],
>              prime := p,
>              tab := [ [  ], [ 6, 1 ], [ 6, 1 ], [ 7, 1 ], [  ],
>                       [ 6, x, 7, y ], [  ], [ 7, 1 ], [ 6, w ] ] );;
gap> L := LiePRingBySCTable(S);
<LiePRing of dimension 7 over prime p with parameters [ w, x, y ]>
gap> ViewPCPresentation(L);
p*l2 = l6
p*l3 = x*l6 + y*l7
[l2,l1] = l6
[l3,l1] = l7
[l4,l2] = l7
[l4,l3] = w*l6
gap> SpecialiseLiePRing(L, 7, [x, y], [0,0]);
<LiePRing of dimension 7 over prime 7>
gap> ViewPCPresentation(last);
7*l2 = l6
[l2,l1] = l6
[l3,l1] = l7
[l4,l2] = l7
[l4,l3] = 3*l6
gap> SpecialiseLiePRing(L, 11, [x, y], [0,10]);
<LiePRing of dimension 7 over prime 11>
gap> ViewPCPresentation(last);
11*l2 = l6
11*l3 = 10*l7
[l2,l1] = l6
[l3,l1] = l7
[l4,l2] = l7
[l4,l3] = 2*l6
gap> Cartesian([0,1],[0,1]);
[ [ 0, 0 ], [ 0, 1 ], [ 1, 0 ], [ 1, 1 ] ]
gap> List(last, v -> SpecialiseLiePRing(L, 2, [x,y], v));
[ <LiePRing of dimension 7 over prime 2>,
  <LiePRing of dimension 7 over prime 2>,
  <LiePRing of dimension 7 over prime 2>,
  <LiePRing of dimension 7 over prime 2> ]
gap> SpecialiseLiePRing(L, p, [x], [0]);
<LiePRing of dimension 7 over prime p with parameters [ y, w ]>
gap> ViewPCPresentation(last);
p*l2 = l6
p*l3 = y*l7
[l2,l1] = l6
[l3,l1] = l7
[l4,l2] = l7
[l4,l3] = w*l6
gap> SpecialiseLiePRing(L, p, [y], [3]);
<LiePRing of dimension 7 over prime p with parameters [ x, w ]>
gap> ViewPCPresentation(last);
p*l2 = l6
p*l3 = x*l6 + 3*l7
[l2,l1] = l6
[l3,l1] = l7
[l4,l2] = l7
[l4,l3] = w*l6
gap> SpecialisePrimeOfLiePRing(L, 29);
<LiePRing of dimension 7 over prime 29 with parameters [ y, x ]>
gap> ViewPCPresentation(last);
29*l2 = l6
29*l3 = x*l6 + y*l7
[l2,l1] = l6
[l3,l1] = l7
[l4,l2] = l7
[l4,l3] = 2*l6
gap>  L := LiePRingsByLibrary(6)[14];
<LiePRing of dimension 6 over prime p with parameters [ x ]>
gap>  K := SpecialisePrimeOfLiePRing(L, 5);
<LiePRing of dimension 6 over prime 5 with parameters [ x ]>
gap> LiePValues(K);
[ [ p, w ], [ 5, 2 ] ]
gap> L := LiePRingsByLibrary(6)[100];
<LiePRing of dimension 6 over prime p>
gap> l := BasisOfLiePRing(L);
[ l1, l2, l3, l4, l5, l6 ]
gap> U := LiePSubring(L, [5*l[1]]);
<LiePRing of dimension 3 over prime p>
gap> BasisOfLiePRing(U);
[ l1, l4, l6 ]
gap>  K := SpecialisePrimeOfLiePRing(L, 5);
<LiePRing of dimension 6 over prime 5>
gap>  b := BasisOfLiePRing(K);
[ l1, l2, l3, l4, l5, l6 ]
gap> LiePSubring(K, [5*b[1]]);
<LiePRing of dimension 2 over prime 5>
gap>  BasisOfLiePRing(last);
[ l4, l6 ]
gap> K := SpecialisePrimeOfLiePRing(L, 7);
<LiePRing of dimension 6 over prime 7>
gap> b := BasisOfLiePRing(K);
[ l1, l2, l3, l4, l5, l6 ]
gap> U := LiePSubring(K, [5*b[1]]);
<LiePRing of dimension 3 over prime 7>
gap> BasisOfLiePRing(U);
[ l1, l4, l6 ]
gap> LiePIdeal(L, [l[1]]);
<LiePRing of dimension 5 over prime p>
gap> BasisOfLiePRing(last);
[ l1, l3, l4, l5, l6 ]
gap> LiePIdeal(K, [b[1]]);
<LiePRing of dimension 5 over prime 7>
gap> LiePIdeal(K, [b[2]]);
<LiePRing of dimension 4 over prime 7>
gap> LiePQuotient(K,last);
<LiePRing of dimension 2 over prime 7>
