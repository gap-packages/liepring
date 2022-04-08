gap> LL := LiePRingsByLibrary(7);;
gap> L := Filtered(LL, x -> Length(ParametersOfLiePRing(x))=2)[1];
<LiePRing of dimension 7 over prime p with parameters [ x, y ]>
gap> NumberOfLiePRingsInFamily(L);
p^2-p
gap> RingInvariants(L);
rec( units := [ x ], zeros := [  ] )
gap> ss := LiePSchurMult(L);
[ rec( norm := [ p ], units := [ x, y ], zeros := [ x*y^2-x*y+1 ] ),
  rec( norm := [ p^2 ], units := [ x ], zeros := [ x*y ] ),
  rec( norm := [ p ], units := [ x, x*y^2-x*y+1, y ], zeros := [  ] ) ]
gap> ElementNumbers(ParametersOfLiePRing(L), ss);
rec( norms := [ [ p^2 ], [ p ] ], numbs := [ p-1, p^2-2*p+1 ] )
gap> L := Filtered(LL, x -> Length(ParametersOfLiePRing(x))=2)[1];
<LiePRing of dimension 7 over prime p with parameters [ x, y ]>
gap> AutGroupDescription(L);
rec( auto := [ [ 1, 0, A13, A14, A15, A16, A17 ],
               [ 0, 1, A23, A24, A25, A26, A27 ] ],
     eqns := [ [  ], [  ] ] )
gap> L := Filtered(LL, x -> Length(ParametersOfLiePRing(x))=2)[2];
<LiePRing of dimension 7 over prime p with parameters [ x, y ]>
gap> AutGroupDescription(L);
rec( auto := [ [ A22^3, 0, A13, A14, A15, A16, A17 ],
               [ 0, A22, A23, A24, A25, A26, A27 ] ],
     eqns := [ [ A22*A24-1/2*A23^2, A22^2*y-y,
                 A22*A23^2*y-2*A24*y, A22^4-1,
                 A23^4*y-4*A24^2*y, A22^3*A23^2-2*A24,
                 A22^2*A23^4-4*A24^2, A22*A23^6-8*A24^3,
                 A23^8-16*A24^4 ] ] )
gap> L := LiePRingsByLibrary(7)[489];
<LiePRing of dimension 7 over prime p with parameters [ x ]>
gap> AutGroupDescription(L);
[ rec( auto := [ [ 1, 0, A13, A14, A15, A16, A17 ],
                 [ 0, 1, A23, A24, A25, A26, A27 ] ],
       comment := "p^8 automorphisms",
       eqns := [ [ A13^2*x-A13*A23+2*A15*x+A14-A25,
              -A13*A23*x+A14*x+A23^2-A25*x-2*A24 ] ] ),
  rec( auto := [ [ 0, A12, A13, A14, A15, A16, A17 ],
                 [ -x, 0, A23, A24, A25, A26, A27 ] ],
      comment := "p^8 automorphisms when x <> 0 mod p",
      eqns := [ [ A12^2*A24*x-A12*A13*A23*x+A12*A13*x^2
                  +2*A12*A15*x^2+A12*A14*x-A13^2*x+A13*x+A15*x-A14,
                  -A12^2*A23*x^3+A12*A13*x^3+A12*A23^2*x-A12*A25*x^2
                  -2*A12*A24*x+A13*A23*x+A13*x^2-A15*x^2+A23*x+A25*x-A24 ],
                [ A12*x+1 ] ] ) ]
