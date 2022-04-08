gap> START_TEST("p12345.tst");

#
gap> ReadPackage("liepring", "tst/testutils.g");
true

#
gap> Ls := LiePRingsByLibrary(2);
[ <LiePRing of dimension 2 over prime p>, 
  <LiePRing of dimension 2 over prime p> ]
gap> List(Ls, LibraryName);
[ "2.0", "2.1" ]
gap> PrintConditions(Ls);
gap> CountFamilies(Ls);
1: 0, 1, 1, 1, 1, 1, 
2: 0, 1, 1, 1, 1, 1, 

#
gap> Ls := LiePRingsByLibrary(3);
[ <LiePRing of dimension 3 over prime p>, 
  <LiePRing of dimension 3 over prime p>, 
  <LiePRing of dimension 3 over prime p>, 
  <LiePRing of dimension 3 over prime p>, 
  <LiePRing of dimension 3 over prime p> ]
gap> List(Ls, LibraryName);
[ "3.0", "3.2", "3.3", "3.4", "3.1" ]
gap> PrintConditions(Ls);
gap> CountFamilies(Ls);
1: 0, 0, 1, 1, 1, 1, 
2: 0, 1, 1, 1, 1, 1, 
3: 0, 1, 1, 1, 1, 1, 
4: 0, 1, 1, 1, 1, 1, 
5: 0, 1, 1, 1, 1, 1, 

#
gap> Ls := LiePRingsByLibrary(4);
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
gap> List(Ls, LibraryName);
[ "4.0", "4.6", "4.7", "4.8", "4.9", "4.10", "4.11", "4.12", "4.13", "4.14", 
  "4.2", "4.3", "4.4", "4.5", "4.1" ]
gap> PrintConditions(Ls);
gap> CountFamilies(Ls);
1: 0, 0, 1, 1, 1, 1, 
2: 0, 1, 1, 1, 1, 1, 
3: 0, 1, 1, 1, 1, 1, 
4: 0, 1, 1, 1, 1, 1, 
5: 0, 0, 1, 1, 1, 1, 
6: 0, 0, 1, 1, 1, 1, 
7: 0, 0, 1, 1, 1, 1, 
8: 0, 0, 1, 1, 1, 1, 
9: 0, 0, 1, 1, 1, 1, 
10: 0, 0, 1, 1, 1, 1, 
11: 0, 1, 1, 1, 1, 1, 
12: 0, 1, 1, 1, 1, 1, 
13: 0, 1, 1, 1, 1, 1, 
14: 0, 1, 1, 1, 1, 1, 
15: 0, 1, 1, 1, 1, 1, 

#
gap> Ls := LiePRingsByLibrary(5);
[ <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p with parameters [ x ]>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p with parameters [ x ]>, 
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
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p with parameters [ x ]>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p with parameters [ x ]>, 
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
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p>, 
  <LiePRing of dimension 5 over prime p> ]
gap> List(Ls, LibraryName);
[ "5.0", "5.37", "5.38", "5.39", "5.40", "5.41", "5.42", "5.43", "5.44", 
  "5.45", "5.46", "5.47", "5.48", "5.49", "5.50", "5.51", "5.52", "5.53", 
  "5.54", "5.55", "5.56", "5.57", "5.58", "5.59", "5.60", "5.61", "5.62", 
  "5.63", "5.64", "5.65", "5.66", "5.67", "5.68", "5.69", "5.70", "5.71", 
  "5.72", "5.73", "5.74", "5.8", "5.9", "5.10", "5.11", "5.12", "5.13", 
  "5.14", "5.15", "5.16", "5.17", "5.18", "5.19", "5.20", "5.21", "5.22", 
  "5.23", "5.24", "5.25", "5.26", "5.27", "5.28", "5.29", "5.30", "5.31", 
  "5.32", "5.33", "5.34", "5.35", "5.36", "5.2", "5.3", "5.4", "5.5", "5.6", 
  "5.7", "5.1" ]
gap> PrintConditions(Ls);
7: [ "x ne 0, x~x^-1", "" ]
11: [ "1+4x not a square", "" ]
28: [ "", "p=1 mod 3" ]
29: [ "", "p=1 mod 3" ]
33: [ "", "p=1 mod 4" ]
34: [ "", "p=1 mod 4" ]
36: [ "", "p=1 mod 3" ]
37: [ "", "p=1 mod 3" ]
52: [ "x ne 0, x~x^-1", "" ]
55: [ "1+4x not a square", "" ]
gap> CountFamilies(Ls);
1: 0, 0, 0, 1, 1, 1, 
2: 0, 1, 1, 1, 1, 1, 
3: 0, 0, 1, 1, 1, 1, 
4: 0, 0, 1, 1, 1, 1, 
5: 0, 0, 1, 1, 1, 1, 
6: 0, 0, 1, 1, 1, 1, 
7: 0, 0, 3, 4, 6, 7, 
8: 0, 0, 1, 1, 1, 1, 
9: 0, 0, 1, 1, 1, 1, 
10: 0, 0, 1, 1, 1, 1, 
11: 0, 0, 2, 3, 5, 6, 
12: 0, 0, 1, 1, 1, 1, 
13: 0, 0, 1, 1, 1, 1, 
14: 0, 0, 1, 1, 1, 1, 
15: 0, 0, 1, 1, 1, 1, 
16: 0, 0, 1, 1, 1, 1, 
17: 0, 0, 1, 1, 1, 1, 
18: 0, 0, 1, 1, 1, 1, 
19: 0, 0, 1, 1, 1, 1, 
20: 0, 0, 1, 1, 1, 1, 
21: 0, 0, 1, 1, 1, 1, 
22: 0, 0, 1, 1, 1, 1, 
23: 0, 0, 1, 1, 1, 1, 
24: 0, 0, 1, 1, 1, 1, 
25: 0, 0, 1, 1, 1, 1, 
26: 0, 0, 1, 1, 1, 1, 
27: 0, 0, 1, 1, 1, 1, 
28: 0, 0, 0, 1, 0, 1, 
29: 0, 0, 0, 1, 0, 1, 
30: 0, 0, 1, 1, 1, 1, 
31: 0, 0, 1, 1, 1, 1, 
32: 0, 0, 1, 1, 1, 1, 
33: 0, 0, 1, 0, 0, 1, 
34: 0, 0, 1, 0, 0, 1, 
35: 0, 0, 1, 1, 1, 1, 
36: 0, 0, 0, 1, 0, 1, 
37: 0, 0, 0, 1, 0, 1, 
38: 0, 0, 1, 1, 1, 1, 
39: 0, 0, 1, 1, 1, 1, 
40: 0, 1, 1, 1, 1, 1, 
41: 0, 1, 1, 1, 1, 1, 
42: 0, 1, 1, 1, 1, 1, 
43: 0, 1, 1, 1, 1, 1, 
44: 0, 1, 1, 1, 1, 1, 
45: 0, 1, 1, 1, 1, 1, 
46: 0, 1, 1, 1, 1, 1, 
47: 0, 1, 1, 1, 1, 1, 
48: 0, 1, 1, 1, 1, 1, 
49: 0, 1, 1, 1, 1, 1, 
50: 0, 1, 1, 1, 1, 1, 
51: 0, 1, 1, 1, 1, 1, 
52: 0, 2, 3, 4, 6, 7, 
53: 0, 1, 1, 1, 1, 1, 
54: 0, 1, 1, 1, 1, 1, 
55: 0, 1, 2, 3, 5, 6, 
56: 0, 0, 1, 1, 1, 1, 
57: 0, 0, 1, 1, 1, 1, 
58: 0, 0, 1, 1, 1, 1, 
59: 0, 0, 1, 1, 1, 1, 
60: 0, 0, 1, 1, 1, 1, 
61: 0, 0, 1, 1, 1, 1, 
62: 0, 0, 1, 1, 1, 1, 
63: 0, 0, 1, 1, 1, 1, 
64: 0, 0, 1, 1, 1, 1, 
65: 0, 0, 1, 1, 1, 1, 
66: 0, 0, 1, 1, 1, 1, 
67: 0, 0, 1, 1, 1, 1, 
68: 0, 0, 1, 1, 1, 1, 
69: 0, 1, 1, 1, 1, 1, 
70: 0, 1, 1, 1, 1, 1, 
71: 0, 1, 1, 1, 1, 1, 
72: 0, 1, 1, 1, 1, 1, 
73: 0, 1, 1, 1, 1, 1, 
74: 0, 1, 1, 1, 1, 1, 
75: 0, 1, 1, 1, 1, 1, 

#
gap> grps := GroupsViaLiePRings(5,7);;
gap> Length(grps) = NrSmallGroups(7^5);
true

#
gap> STOP_TEST("p12345.tst", 1);
