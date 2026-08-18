#############################################################################
##  
SetPackageInfo( rec(
PackageName := "LiePRing",
Subtitle := "Database and algorithms for Lie p-rings",
Version := "2.9.3",
Date := "18/08/2026",
License := "GPL-2.0-or-later",

Persons := [
  rec( 
    LastName      := "Eick",
    FirstNames    := "Bettina",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "beick@tu-bs.de",
    WWWHome       := "http://www.iaa.tu-bs.de/beick",
    Place         := "TU Braunschweig" ),
  rec( 
    LastName      := "Vaughan-Lee",
    FirstNames    := "Michael",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "michael.vaughan-lee@chch.ox.ac.uk",
    WWWHome       := "http://users.ox.ac.uk/~vlee",
    Place         := "Oxford"),
],

Status           := "accepted",
CommunicatedBy   := "Leonard Soicher (London)",
AcceptDate       := "09/2014",

PackageWWWHome  := "https://gap-packages.github.io/liepring/",
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/liepring",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/liepring-", ~.Version ),
ArchiveFormats := ".tar.gz",

AbstractHTML := "The <span class=\"pkgname\">LiePRing</span> package provides \
access to a database of nilpotent Lie rings of order p^n for p > 2 and n &lt;= 7, \
together with algorithms for computing with them.",

PackageDoc := rec(
  BookName  := "LiePRing",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "LiePRing Package",
),

AvailabilityTest := ReturnTrue,

Dependencies := rec(
  GAP := "4.8",
  NeededOtherPackages := [["LieRing", ">=2.1"]],
  SuggestedOtherPackages := [["Singular", ">=10"]],
  TestPackages := [["SmallGrp", ">=1.0"]],
  ExternalConditions := []
),

BannerString := Concatenation( 
    "----------------------------------------------------------------\n",
    "Loading  LiePRing ", ~.Version, "\n",
    "by Bettina Eick and Michael Vaughan-Lee \n",
    "----------------------------------------------------------------\n" ),

Keywords := ["Lie ring", "Lie p-ring", "p-group", "classification"],

TestFile := "tst/testall.g",

AutoDoc := rec(
  TitlePage := rec(
    Abstract := [
      "&LiePRing; gives access to the database of Lie <M>p</M>-rings of order",
      "at most <M>p^7</M> as determined by Mike Newman, Eamonn O'Brien and",
      "Michael Vaughan-Lee, see <Cite Key=\"NOV04\"/> and <Cite Key=\"OVL05\"/>,",
      "and it provides some functionality to work with these Lie <M>p</M>-rings.",
      "<P/>",
      "If you use &LiePRing;, then please cite it as:",
      "<E>Bettina Eick and Michael Vaughan-Lee</E>, LiePRing -- A GAP Package",
      "for computing with nilpotent Lie rings of prime-power order (2014), see",
      "<URL>https://www.gap-system.org/Packages/liepring.html</URL>",
    ],
    Copyright := [
      "&LiePRing; is free software; you can redistribute it under the terms of",
      "the <URL Text=\"GNU General Public License\">https://www.fsf.org/licenses/gpl.html</URL>",
      "as published by the Free Software Foundation; either version 2 of the",
      "License, or (at your option) any later version. &LiePRing; is distributed",
      "in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even",
      "the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR",
      "PURPOSE. See the GNU General Public License for more details.",
    ],
    Acknowledgements := [
      "The Lazard correspondence induces a one-to-one correspondence between the",
      "Lie <M>p</M>-rings of order <M>p^n</M> and class less than <M>p</M> and",
      "the <M>p</M>-groups of order <M>p^n</M> and class less than <M>p</M>.",
      "&LiePRing; provides a function to evaluate this correspondence; this",
      "function has been implemented and given to us by Willem de Graaf.",
    ],
  ),
),

));

