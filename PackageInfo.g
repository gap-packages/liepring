#############################################################################
##  
SetPackageInfo( rec(
PackageName := "LiePRing",
Subtitle := "Database and algorithms for Lie p-rings",
Version := "2.9.1",
Date := "11/06/2024",
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
    place         := "Oxford"),
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

AbstractHTML := "",

PackageDoc := rec(
  BookName  := "LiePRing",
  ArchiveURLSubset := ["doc","htm"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "LiePRing Package",
),

AvailabilityTest := ReturnTrue,

Dependencies := rec(
  GAP := "4.8",
  NeededOtherPackages := [["LieRing", ">=2.1"]],
  SuggestedOtherPackages := [["Singular", ">=10"]],
  ExternalConditions := []
),

BannerString := Concatenation( 
    "----------------------------------------------------------------\n",
    "Loading  LiePRing ", ~.Version, "\n",
    "by Bettina Eick and Michael Vaughan-Lee \n",
    "----------------------------------------------------------------\n" ),

Keywords := ["", "", ""],

TestFile := "tst/testall.g" 

));

