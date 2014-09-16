#############################################################################
##  
SetPackageInfo( rec(
PackageName := "LiePRing",
Subtitle := "Database and algorithms for Lie p-rings",
Version := "1.6",
Date := "08/11/2013",

PackageWWWHome 
  := "http://www.icm.tu-bs.de/~beick/soft/liepring/",
ArchiveURL 
  := Concatenation( ~.PackageWWWHome, "liepring-", ~.Version ),
ArchiveFormats := ".tar.gz",
Persons := [
  rec( 
    LastName      := "Vaughan-Lee",
    FirstNames    := "Michael",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "michael.vaughan-lee@chch.ox.ac.uk",
    WWWHome       := "users.ox.ac.uk/~vlee",
    place         := "Oxford"),
  rec( 
    LastName      := "Eick",
    FirstNames    := "Bettina",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "beick@tu-bs.de",
    WWWHome       := "www.icm.tu-bs.de/~beick",
    Place         := "TU Braunschweig" ),
],

Status           := "accepted",
CommunicatedBy   := "Leonard Soicher (London)",
AcceptDate       := "09/2014",

README_URL := 
  Concatenation( ~.PackageWWWHome, "README" ),
PackageInfoURL := 
  Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),

AbstractHTML := "",

PackageDoc := rec(
  BookName  := "LiePRing",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "LiePRing Package",
),

AvailabilityTest := ReturnTrue,

Dependencies := rec(
  GAP := "4.5.3",
  NeededOtherPackages := [["LieRing", ">=2.1"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

BannerString := Concatenation( 
    "----------------------------------------------------------------\n",
    "Loading  LiePRing ", ~.Version, "\n",
    "by Michael Vaughan-Lee and Bettina Eick \n",
    "----------------------------------------------------------------\n" ),

Keywords := ["", "", ""]

));

