#############################################################################
##  
SetPackageInfo( rec(
PackageName := "LiePRing",
Subtitle := "Database and algorithms for Lie p-rings",
Version := "1.9.2",
Date := "11/10/2018",

Persons := [
  rec( 
    LastName      := "Vaughan-Lee",
    FirstNames    := "Michael",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "michael.vaughan-lee@chch.ox.ac.uk",
    WWWHome       := "http://users.ox.ac.uk/~vlee",
    Place         := "Oxford",
    Institution   := "Oxford University"),
  rec( 
    LastName      := "Eick",
    FirstNames    := "Bettina",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "beick@tu-bs.de",
    WWWHome       := "http://www.icm.tu-bs.de/~beick",
    Place         := "Braunschweig",
    Institution   := "TU Braunschweig"),
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
  ArchiveURLSubset := ["doc"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "LiePRing Package",
),

AvailabilityTest := ReturnTrue,

Dependencies := rec(
  GAP := "4.8",
  NeededOtherPackages := [["LieRing", ">=2.1"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

BannerString := Concatenation( 
    "----------------------------------------------------------------\n",
    "Loading  LiePRing ", ~.Version, "\n",
    "by Michael Vaughan-Lee and Bettina Eick \n",
    "----------------------------------------------------------------\n" ),

TestFile := "tst/testall.g",

Keywords := ["", "", ""]

));

