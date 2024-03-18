#############################################################################
##  
SetPackageInfo( rec(
PackageName := "LiePRing",
Subtitle := "Database and algorithms for Lie p-rings",
Version := "2.8",
Date := "21/10/2022",
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
  ArchiveURLSubset := ["doc"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "LiePRing Package",
),

AvailabilityTest := function()
  # When LiePRing gets loaded, the Singular package will be loaded,
  # and then 'StartSingular()' will be called in read.g,
  # without modifying the variable 'sing_exec'.
  # Thus we use the code from Singular's
  # 'CheckSingularExecutableAndTempDir' in order to determine
  # whether a Singular executable will be found by the Singular package.
  # If not then LiePRing cannot be loaded,
  # because otherwise the call of 'StartSingular()' would run into an error.
  local IsExec, sing_exec;

  IsExec:= path -> IsString( path ) and not IsDirectoryPath( path )
                   and IsExecutableFile( path );

  if IsBoundGlobal( "sing_exec" ) then
    # The Singular package has been loaded.
    sing_exec:= ValueGlobal( "sing_exec" );
  else
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
  return IsExec( sing_exec );
end,

Dependencies := rec(
  GAP := "4.8",
  NeededOtherPackages := [["LieRing", ">=2.1"], ["Singular", ">=10"]],
  SuggestedOtherPackages := [],
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

