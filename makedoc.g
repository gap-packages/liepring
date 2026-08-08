#############################################################################
##
##  makedoc.g
##
##  Builds the package documentation with AutoDoc/GAPDoc.
##
#############################################################################

LoadPackage("AutoDoc");

# Run this from the package's root directory: gap makedoc.g
AutoDoc(rec(
    autodoc := rec(scan_dirs := []),
    gapdoc := rec(main := "main", files := []),
    extract_examples := true,
    scaffold := rec(
        includes := [
            "preamble.xml",
            "intro.xml",
            "lierings.xml",
            "database.xml",
            "advan.xml"
        ],
        entities := rec(
            LiePRing := "<Package>LiePRing</Package>",
        ),
        bib := "liepring.bib",
    ),
));

QuitGap();
