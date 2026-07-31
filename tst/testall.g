LoadPackage("liepring");
LoadPackage("smallgrp");   # used by tst/p12345.tst
dirs := DirectoriesPackageLibrary("liepring", "tst");
TestDirectory(dirs, rec(exitGAP := true,
        testOptions := rec(compareFunction:="uptowhitespace")));
FORCE_QUIT_GAP(1);
