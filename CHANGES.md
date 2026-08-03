This file describes changes in the LiePRing package.

## 2.9.2 (2026-08-03)

  - Make the SmallGrp package an optional instead of a required dependency
  - Various janitorial changes

## 2.9.1 (2024-06-11)

  - Ensure that unit polynomials have leading coefficient 1

## 2.9 (2024-04-29)

  - Make the package independent of the Singular package by using GAP's
    own Groebner basis code
  - Restore the `DepthVector` helper for backwards compatibility
  - Use GAP's `SizeGL` in place of a package-local implementation
  - Various janitorial changes

## 2.8 (2022-10-21)

  - Make the helpers `MyBaseMat` and `FactorSpace` private
  - Move consistency check code into a separate file

## 2.7 (2022-08-05)

  - Protect global functions from accidental redefinition, and fix
    warnings about `PGroupByLiePRing` being defined twice
  - Replace the package-local `DepthVector` and `Pos` helpers by GAP's
    `PositionNonZero`
  - Correct several minor issues in the manual and update DOI links to
    use https

## 2.6 (2022-04-11)

  - Add a library of Lie p-rings of dimension 8 and maximal class,
    together with functions to access it
  - Split the introductory manual chapter into two chapters and fix
    various documentation issues
  - Add a test file and register it in `PackageInfo.g`
  - Move the test suite from Travis CI to GitHub Actions
  - Add GPL-2.0-or-later license metadata and remove CVS leftovers

## 2.5 (2021-02-02)

  - Reorganize the package sources into subdirectories (`basic`, `rings`,
    `class`, `evals`, `advan`, `check`)
  - Add `AutGrpDescription` for describing automorphism groups of generic
    Lie p-rings
  - Add functions for polynomial arithmetic, units and zeros, and for
    computing series of Lie p-rings
  - Add a new manual chapter on Lie p-rings

## 2.0 (2019-10-11)

  - Add an updated library of Lie p-rings of dimension at most 7 and
    adjust the code to the changed database
  - Add `LiePSchurMult` for computing the Schur multipliers of the
    Lie p-rings in a family, plus `ElementNumber` and `ElementNumbers`
    for counting the parameters realizing them
  - Add `FindAutos` for describing the automorphism group of a generic
    Lie p-ring
  - Add supporting code for ring theory, number theory and linear algebra
  - Fix the citation of Khu98

## 1.9.2 (2018-10-11)

  - Add further tests

## 1.9.1 (2018-10-07)

  - Change `NumberOfLiePRingsInFamilyNC` to return `fail` as a last resort
  - Simplify `EchelonForm` and the use of `QuotRemInt`
  - Add many more tests, and declare the existing test suite in
    `PackageInfo.g`
  - Convert the README to markdown

## 1.9 (2018-03-11)

  - Modernize the package for newer GAP versions by no longer using
    `RequirePackage` and `MutableNullMat`
  - Turn `IsLiePRing` and `IsParentLiePRing` into properties
  - Fix the HTML version of the manual and tweak the manual build
  - Add tests based on the manual examples
  - Add continuous integration and code-coverage support
  - Move the package to its current GitHub repository

## 1.8 (2014-12-02)

  - Extend `LiePRingsInFamily` by an optional flag to return the
    corresponding p-groups or their codes instead of Lie p-rings

## 1.6 (2014-09-16)

  - Add code for computing PORC polynomials counting the Lie p-rings in
    a family
  - Support the additional parameters `(p-1,n)`, `(p+1,3)` and
    `(p^2-1,16)` used by the library
  - Add consistency checks and revise the manual

## 1.5 (2013-11-28)

  - Correct the tables of evaluated values used by the library

## 1.4 (2013-11-25)

  - Fix the path to the start page of the HTML manual

## 1.3 (2013-11-25)

  - Fix the name of the start page of the HTML manual

## 1.2 (2013-11-25)

  - Fix the start page of the HTML manual

## 1.1 (2013-11-25)

  - Add a README and fix the manual build script

## 1.0 (2013-11-13)

  - Initial release
