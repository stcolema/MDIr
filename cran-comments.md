## Resubmission

This is a resubmission of a package archived on 2023-05-31 for check problems
(installation failures on macOS and Fedora/clang, and a NOTE that GNU make is a
SystemRequirements). Both came from the build system and are removed:

* `std::execution::par` (C++17 parallel algorithms, unsupported by libc++ on
  macOS), RcppParallel, TBB and OpenMP have been removed. All loops are serial.
* `src/Makevars` and `src/Makevars.win` no longer call `Rscript` through
  `$(shell ...)`, and `SystemRequirements: GNU make` has been dropped.
* C++17 is no longer forced.

## Test environments

* Ubuntu 24.04, R 4.3.3, gcc 13: `R CMD check --as-cran` (examples with
  `--run-donttest`, tests and vignette included)
* Ubuntu 24.04, clang 18 with libc++: package installs and runs
  (an approximation of the macOS toolchain, not a macOS test)

## R CMD check results

0 errors | 0 warnings | 1 note (in a normal locale)

* NOTE: "New submission / Package was archived on CRAN": see above.

Not yet tested on macOS, Windows or R-devel; please check with win-builder and
`rhub` before submission.
