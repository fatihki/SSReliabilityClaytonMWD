## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Resubmission

This is a re-submission. In this version, we addressed the issues raised by CRAN:

- Added missing `\value{}` sections to Rd files, including for the S3 method `print.SSRfit`, to properly document return values and output structure.

- Removed calls to `rm(list = ls())` from vignette files to avoid modifying the user’s global environment.

No other changes were made.