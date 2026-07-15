## Resubmission

This is a resubmission (version 2.2.1) addressing the points raised by
CRAN review on the 2.2.0 submission:

* Added \value tags to the .Rd files for all exported methods
  (matrix_i2, ml_cvh1, ml_cvh2, ml_m1, ml_m2, pred_interval_esmeans,
  R2_calc, ratio_i2), documenting the structure (class) and meaning of
  the returned objects.
* Removed the example for the unexported function transform_data()
  (the function remains an internal helper).
* Replaced all \dontrun{} wrappers with \donttest{} and made the examples
  executable (the affected examples fit small metafor models and run in a
  few seconds each).

## Test environments

* local macOS, R 4.5.3
* win-builder (R-devel and R-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* The single NOTE is "New submission", which is expected for a first-time
  submission to CRAN.

## Notes for the maintainer

* The DOI in inst/CITATION (<doi:10.1111/2041-210X.14152>) is the published
  Methods in Ecology and Evolution article. It is valid; some automated
  checkers may report a 403 because the publisher blocks non-browser requests.
