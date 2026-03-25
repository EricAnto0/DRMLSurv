# Changelog

## DRMLSurv (development version)

### 0.1.0

#### New features

- Added core functions for censoring-score estimation, matching-based
  imputation, and policy learning.
- Added
  [`ComputeScores()`](https://ericanto0.github.io/DRMLSurv/reference/ComputeScores.md),
  [`get_doublescores()`](https://ericanto0.github.io/DRMLSurv/reference/get_doublescores.md),
  and
  [`impute_censored_outcomes()`](https://ericanto0.github.io/DRMLSurv/reference/impute_censored_outcomes.md)
  workflows.
- Added
  [`rfdtr()`](https://ericanto0.github.io/DRMLSurv/reference/rfdtr.md)
  and
  [`Drmatch()`](https://ericanto0.github.io/DRMLSurv/reference/Drmatch.md)
  for dynamic treatment regime learning.

#### Documentation

- Added roxygen2 documentation.
- Added README and package badges.

#### Infrastructure

- Added test suite with `testthat`.
- Added GitHub Actions for R CMD check and coverage.
