
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DRMLSurv

<!-- badges: start -->

[![R-CMD-check](https://github.com/EricAnto0/DRMLSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EricAnto0/DRMLSurv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

DRMLSurv is a two-stage survival analysis toolkit that imputes censored
stage-specific times by donor matching, builds matched counterfactual
outcomes for static and dynamic treatment regimes, and learns optimal
dynamic treatment rules with machine learning, cross-fitting, and
balance diagnostics.

## Installation

You can install the development version of DRMLSurv from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("EricAnto0/DRMLSurv")
```

## Overview

DRMLSurv implements a three-component pipeline for two-stage survival
data with dynamic treatment regimes (DTRs):

1.  **Censoring imputation** — imputes censored stage-specific survival
    times via constrained donor matching
    (`impute_censored_outcomes()`).
2.  **Counterfactual construction** — builds matched pseudo-outcomes
    under both treatment arms using nearest-neighbor or full matching
    (`matchpotential_DTR()`).
3.  **Policy learning** — estimates optimal dynamic treatment rules with
    random-forest classifiers and cross-fitting (`Drmatch()`,
    `rfdtr()`).

Key features:

- Nearest-neighbor and full matching with optional exact-matching
  constraints
- Balance diagnostics via
  [cobalt](https://ngreifer.github.io/cobalt/) Love plots
  (`plotbalance = TRUE`)
- Propensity and prognostic score estimation via SuperLearner or
  glm/glmnet (`ComputeScores()`)
- S3 `predict()` and `summary()` methods for fitted `Drmatch` objects

## Example

``` r
library(DRMLSurv)

# Construct matched counterfactual outcomes for a single treatment stage.
# 'dat' must contain columns for treatment (tx), outcome (Y), ID, and
# matching covariates (x1, x2).  Set plotbalance = TRUE to display a
# cobalt Love plot of covariate balance after matching.
matched <- matchpotential_DTR(
  dat         = dat,
  txgroup     = "tx",
  exact_vars  = NULL,
  compY       = "Y",
  vec         = list(c("x1", "x2"), c("x1", "x2")),  # ATT / ATC covariate sets
  Id          = "id",
  method      = "nearest",
  k           = 3,
  replace     = TRUE,
  distance    = "mahalanobis",
  plotbalance = TRUE   # requires the 'cobalt' package
)

# Fit the full two-stage DTR pipeline.
# 'mydata' must contain stage indicators, survival times, event indicators,
# treatment variables, and covariate columns listed in names.var1 / names.var2.
fit <- Drmatch(
  data        = mydata,
  id.var      = "id",
  eta2.var    = "eta2",
  Y1.var      = "Y1",
  Y2.var      = "Y2",
  delta.var   = "delta",
  OY.var      = "OY",
  A1.var      = "A1",
  A2.var      = "A2",
  names.var1  = c("age", "sex", "ecog"),
  names.var2  = c("age", "sex", "ecog", "Y1"),
  useds       = TRUE,
  superLearn  = FALSE,   # set TRUE to use SuperLearner ensemble
  plotbalance = FALSE,
  cap_months  = 36
)

# Predict optimal treatments for new patients
pred <- predict(fit, newdata = newdata, stage = "both")
head(pred)

# Summarise fitted DTR performance
summary(fit)
```

## Key functions

| Function | Purpose |
|---|---|
| `Drmatch()` | Fit the full two-stage DTR pipeline |
| `predict.Drmatch()` | Predict optimal treatments from a fitted model |
| `summary.Drmatch()` | Performance metrics for a fitted DTR |
| `impute_censored_outcomes()` | Two-stage censoring imputation |
| `matchpotential_DTR()` | Matched counterfactual construction |
| `rfdtr()` | Random-forest DTR policy learner with CV tuning |
| `ComputeScores()` | Propensity and prognostic score estimation |
| `propensityplot()` | Propensity score overlap plot |
| `policy_summary_metrics()` | Classification and value metrics |
