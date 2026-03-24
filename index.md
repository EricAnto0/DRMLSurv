# DRMLSurv

![](reference/figures/logo.png)

DRMLSurv is an R package for two-stage survival analysis with censoring,
matching-based imputation, counterfactual outcome construction, and
machine-learning estimation of dynamic treatment rules.

Main features include:

- donor-based imputation of censored stage-1 and stage-2 survival times
- matched counterfactual outcome construction under alternative
  treatment paths
- policy learning with random forests and cross-validation
- policy summary metrics and balance diagnostics

## Installation

``` r
# install.packages("remotes")
remotes::install_github("EricAnto0/DRMLSurv")
```

# load package

``` r
library(DRMLSurv)
```
