
# RMSurv <img src="man/figures/RMSurvLogo.png" align="right" height="180" width="220" /></a>

<!-- badges: start -->

![](https://cranlogs.r-pkg.org/badges/grand-total/RMSurv)
[![R-CMD-check](https://github.com/EricAnto0/RMSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EricAnto0/RMSurv/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![pkgdown](https://img.shields.io/badge/pkgdown-site-blue)](https://EricAnto0.github.io/RMSurv/)
[![codecov](https://codecov.io/gh/EricAnto0/RMSurv/branch/main/graph/badge.svg)](https://codecov.io/gh/EricAnto0/RMSurv)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

------------------------------------------------------------------------

RMSurv is an R package for two-stage survival analysis with censoring,
matching-based imputation, counterfactual outcome construction, and
machine-learning estimation of dynamic treatment rules.

Main features include:

- estimate propensity and prognostic (and censoring) scores with
  SuperLearner
- donor-based imputation of censored stage-1 and stage-2 survival times
- matched counterfactual outcome construction under alternative
  treatment paths
- policy learning with random forests and cross-validation
- obtaining optimized regimes and policy summary metrics

## Installation

``` r
# install.packages("remotes")
remotes::install_github("EricAnto0/RMSurv")
```

## Load package

``` r
library(RMSurv)
```

## Example workflow

``` r
data("DATASET", package = "RMSurv")
dat <- DATASET

set.seed(123)
folds <- caret::createFolds(dat$txgroup1L.sd, k = 3)

test_idx  <- folds[[1]]
train_idx <- setdiff(seq_len(nrow(dat)), test_idx)

train_data <- dat[train_idx, , drop = FALSE]
test_data  <- dat[test_idx, , drop = FALSE]

 trainmod = Drmatch(
    data                       = train_data,
    id.var                     = 'patientid',
    eta2.var                   = 'eta2',
    Y1.var                     = 'OS_time.1L',
    Y2.var                     = 'OS_time.2L',
    delta.var                  = 'deathInd.raw',
    OY.var                     = 'OS_time',
    A1.var                     = 'txgroup1L.sd',
    A2.var                     = 'txgroup2L0.sd',
    names.var1                 = c("ageAt1L.sd","gender.sd","Albumin1st.sd","Lymphocyte1st.sd",
                                   "ECOG1st0.sd","ECOG1st1.sd", 'firstLineStartTime.sd'),
    names.var2                 =  c("ageAt1L.sd","OS_time.1L","Albumin2nd.sd",
                                    "gender.sd","Lymphocyte2nd.sd","ECOG2nd0.sd",
                                    "ECOG2nd1.sd", "txgroup1L.sd", 'firstLineStartTime.sd'),
    cores                      = 4,
    sl.seed                    = 1234,
    A.SL.library              = list(
      "SL.ranger", 
      "SL.glm", 
      "SL.glmnet"), 
    Y.SL.library              = c("LIB_COXlasso", "LIB_COXall",  "LIB_COXen"
    ),
    A.method                   = "method.NNloglik",
    Y.method                   = "ibll",
    plotps                     = FALSE,
    ngrid                      = 5000,
    param.tune                 =  list(
      LIB_COXlasso = list(lambda=seq(0.001, 0.25, length.out = 10)),
      LIB_COXall = NULL,
      LIB_COXen = list(alpha=seq(.1, .9, length.out = 10),
                       lambda=seq(0.001, 0.1, length.out = 10)) #NULL,
      # LIB_RSF = list(
      #   mtry = unique(round(c(sqrt(length(names.var2)), seq(1, floor(length(names.var2)-1), length.out = 10)))),#round(nX/3), # Number of variables to consider at each split
      #   nodesize = 5,
      #   ntree = c(1000)
      # )
    ),
    maxit                      = 10000,
    runseed                    = 2025,
    useds                      = TRUE,
    modeltype                  = "ranger",
    usecv                      = TRUE,
    doublepg                   = TRUE,
    model.pg                   = "cox",
    superLearn                 = TRUE,
    distance                   = 'mahalanobis',
    method                     = 'nearest',
    K                          = 3,
    replacement                = TRUE,
    cap_months                 = 24
  )
Treatment PS: CV.SuperLearner: 3.806 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Treatment PG1: survivalSL: 24.935 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Treatment PG0: survivalSL: 11.697 sec elapsed
Censoring PS: CV.SuperLearner: 3.695 sec elapsed
Time to compute stage-2 scores: 45.774 sec elapsed
Treatment PS: CV.SuperLearner: 13.16 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Treatment PG1: survivalSL: 95.952 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Treatment PG0: survivalSL: 62.751 sec elapsed
Censoring PS: CV.SuperLearner: 12.223 sec elapsed
Time to compute stage-1 scores: 186.957 sec elapsed
obtain the double scores for training fold: 232.769 sec elapsed
Stage 2 matching: 2.975 sec elapsed
Stage 1 matching: 23.04 sec elapsed
Imputation of censored time: 26.037 sec elapsed
Best set of tuning parameters and metrics overall

 ntree = 1000
mtry = 2
nodesize = 10
CCR = 0.592706405035172
OOB = 0.441937444480808
Score = 7.75173579249977 
Best set of tuning parameters and metrics overall

 ntree = 1000
mtry = 3
nodesize = 2
CCR = 0.566101254680702
OOB = 0.441964526192555
Score = 12.0718535827058 
pred = predict(trainmod, newdata = test_data); head(pred)
  row_id A1.opt A2.opt
1      1      1     -1
2      2      1     NA
3      3      1     NA
4      4     -1     NA
5      5      1     NA
6      6     -1     NA
```
