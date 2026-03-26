
# DRMLSurv

<!-- badges: start -->

[![R-CMD-check](https://github.com/EricAnto0/DRMLSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EricAnto0/DRMLSurv/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![pkgdown](https://img.shields.io/badge/pkgdown-site-blue)](https://EricAnto0.github.io/DRMLSurv/)
[![codecov](https://codecov.io/gh/EricAnto0/DRMLSurv/branch/main/graph/badge.svg)](https://codecov.io/gh/EricAnto0/DRMLSurv)
<!-- badges: end -->

<img src="man/figures/DRMLSurvLogo.png" align="right" height="180" width="220" />

DRMLSurv is an R package for two-stage survival analysis with censoring,
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
remotes::install_github("EricAnto0/DRMLSurv")
```

## Load package

``` r
library(DRMLSurv)
```

## Example workflow

``` r
data("DATASET", package = "DRMLSurv")
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
                       lambda=seq(0.001, 0.1, length.out = 10))
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
Scores: censoring scores (SL): 3.745 sec elapsed
Stage 2 matching: 1.96 sec elapsed
Scores: censoring scores (SL): 16.049 sec elapsed
Stage 1 matching: 13.972 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Scores: treatment PS + pg0/pg1 (SL): 32.946 sec elapsed
Stage 2 DoubleScore: 32.947 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Scores: treatment PS + pg0/pg1 (SL): 110.1 sec elapsed
Stage 1 DoubleScore: 110.101 sec elapsed
Scores: censoring scores (SL): 3.672 sec elapsed
Stage 2 matching: 1.765 sec elapsed
Scores: censoring scores (SL): 15.574 sec elapsed
Stage 1 matching: 14.149 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Scores: treatment PS + pg0/pg1 (SL): 33.143 sec elapsed
Stage 2 DoubleScore: 33.144 sec elapsed
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
  |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |=========                                                             |  13%  |                                                                              |==============                                                        |  20%  |                                                                              |===================                                                   |  27%  |                                                                              |=======================                                               |  33%  |                                                                              |============================                                          |  40%  |                                                                              |=================================                                     |  47%  |                                                                              |=====================================                                 |  53%  |                                                                              |==========================================                            |  60%  |                                                                              |===============================================                       |  67%  |                                                                              |===================================================                   |  73%  |                                                                              |========================================================              |  80%  |                                                                              |=============================================================         |  87%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
Scores: treatment PS + pg0/pg1 (SL): 111.273 sec elapsed
Stage 1 DoubleScore: 111.275 sec elapsed
Best set of tuning parameters and metrics overall

 ntree = 1000
mtry = 2
nodesize = 2
CCR = 0.557317248755605
OOB = 0.442692123836303
Score = 7.84107735199969 
Best set of tuning parameters and metrics overall

 ntree = 1000
mtry = 2
nodesize = 10
CCR = 0.546726935592095
OOB = 0.459250100944479
Score = 12.4968720062604 
Score = 12.4992031699597 
#pred = predict(trainmod, newdata = test_data)
res <- summary(trainmod, newdata = test_data)
res
               Stage 1     Stage 2     Overall
Accuracy    0.53151101  0.45853659  0.20000000
MCC         0.07769275 -0.01458439 -0.01281208
Sensitivity 0.60979730  0.31147541  0.23973930
Specificity 0.46758621  0.67469880  0.74672252
PPV         0.48326640  0.58461538  0.23776803
NPV         0.59473684  0.40000000  0.74674813
F1          0.53920836  0.40641711  0.20188035
```
