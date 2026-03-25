# Dynamic Treatment Regime Learning with DRMLSurv

``` r
library(DRMLSurv)
```

## Introduction

`DRMLSurv` is an R package for dynamic treatment regime learning in
two-stage survival settings with censoring. The package is designed for
analyses in which treatment decisions are made sequentially, some event
times are censored, and the inferential target is a treatment rule
rather than a single static treatment contrast.

In many clinical applications, a subject receives an initial treatment
at stage 1, may or may not proceed to stage 2, and has an overall
survival outcome that may be incompletely observed. In such settings,
valid regime learning requires more than a direct comparison of observed
outcomes. One must account for censoring, construct appropriate
stage-specific and overall counterfactual outcomes, and then estimate
treatment rules using the completed or matched pseudo-outcomes.

The `DRMLSurv` workflow combines four main components:

1.  nuisance-score estimation, including treatment, prognostic, and
    censoring-related scores;
2.  matching-based imputation of censored stage-specific outcomes;
3.  construction of counterfactual regime-specific outcomes;
4.  policy learning using random forests or related machine-learning
    procedures.

This vignette gives a single end-to-end overview of the package workflow
and explains how the main functions fit together. The code examples are
intentionally schematic and are written to clarify the workflow. In
practice, the exact variable names, covariate sets, matching
restrictions, learner libraries, and tuning grids should reflect the
scientific problem and the structure of the available data.

## Conceptual setup

Consider a two-stage treatment setting. Let

- $A_{1} \in \{ - 1,1\}$ denote the stage-1 treatment,
- $A_{2} \in \{ - 1,1\}$ denote the stage-2 treatment among those who
  enter stage 2,
- $\eta_{2} \in \{ 0,1\}$ denote stage-2 entry,
- $Y_{1}$ denote stage-1 survival time,
- $Y_{2}$ denote stage-2 survival time,
- $Y$ denote overall survival time,
- $\Delta \in \{ 0,1\}$ denote the event indicator.

The target is a dynamic treatment regime
$$d = \left( d_{1},d_{2} \right),$$ where $d_{1}(\mathbf{x})$ assigns
stage-1 treatment based on baseline covariates and $d_{2}(\mathbf{x})$
assigns stage-2 treatment based on the relevant second-stage
information.

The general statistical objective is to identify a treatment rule
$d^{\star}$ that optimizes a policy value such as
$$V(d) = \mathbb{E}\{ Y^{d}\},$$ where $Y^{d}$ is the counterfactual
overall survival outcome under regime $d$.

Because $Y^{d}$ is not directly observed for all treatment paths, the
package builds completed regime-specific outcomes through matching and
imputation.

## Data requirements

A typical dataset for `DRMLSurv` contains:

- an individual identifier;
- a stage-2 entry indicator;
- stage-1 and stage-2 treatment variables;
- stage-specific and overall survival times;
- a censoring/event indicator;
- stage-specific covariates.

A common variable layout is:

- `patientid`: subject identifier;
- `eta2`: stage-2 entry indicator;
- `deathInd.raw`: event indicator;
- `OS_time`: observed overall survival time;
- `OS_time.1L`: stage-1 survival time;
- `OS_time.2L`: stage-2 survival time;
- `txgroup1L.sd`: stage-1 treatment;
- `txgroup2L0.sd`: stage-2 treatment.

If you have a package dataset, load it first.

``` r
data("DATASET", package = "DRMLSurv")
dat <- DATASET
```

You may then verify that the required variables are present.

``` r
required_vars <- c(
  "patientid", "eta2", "deathInd.raw",
  "OS_time", "OS_time.1L", "OS_time.2L",
  "txgroup1L.sd", "txgroup2L0.sd"
)

setdiff(required_vars, names(dat))
```

## Nuisance-score estimation

A central feature of the package is the estimation of nuisance
quantities used later for imputation, matching, and policy learning.
These include treatment scores, censoring scores, and prognostic scores.

At the lowest level, the function
[`ComputeScores()`](https://ericanto0.github.io/DRMLSurv/reference/ComputeScores.md)
estimates such quantities for a specified outcome, treatment, event
indicator, and covariate matrix.

For example, if one wants treatment and prognostic scores for a stage-1
survival outcome, one may use:

``` r
score_out <- ComputeScores(
  data       = dat,
  id         = "patientid",
  Y          = "OS_time.1L",
  event      = "deathInd.raw",
  X          = c("ageAt1L.sd", "gender.sd", "Albumin1st.sd", "Lymphocyte1st.sd"),
  A          = "txgroup1L.sd",
  doublepg   = TRUE,
  censmod    = FALSE,
  outer_CV   = 5,
  cores      = 1,
  superLearn = TRUE,
  model.pg   = "cox"
)

head(score_out)
```

In many applications, it is more convenient to compute stage-1 and
stage-2 scores together. This is what
[`get_doublescores()`](https://ericanto0.github.io/DRMLSurv/reference/get_doublescores.md)
is designed to do.

``` r
dat_scored <- get_doublescores(
  data          = dat,
  id.var        = "patientid",
  eta2.var      = "eta2",
  Y1.var        = "OS_time.1L",
  Y2.var        = "OS_time.2L",
  delta.var     = "deathInd.raw",
  OY.var        = "OS_time",
  A1.var        = "txgroup1L.sd",
  A2.var        = "txgroup2L0.sd",
  names.var1    = c("ageAt1L.sd", "gender.sd", "Albumin1st.sd", "Lymphocyte1st.sd"),
  names.var2    = c("ageAt1L.sd", "gender.sd", "Albumin2nd.sd", "Lymphocyte2nd.sd", "OS_time.1L"),
  Xtrt1         = NULL,
  Xtrt2         = NULL,
  useds         = TRUE,
  cores         = 1,
  tau           = NULL,
  sl.seed       = 123,
  A.SL.library1 = c("SL.mean", "SL.glm"),
  A.SL.library2 = c("SL.mean", "SL.glm"),
  Y.SL.library  = c("LIB_COXall"),
  A.method      = "method.AUC",
  Y.method      = "ibll",
  ngrid         = 200,
  censmod       = FALSE,
  doublepg      = TRUE,
  superLearn    = TRUE
)

names(dat_scored)
```

The output contains the subject identifier together with estimated score
variables, which are then used downstream as matching covariates.

## Imputation of censored outcomes

In two-stage survival analyses, some outcomes are only partially
observed because of censoring. The function
[`impute_censored_outcomes()`](https://ericanto0.github.io/DRMLSurv/reference/impute_censored_outcomes.md)
addresses this by using stage-specific donor matching to complete
censored outcomes.

Conceptually, the function performs two tasks. First, it imputes stage-2
outcomes among those who entered stage 2 but are censored. Second, it
uses completed stage-2 information to help construct completed overall
outcomes at stage 1.

A typical call is:

``` r
dat_imp <- impute_censored_outcomes(
  data          = dat,
  id.var        = "patientid",
  eta2.var      = "eta2",
  Y1.var        = "OS_time.1L",
  Y2.var        = "OS_time.2L",
  delta.var     = "deathInd.raw",
  OY.var        = "OS_time",
  A1.var        = "txgroup1L.sd",
  A2.var        = "txgroup2L0.sd",
  names.var1    = c("ageAt1L.sd", "gender.sd", "Albumin1st.sd", "Lymphocyte1st.sd"),
  names.var2    = c("ageAt1L.sd", "gender.sd", "Albumin2nd.sd", "Lymphocyte2nd.sd", "OS_time.1L"),
  exact1.vars   = NULL,
  exact2.vars   = NULL,
  usecov        = FALSE,
  useds         = FALSE,
  cores         = 1,
  tau           = NULL,
  sl.seed       = 123,
  A.SL.library1 = c("SL.mean", "SL.glm"),
  A.SL.library2 = c("SL.mean", "SL.glm"),
  Y.SL.library  = c("LIB_COXall"),
  A.method      = "method.AUC",
  Y.method      = "ibll",
  ngrid         = 200,
  pscens        = TRUE,
  pgcens        = TRUE,
  superLearn    = TRUE,
  distance      = "mahalanobis",
  method        = "nearest",
  K             = 1,
  replacement   = TRUE
)
```

The resulting dataset contains completed quantities needed for later
regime construction.

## Stage-specific censoring matching

If a more modular workflow is desired, the stage-1 and stage-2 censoring
completion steps can be run directly.

For stage 2, the function `match_censored_stage2()` uses donor matching
among subjects who entered stage 2.

``` r
res2 <- match_censored_stage2(
  dat        = dat[dat$eta2 == 1, , drop = FALSE],
  id.var     = "patientid",
  delta.var  = "deathInd.raw",
  OY.var     = "OS_time",
  Y2.var     = "OS_time.2L",
  formula2   = (deathInd.raw == 0) ~ pg02 + pg12 + ps2,
  exact.vars = NULL,
  method     = "nearest",
  distance   = "mahalanobis",
  k          = 1,
  replace    = TRUE,
  aggregate  = "mean"
)

res2$n_imputed
head(res2$data_merged)
```

For stage 1, the analogous function is `match_censored_stage1()`.

``` r
res1 <- match_censored_stage1(
  dat        = dat,
  id.var     = "patientid",
  exact_vars = NULL,
  formula    = (DeathInd1 == 0) ~ pg01 + pg11 + ps1,
  death1     = "DeathInd1",
  OY         = "OS_time",
  y1         = "OS_time.1L",
  distance   = "mahalanobis",
  method     = "nearest",
  k          = 1,
  replace    = TRUE,
  y_cols     = c("CompOS_time", "CompOS_2LChemo1Lobs", "CompOS_2LCheMon1Lobs"),
  aggregate  = "mean"
)

res1$n_imputed
head(res1$data_merged)
```

These functions are helpful when one wants direct control over the
censoring imputation stage rather than relying on the higher-level
wrapper.

## Stage-2 counterfactual construction

After stage-2 outcomes are completed, the next task is to construct
opposite-arm matched outcomes under alternative second-stage treatments.
This is done by
[`matchpotential_DTR()`](https://ericanto0.github.io/DRMLSurv/reference/matchpotential_DTR.md).

A typical stage-2 counterfactual matching call is:

``` r
dat_2L <- matchpotential_DTR(dat = dat_2L, txgroup = "txgroup1L.sd", exact_vars = NULL, 
                   compY = 'compY2', vec = c("pg02", "pg12", "ps2"), Id = 'patientid',
                   method = "nearest", k =1, 
                   replace = TRUE, caliper = NULL,
                   distance = "mahalanobis",
                   compW = "ipcw.R")
```

The returned object includes matched opposite-arm outcomes such as
paired stage-2 survival or completed stage-2 outcomes. These are then
used to evaluate second-stage treatment rules.

## Stage-1 counterfactual construction

Stage-1 counterfactuals are built similarly. The function
[`matchpotential_DTR()`](https://ericanto0.github.io/DRMLSurv/reference/matchpotential_DTR.md)
constructs paired stage-1 regime-specific outcome quantities.

A typical stage-1 matching call is:

``` r
dat_1L_final <- matchpotential_DTR(dat = dat_1L, txgroup = "txgroup1L.sd", exact_vars = NULL, 
                   compY = 'compOY', vec = c("pg01", "pg11", "ps1"), Id = 'patientid',
                   method = "nearest", k =1, 
                   replace = TRUE, caliper = NULL,
                   distance = "mahalanobis",
                   compW = "ipcw.R")
```

This matched outcomes form the pseudo-responses needed for
policy-learning algorithms.

## Policy learning with random forests

The package includes
[`rfdtr()`](https://ericanto0.github.io/DRMLSurv/reference/rfdtr.md) for
random-forest-based policy learning. This function uses matched labels,
matched outcomes, and optional hyperparameter tuning to learn treatment
rules.

For a stage-2 analysis, one may form a learning dataset such as:

``` r
obs_stage2 <- data.frame(
  A = dat_2L$newTxt,
  dat_2L[, c("ageAt1L.sd", "gender.sd", "Albumin2nd.sd", "Lymphocyte2nd.sd"), drop = FALSE]
)

gridpar <- expand.grid(
  mtry     = c(1, 2, 3),
  ntree    = c(500),
  nodesize = c(2, 5)
)

rf_fit <- rfdtr(
  modeltype = "ranger",
  usecv     = TRUE,
  sl.seed   = 123,
  obs       = obs_stage2,
  W         = dat_2L$match.weight,
  gridpar   = gridpar,
  metric    = "policyval",
  A.obs     = dat_2L$txgroup2L0.sd,
  Q.obs     = dat_2L$compOS2L,
  Q.match   = dat_2L$pairedCompOS2L
)

rf_fit$best
```

The returned object contains the fitted model, estimated optimal
treatments for the observed sample, a tuning summary, and the selected
best tuning parameters.

## Integrated fitting with `Drmatch()`

For end-to-end analysis, the package provides
[`Drmatch()`](https://ericanto0.github.io/DRMLSurv/reference/Drmatch.md),
which wraps imputation, matching, score construction, and policy
learning in a single pipeline.

A typical call is:

``` r
fit <-  trainmod = Drmatch(
    data                       = data,
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
```

## Summary and prediction

Once a `Drmatch` object has been fitted, the fitted rule can be
summarized using the S3 summary method.

``` r
summary(fit, newdata = dat)
```

Predictions for new data may be obtained with the corresponding
prediction method.

``` r
pred <- predict(fit, newdata = dat)
head(pred)
```

Policy performance on held-out or validation data can then be summarized
with
[`policy_summary_metrics()`](https://ericanto0.github.io/DRMLSurv/reference/policy_summary_metrics.md).

``` r
metrics <- policy_summary_metrics(
  tmpData = dat,
  estA1   = pred$estA1,
  estA2   = pred$estA2
)

metrics
```

## Practical recommendations

The full `DRMLSurv` workflow can be computationally demanding. This is
especially true when large SuperLearner libraries, large tuning grids,
repeated cross-validation, or large datasets are used. During
development and testing, it is often sensible to begin with:

- a smaller subset of the data;
- a reduced learner library;
- a small tuning grid;
- serial computation with `cores = 1`.

Special attention should also be given to exact-matching variables,
stage-specific covariate sets, censoring-model formulas, and the
interpretation of completed regime-specific outcomes.

## Development workflow

During package development, a practical sequence is:

``` r
devtools::document()
devtools::test()
devtools::check()
```

When the package README is updated, it can be rebuilt with:

``` r
devtools::build_readme()
```

## Conclusion

`DRMLSurv` provides a flexible framework for dynamic treatment regime
learning in two-stage survival settings with censoring. The package is
built around three key ideas: donor-based completion of partially
observed outcomes, matching-based construction of regime-specific
counterfactuals, and machine-learning-based estimation of treatment
rules. The component functions allow each analytical stage to be
inspected separately, while
[`Drmatch()`](https://ericanto0.github.io/DRMLSurv/reference/Drmatch.md)
offers a practical integrated entry point for end-to-end analysis.

In applied work, the most important modeling decisions concern the
choice of covariates, the matching strategy, the nuisance-score
estimation procedure, and the policy-learning objective. The package is
designed to support those decisions while keeping the overall workflow
transparent and reproducible.
