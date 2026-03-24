# Two-stage DTR learning/evaluation pipeline with cross-fitting

`simrun()` runs one replicate (indexed by `jj`) of a two-stage dynamic
treatment regime (DTR) learning and evaluation pipeline using outer
cross-fitting (default 5 folds). Within each fold, the function imputes
censoring-related composite outcomes on the training partition,
constructs double scores (propensity/prognostic scores) when requested,
performs matching-based pseudo-outcome construction for stage 2 and
stage 1, fits random-forest policies via [`rfdtr`](rfdtr.md), applies
the learned rules to the held-out fold, and returns fold-level value
summaries and benchmark comparisons.

## Usage

``` r
Drmatch(
  data = data,
  id.var = "PatientID",
  eta2.var = "eta2",
  Y1.var = "OS_time.1L",
  Y2.var = "OS_time.2L",
  delta.var = "deathInd.raw",
  OY.var = "OS_time",
  A1.var = "txgroup1L.sd",
  A2.var = "txgroup2L0.sd",
  exact1.vars = c("txgroup1L.sd", "gender.sd", "ECOG1st0.sd", "ECOG1st1.sd"),
  exact2.vars = c("txgroup1L.sd", "gender.sd", "ECOG2nd0.sd", "ECOG2nd1.sd",
    "txgroup2L0.sd"),
  names.var1 = c("ageAt1L", "gender.sd", "ECOG1st0.sd", "ECOG1st1.sd", "Albumin1st",
    "Lymphocyte1st"),
  names.var2 = c("ageAt1L", "gender.sd", "ECOG2nd0.sd", "Lymphocyte2nd", "OS_time.1L",
    "ECOG2nd1.sd", "Albumin2nd", "Lymphocyte2nd"),
  usecov = FALSE,
  cores = 5,
  tau = NULL,
  sl.seed = 1234,
  A.SL.library = c("SL.glm", "SL.glmnet", "SL.ranger"),
  Y.SL.library = c("LIB_COXlasso", "LIB_COXen", "LIB_AFTggamma", "LIB_RSF"),
  A.method = "method.NNloglik",
  Y.method = "ibll",
  plotps = FALSE,
  ngrid = 1000,
  param.tune = NULL,
  param.weights.fix = NULL,
  param.weights.init = NULL,
  optim.method = "Nelder-Mead",
  maxit = 10000,
  usepenalty = FALSE,
  runseed = 2025,
  useds = TRUE,
  Taus = NULL,
  adjustdelta1 = FALSE,
  modeltype = "ranger",
  usecv = TRUE,
  doublepg = TRUE,
  model.pg = "cox",
  standardize = FALSE,
  superLearn = TRUE,
  pslink = "logit",
  pglink = "lognormal",
  distance = "mahalanobis",
  method = "nearest",
  K = 3,
  replacement = TRUE,
  cap_months = 36
)
```

## Arguments

- data:

  A `data.frame` containing all required columns used by `Get_data()`,
  matching, scoring, and policy learning.

- id.var:

  Name of unique subject identifier column.

- eta2.var:

  Name of stage-2 entry indicator (1 = enters stage 2).

- Y1.var:

  Name of stage-1 observed time/outcome component.

- Y2.var:

  Name of stage-2 observed time/outcome component.

- delta.var:

  Name of event indicator (1 = event, 0 = censored).

- OY.var:

  Name of overall observed time (used as a boundary in donor filtering).

- A1.var:

  Name of stage-1 observed treatment/action.

- A2.var:

  Name of stage-2 observed treatment/action.

- exact1.vars:

  Character vector of variables used for exact matching at stage 1. Use
  `NULL` or `character(0)` to disable exact matching.

- exact2.vars:

  Character vector of variables used for exact matching at stage 2. Use
  `NULL` or `character(0)` to disable exact matching.

- names.var1:

  Character vector of covariate names used at stage 1
  (scores/matching/model).

- names.var2:

  Character vector of covariate names used at stage 2
  (scores/matching/model).

- usecov:

  Logical; if `TRUE`, matching formulas are built from
  `names.var1`/`names.var2`. If `FALSE`, matching may rely on score
  variables (e.g., pg/ps) produced by the pipeline.

- cores:

  Integer number of workers used for internal parallel steps.

- tau:

  Optional horizon(s) passed to score construction; in the shown
  implementation, `rtau` is set to `cap_months` for both stages.

- sl.seed:

  Seed passed to SuperLearner / score construction components.

- A.SL.library:

  Candidate learners for treatment (propensity) score models.

- Y.SL.library:

  Candidate learners for outcome/censoring models used in
  [`ComputeScores()`](ComputeScores.md).

- A.method:

  Optional optimization method for propensity fitting in score
  construction.

- Y.method:

  Optional optimization method for outcome fitting in score
  construction.

- plotps:

  Logical; if `TRUE`, generate propensity score overlap plots where
  implemented.

- ngrid:

  Integer; grid size for internal optimization routines used by score
  construction.

- param.tune:

  Optional tuning structure passed into score construction.

- param.weights.fix, param.weights.init:

  Optional fixed/initial weights passed to score construction.

- optim.method:

  Optimization method used in score construction (e.g.,
  `"Nelder-Mead"`).

- maxit:

  Maximum iterations for score construction optimizers.

- usepenalty:

  Logical; if `TRUE`, define `SL.glmnet.tune1`/`SL.glmnet.tune2`
  wrappers with `penalty.factor` for stage-specific propensity
  estimation.

- runseed:

  Integer vector of seeds; replicate `jj` uses `runseed[jj]`.

- useds:

  Logical; if `TRUE`, compute/use double scores via
  [`get_doublescores`](get_doublescores.md).

- Taus:

  Optional matrix/data.frame of candidate horizons
  (legacy/compatibility).

- adjustdelta1:

  Logical; if `TRUE`, modifies stage-1 censoring indicator as in
  [`impute_censored_outcomes()`](impute_censored_outcomes.md) and
  [`get_doublescores()`](get_doublescores.md).

- modeltype:

  Character; model type for policy learning in [`rfdtr`](rfdtr.md)
  (e.g., `"ranger"`).

- usecv:

  Logical; if `TRUE`, tune policy learner with internal CV in
  [`rfdtr`](rfdtr.md).

- doublepg:

  Logical; if `TRUE`, use arm-specific prognostic scores (pg0/pg1) where
  applicable.

- model.pg:

  Character; prognostic model family (e.g., `"cox"` or `"aft"`) for
  score construction.

- standardize:

  Logical; if `TRUE`, standardize covariates for `glmnet` components.

- superLearn:

  Logical; if `TRUE`, use SuperLearner; if `FALSE`, use glmnet-only path
  where implemented by [`ComputeScores()`](ComputeScores.md).

- pslink:

  Character; link for propensity estimation in score construction (e.g.,
  `"logit"`).

- pglink:

  Character; distribution/link for prognostic score modeling
  (project-specific; e.g., `"lognormal"`).

- distance:

  Matching distance used by
  [`matchpotential_DTR`](matchpotential_DTR.md) (e.g., `"mahalanobis"`).

- method:

  Matching method used in [`matchpotential_DTR`](matchpotential_DTR.md)
  (e.g., `"nearest"` or `"full"`).

- K:

  Integer donor ratio for matching (number of donors per anchor).

- replacement:

  Logical; whether to allow donor reuse in nearest-neighbor matching.

- cap_months:

  Numeric truncation cap (months). Value summaries are computed after
  capping at this horizon.

## Value

A `data.frame` with one row per outer fold (default 5). Columns include:

- Fold-level value summaries (means and medians) under the learned DTR
  and under several alternative constructions produced inside the
  pipeline (e.g., `Ttot`, `T2`, `T2.rf`, `TtotOSml`, `Ttot1LOSmed`,
  etc.).

- Agreement metrics: `CCR1` (stage-1 agreement) and `CCR2` (stage-2
  agreement among `eta2==1`).

- Observed and uniform-regime benchmarks (means/medians) at stage 2 and
  overall.

- Diagnostics stored as list-columns, including RF best hyperparameters
  and predicted optimal actions.

The returned object also includes `iter` (replicate index) and `fold`.
If a fatal error occurs, `NULL` is returned.

## Details

The outer loop is a 5-fold split created by
[`caret::createFolds()`](https://rdrr.io/pkg/caret/man/createDataPartition.html)
on the stage-1 observed treatment `A1.var`. For each fold:

1.  Split into training and test partitions.

2.  Standardize selected time variables in training (`OY.sd`, `Y1.sd`,
    `Y2.sd`).

3.  Impute censoring-related composite outcomes on the training fold via
    [`impute_censored_outcomes`](impute_censored_outcomes.md)
    (internally calling stage-2 then stage-1 imputation).

4.  Optionally compute double scores on the imputed training fold via
    [`get_doublescores`](get_doublescores.md).

5.  Stage 2 (among `eta2==1`): build matched opposite-arm
    pseudo-outcomes via [`matchpotential_DTR`](matchpotential_DTR.md),
    tune/fit a random-forest policy via [`rfdtr`](rfdtr.md), and compute
    an estimated stage-2 optimal value by substituting matched outcomes
    for units assigned to a non-optimal observed treatment.

6.  Stage 1: construct a pseudo-outcome `Rtilde` (including the
    estimated stage-2 optimal value for stage-2 entrants), perform
    matching via [`matchpotential_DTR`](matchpotential_DTR.md), and
    tune/fit the stage-1 policy via [`rfdtr`](rfdtr.md).

7.  Apply learned policies to the held-out fold and compute fold-level
    summaries: value estimates (means and medians, capped at
    `cap_months`), agreement measures (`CCR1`, `CCR2`), and
    observed/uniform benchmark regimes.

## See also

`Get_data` `Get_data`,
[`impute_censored_outcomes`](impute_censored_outcomes.md),
[`get_doublescores`](get_doublescores.md),
[`matchpotential_DTR`](matchpotential_DTR.md), [`rfdtr`](rfdtr.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires project-specific helpers: Get_data(), impute_censored_outcomes(),
# get_doublescores(), matchpotential_DTR(), rfdtr(), and my_score.Surv().
#
# res <- simrun(
#   data = your_data,
#   runseed = 2025,
#   cores = 2
# )
# head(res)
} # }
```
