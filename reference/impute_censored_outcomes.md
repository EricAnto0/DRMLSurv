# Impute censored stage-specific outcomes via matching with optional learned censoring scores

Performs two-stage imputation of censored outcomes in a sequential
treatment / DTR setting. The function (i) constructs a stage-2 composite
outcome among subjects who enter stage 2 (`eta2==1`) by matching
censored subjects (`delta2==0`) to eligible uncensored donors, and (ii)
constructs a stage-1 composite outcome for all subjects by matching
stage-1 censored subjects (`delta1==0`) to eligible donors, optionally
leveraging learned censoring propensity and prognostic scores in the
matching distance.

Matching is performed by helper routines
[`impute_censored_stage2()`](https://ericanto0.github.io/DRMLSurv/reference/impute_censored_stage2.md)
and
[`impute_censored_stage1()`](https://ericanto0.github.io/DRMLSurv/reference/impute_censored_stage1.md),
which implement nearest-neighbor or optimal matching using Mahalanobis
or other distances, with optional exact matching constraints. When
`usecov=FALSE`, matching covariates are replaced by learned censoring
scores computed by
[`ComputeScores()`](https://ericanto0.github.io/DRMLSurv/reference/ComputeScores.md):
`pscens` (a censoring propensity score) and/or `pgcens` (a censoring
prognostic score).

## Usage

``` r
impute_censored_outcomes(
  data,
  id.var,
  eta2.var,
  Y1.var,
  Y2.var,
  delta.var,
  OY.var,
  A1.var,
  A2.var,
  names.var1,
  names.var2,
  exact1.vars = character(0),
  exact2.vars = character(0),
  usecov = TRUE,
  useds = FALSE,
  adjustdelta1 = FALSE,
  cores = 1,
  tau,
  sl.seed = 123,
  A.SL.library1,
  A.SL.library2,
  Y.SL.library,
  A.method = NULL,
  Y.method = NULL,
  param.weights.fix = NULL,
  param.weights.init = NULL,
  optim.method = NULL,
  maxit = 1000,
  penalty1 = NULL,
  penalty2 = NULL,
  ngrid = 50,
  pscens = TRUE,
  pgcens = TRUE,
  param.tune = NULL,
  plotps = FALSE,
  model.pg = "cox",
  standardize = FALSE,
  superLearn = TRUE,
  pslink = "logit",
  distance = "mahalanobis",
  method = "nearest",
  K = 3,
  replacement = TRUE
)
```

## Arguments

- data:

  A data.frame containing all variables needed for stage-1 and stage-2
  processing.

- id.var:

  Character scalar. Subject identifier column name.

- eta2.var:

  Character scalar. Stage-2 entry indicator column name (1=entered stage
  2, 0=did not).

- Y1.var:

  Character scalar. Stage-1 time/outcome component used in the composite
  outcome.

- Y2.var:

  Character scalar. Stage-2 time/outcome component to be imputed for
  censored stage-2 subjects.

- delta.var:

  Character scalar. Stage-2 event indicator column name (1=observed,
  0=censored).

- OY.var:

  Character scalar. Overall outcome column name (used when `eta2==0` and
  observed).

- A1.var:

  Character scalar. Stage-1 treatment column name.

- A2.var:

  Character scalar. Stage-2 treatment column name.

- names.var1:

  Character vector. Covariate names available at stage 1 (used when
  `usecov=TRUE`).

- names.var2:

  Character vector. Covariate names available at stage 2 (used when
  `usecov=TRUE`).

- exact1.vars, exact2.vars:

  Character vectors. Variables used for exact matching at stage 1 /
  stage 2. Default is none.

- usecov:

  Logical. If TRUE, matching uses the covariates in `names.var1` /
  `names.var2`. If FALSE, matching uses learned censoring scores
  (`pscens` and/or `pgcens`).

- useds:

  Logical. Reserved for future use (currently not used in the provided
  implementation).

- adjustdelta1:

  Logical. If TRUE, modifies the stage-1 event indicator to treat
  certain stage-2 events as censored by creating `deltaadj`. (Current
  code sets `deltaadj=0` for `eta2==1 & delta==1`.)

- cores:

  Integer. Requested number of cores for downstream scoring routines (if
  supported).

- tau:

  Optional numeric. Truncation horizon passed to scoring routines for
  mean survival calculations.

- sl.seed:

  Integer. RNG seed passed to SuperLearner-based scoring.

- A.SL.library1, A.SL.library2:

  Character vectors. SuperLearner libraries for stage-1 and stage-2
  censoring propensity models (passed to `ComputeScores`).

- Y.SL.library:

  Character vector. Learners for survivalSL prognostic modeling (passed
  through).

- A.method, Y.method:

  Optional characters. Risk/metric identifiers for SuperLearner /
  survivalSL.

- param.weights.fix, param.weights.init, optim.method, maxit, penalty1,
  penalty2, param.tune:

  Tuning and optimization controls passed to `ComputeScores` /
  survivalSL scoring.

- ngrid:

  Integer. Number of grid points for survival-curve integration in
  prognostic scoring.

- pscens, pgcens:

  Logical. Whether to compute and use censoring propensity (`pscens`)
  and/or censoring prognostic (`pgcens`) scores when `usecov=FALSE`.

- plotps:

  Logical. If TRUE, produces diagnostic propensity plots for censoring
  scores.

- model.pg:

  Character. Prognostic modeling family when `superLearn=FALSE` ("cox"
  or "aft").

- standardize:

  Logical. Whether to standardize covariates for glmnet when
  `superLearn=FALSE`.

- superLearn:

  Logical. If TRUE, uses SuperLearner-based scoring; otherwise uses
  glm/glmnet-based scoring.

- pslink:

  Character. Link for binomial models ("logit" or "probit").

- distance:

  Character. Distance type used by matching routines (e.g.,
  "mahalanobis").

- method:

  Character. Matching method ("nearest" or "optimal").

- K:

  Integer. Donor ratio (number of matched donors per censored subject).

- replacement:

  Logical. Whether donors can be reused across matches.

## Value

A data.frame containing the original data augmented with
imputed/composite outcomes. At minimum, the output includes:

- `compY2`: imputed/observed stage-2 component

- `compOY`: composite overall outcome used for stage-1 imputation

- `delta1`, `delta2`: derived indicators of composite outcome
  observability

If `usecov=FALSE`, additional columns containing raw and standardized
censoring scores may also be present, depending on `pscens` and
`pgcens`.

## Details

**Stage-2 imputation.** The function subsets to `eta2==1` and imputes
`Y2` for censored subjects using donors with observed stage-2 outcomes,
subject to donor eligibility constraints encoded in
[`impute_censored_stage2()`](https://ericanto0.github.io/DRMLSurv/reference/impute_censored_stage2.md)
(e.g., donor observed time exceeding the recipient’s observed time). The
result is stored as `compY2`. For uncensored subjects, `compY2` is set
to the observed `Y2`.

**Stage-1 composite outcome.** A composite outcome `compOY` is
constructed as `compOY = Y1 + compY2` for stage-2 entrants; for
non-entrants (`eta2==0`) with observed final outcome, `compOY` is set to
the observed overall outcome `OY`. Stage-1 censoring indicators `delta1`
and `delta2` are defined from whether `compOY` is observed. If
`adjustdelta1=TRUE`, the stage-1 event indicator is modified to treat
certain stage-2 events as censored (see the `adjustdelta1` section
below).

**Censoring-score matching (optional).** If `usecov=FALSE`, censoring
scores are computed via `ComputeScores(censmod=TRUE, doublepg=FALSE)`
separately for stage 2 and stage 1. The raw probability scores are
optionally logit-transformed and standardized (z-scored) prior to
matching. These scores are then used as matching covariates in place of
the original covariates.

**Important assumption.** This function assumes `delta==1` denotes an
observed event/outcome (i.e., not censored) and `delta==0` denotes
censoring, consistent with `Surv(time, delta)`. If your data use the
opposite convention, you must recode before calling this function.

## See also

[`ComputeScores`](https://ericanto0.github.io/DRMLSurv/reference/ComputeScores.md),
`impute_censored_stage2`, `impute_censored_stage1`
