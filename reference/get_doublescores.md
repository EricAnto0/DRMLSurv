# Compute stage-specific treatment, prognostic, and optional censoring scores

Computes and attaches stage-specific score summaries for a two-stage
treatment setting by calling
[`ComputeScores`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
separately at stage 2 and stage 1.

The function is a wrapper that:

1.  restricts to subjects with `eta2 == 1` and computes stage-2 scores
    using `Y2.var`, `A2.var`, `names.var2`, and `Xtrt2`;

2.  computes stage-1 scores on the full cohort using either `OY.var` or
    `Y1.var` depending on `adjustdelta1`, together with `A1.var`,
    `names.var1`, and `Xtrt1`;

3.  renames the outputs from
    [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
    into stage-specific columns and merges them back into the original
    dataset by subject ID.

## Usage

``` r
get_doublescores(
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
  Xtrt1 = NULL,
  Xtrt2 = NULL,
  useds = FALSE,
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
  stratifyCV = TRUE,
  maxit = 1000,
  penalty1 = NULL,
  penalty2 = NULL,
  ngrid = 50,
  censmod = TRUE,
  pscens = TRUE,
  pgcens = TRUE,
  doublepg = TRUE,
  param.tune = NULL,
  adjustdelta1 = FALSE,
  plotps = FALSE,
  model.pg = "cox",
  standardize = FALSE,
  superLearn = TRUE,
  pslink = "logit",
  pglink = "lognormal",
  sl_parallel = c("multicore", "seq")
)
```

## Arguments

- data:

  A `data.frame` containing subject identifiers, stage indicators,
  outcomes, treatment variables, and covariates required for stage-1 and
  stage-2 score estimation.

- id.var:

  Character scalar. Name of the subject identifier column.

- eta2.var:

  Character scalar. Name of the stage-2 entry indicator column, where
  `1` denotes entry into stage 2 and `0` denotes no entry.

- Y1.var:

  Character scalar. Name of the stage-1 outcome or time variable used
  when `adjustdelta1 = TRUE`.

- Y2.var:

  Character scalar. Name of the stage-2 outcome or time variable.

- delta.var:

  Character scalar. Name of the event indicator variable used in the
  survival or censoring models.

- OY.var:

  Character scalar. Name of the overall outcome or time variable used
  for stage-1 scoring when `adjustdelta1 = FALSE`.

- A1.var:

  Character scalar. Name of the stage-1 treatment indicator variable.

- A2.var:

  Character scalar. Name of the stage-2 treatment indicator variable.

- names.var1:

  Character vector. Covariate names used in the stage-1 prognostic score
  model.

- names.var2:

  Character vector. Covariate names used in the stage-2 prognostic score
  model.

- Xtrt1:

  Character vector or `NULL`. Covariate names used in the stage-1
  treatment propensity model. If `NULL`,
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  uses `names.var1`.

- Xtrt2:

  Character vector or `NULL`. Covariate names used in the stage-2
  treatment propensity model. If `NULL`,
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  uses `names.var2`.

- useds:

  Logical. If `TRUE`, compute and merge the stage-specific scores. If
  `FALSE`, return `data` unchanged.

- cores:

  Integer. Number of cores passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for model fitting.

- tau:

  Optional numeric truncation horizon passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for restricted mean prediction or time-grid construction.

- sl.seed:

  Integer. Random seed passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).

- A.SL.library1:

  Character vector. SuperLearner library for the stage-1 treatment
  propensity model.

- A.SL.library2:

  Character vector. SuperLearner library for the stage-2 treatment
  propensity model.

- Y.SL.library:

  Character vector. Learners used for prognostic survival modeling
  inside
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).

- A.method:

  Optional character scalar. Performance metric passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for treatment propensity estimation.

- Y.method:

  Optional character scalar. Performance metric passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for prognostic survival estimation.

- param.weights.fix:

  Optional numeric vector. Fixed ensemble weights passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  when supported by the underlying learner.

- param.weights.init:

  Optional numeric vector. Initial ensemble weights passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  when supported by the underlying learner.

- optim.method:

  Character scalar or `NULL`. Optimization method forwarded to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).

- stratifyCV:

  Logical. Passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).
  If `TRUE`, cross-validation folds are stratified when supported by the
  underlying fitting procedure.

- maxit:

  Integer. Maximum number of optimization iterations passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).

- penalty1:

  Optional tuning parameter or penalty value passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for stage-1 prognostic estimation.

- penalty2:

  Optional tuning parameter or penalty value passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for stage-2 prognostic estimation.

- ngrid:

  Integer. Number of grid points used by
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  when approximating restricted means or evaluating predicted survival
  curves.

- censmod:

  Logical. If `TRUE`, request censoring-related scores from
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  in addition to treatment propensity and treatment prognostic scores.

- pscens:

  Logical. If `TRUE` and `censmod = TRUE`, estimate censoring propensity
  scores within
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).

- pgcens:

  Logical. If `TRUE` and `censmod = TRUE`, estimate censoring prognostic
  scores within
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).

- doublepg:

  Logical. Passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).
  If `TRUE`, estimate treatment-specific prognostic scores separately by
  treatment arm. If `FALSE`, the returned prognostic components may be
  partially unestimated and therefore remain `NA`.

- param.tune:

  Optional list or tuning object passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for learner-specific tuning.

- adjustdelta1:

  Logical. If `TRUE`, construct an adjusted stage-1 event indicator
  `deltaadj` and use `Y1.var` instead of `OY.var` in the stage-1 scoring
  call.

- plotps:

  Logical. If `TRUE`, plot the raw treatment propensity score
  distribution at each stage using
  [`propensityplot`](https://ericanto0.github.io/RMSurv/reference/propensityplot.md),
  when available.

- model.pg:

  Character scalar. Prognostic model type passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).
  Currently intended values are `"cox"` and `"aft"`.

- standardize:

  Logical. Passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).
  If `TRUE`, standardize covariates for penalized regression fits when
  applicable.

- superLearn:

  Logical. Passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).
  If `TRUE`, use SuperLearner-based fitting; otherwise use the
  parametric or penalized alternatives implemented there.

- pslink:

  Character scalar. Link function for binomial propensity models passed
  to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md),
  typically `"logit"` or `"probit"`.

- pglink:

  Character scalar. Distribution used when `model.pg = "aft"` inside
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md),
  for example `"exponential"`, `"weibull"`, `"lognormal"`, or
  `"loglogistic"`.

- sl_parallel:

  Character scalar. Parallel mode passed to
  [`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
  for SuperLearner-based fitting. Must be one of `"multicore"` or
  `"seq"`.

## Value

A `data.frame` equal to `data` augmented with stage-specific score
columns. If `useds = FALSE`, the original `data` is returned unchanged.

The following columns are attached for stage 1:

- `probps1`: raw treatment propensity score,

- `prog01`, `prog11`: raw treatment-specific prognostic scores,

- `probcens1`: raw censoring propensity score,

- `progcens1`: raw censoring prognostic score,

- `ps1`, `pg01`, `pg11`, `pscens1`, `pgcens1`: scaled versions of the
  above scores.

The analogous columns `probps2`, `prog02`, `prog12`, `probcens2`,
`progcens2`, `ps2`, `pg02`, `pg12`, `pscens2`, and `pgcens2` are
attached for stage 2. Subjects with `eta2 == 0` receive `NA` for all
stage-2 score columns.

If `doublepg = TRUE`, the convenience columns `pg1ct`, `pg1tc`, `pg2ct`,
and `pg2tc` are also added.

## Details

**Scores returned by
[`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).**

For each call,
[`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
returns a fixed set of score columns:

- `ps`: treatment propensity score,

- `pg0`, `pg1`: treatment-specific prognostic scores,

- `pscens`: censoring propensity score,

- `pgcens`: censoring prognostic score,

- `ps_sc`, `pg0_sc`, `pg1_sc`, `pscens_sc`, `pgcens_sc`: scaled versions
  of the corresponding raw scores.

This wrapper renames those outputs to stage-specific names and attaches
them to `data`. For stage 1, the suffix `1` is used; for stage 2, the
suffix `2` is used.

**Stage 2.**

Stage-2 scores are computed only among subjects satisfying
`data[[eta2.var]] == 1`. These scores are then merged back into the full
dataset. Subjects who do not enter stage 2 receive `NA` for all stage-2
score columns.

**Stage 1.**

Stage-1 scores are computed on the full dataset. By default, the stage-1
scoring call uses `Y = OY.var` and `event = delta.var`. If
`adjustdelta1 = TRUE`, the function instead uses `Y = Y1.var` and a
modified event indicator `deltaadj`, where `deltaadj` is initialized as
`delta.var` and then set to 0 for subjects with `eta2 == 1` and
`delta == 1`.

**Treatment and censoring covariates.**

The prognostic model covariates are supplied through `names.var1` and
`names.var2`. The treatment propensity model covariates are supplied
separately through `Xtrt1` and `Xtrt2`. If `Xtrt1` or `Xtrt2` is `NULL`,
then
[`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
uses the corresponding prognostic covariates.

**Censoring-related scores.**

If `censmod = TRUE`, the wrapper also requests censoring-related scores
from
[`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md).
The arguments `pscens` and `pgcens` determine whether censoring
propensity and censoring prognostic scores are actively estimated. Since
[`ComputeScores()`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md)
returns a fixed output structure, the corresponding columns are still
present in the returned data even when those components are not
estimated; in such cases they are typically `NA`.

**Treatment-specific prognostic scores.**

If `doublepg = TRUE`, the wrapper additionally creates convenience
variables comparing the observed-treatment and opposite-treatment
prognostic scores:

- `pg1ct`, `pg1tc` for stage 1,

- `pg2ct`, `pg2tc` for stage 2.

These are constructed from the standardized prognostic scores: `pg01`,
`pg11`, `pg02`, and `pg12`.

**No-op behavior.**

If `useds = FALSE`, the function returns `data` unchanged.

## See also

[`ComputeScores`](https://ericanto0.github.io/RMSurv/reference/ComputeScores.md),
[`propensityplot`](https://ericanto0.github.io/RMSurv/reference/propensityplot.md)
