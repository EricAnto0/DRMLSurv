# Compute stage-specific treatment propensity and prognostic “double scores”

Computes and attaches stage-specific *double scores* for a two-stage
treatment setting. The function is a thin orchestrator around
[`ComputeScores`](ComputeScores.md) that: (i) restricts to stage-2
entrants (`eta2==1`) to compute stage-2 scores using `A2` and `Y2`, then
(ii) computes stage-1 scores on the full cohort using `A1` and either
`OY` (overall outcome) or `Y1` (stage-1 time) depending on
`adjustdelta1`.

The output is the original dataset augmented with both *raw* score
columns (propensities and prognostic scores) and *standardized* score
columns intended for distance-based matching or downstream modeling.

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
  Xtrt1,
  Xtrt2,
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
  maxit = 1000,
  penalty1 = NULL,
  penalty2 = NULL,
  ngrid = 50,
  censmod = TRUE,
  doublepg = TRUE,
  param.tune = NULL,
  adjustdelta1 = FALSE,
  plotps = FALSE,
  model.pg = "cox",
  standardize = FALSE,
  superLearn = TRUE,
  pslink = "logit",
  pglink = NULL
)
```

## Arguments

- data:

  A data.frame containing all required stage-1 and stage-2 variables.

- id.var:

  Character scalar. Subject identifier column name.

- eta2.var:

  Character scalar. Stage-2 entry indicator column name (1=entered stage
  2, 0=did not).

- Y1.var:

  Character scalar. Stage-1 time/outcome component (used only when
  `adjustdelta1=TRUE`).

- Y2.var:

  Character scalar. Stage-2 outcome/time column used for stage-2
  prognostic scoring.

- delta.var:

  Character scalar. Event indicator column name used for survival
  modeling.

- OY.var:

  Character scalar. Overall outcome/time column used for stage-1
  prognostic scoring when `adjustdelta1=FALSE`.

- A1.var:

  Character scalar. Stage-1 treatment indicator column name.

- A2.var:

  Character scalar. Stage-2 treatment indicator column name.

- names.var1:

  Character vector. Covariate names for the stage-1 prognostic model.

- names.var2:

  Character vector. Covariate names for the stage-2 prognostic model
  (stage-2 entrants only).

- Xtrt1:

  Character vector. Covariate names for the stage-1 treatment propensity
  model (if different from `names.var1`).

- Xtrt2:

  Character vector. Covariate names for the stage-2 treatment propensity
  model (if different from `names.var2`).

- useds:

  Logical. If TRUE, compute and merge scores. If FALSE, return `data`
  unchanged.

- cores:

  Integer. Number of cores passed to `ComputeScores` (if supported by
  the backend).

- tau:

  Optional numeric. Truncation horizon used in prognostic mean
  calculations inside `ComputeScores`.

- sl.seed:

  Integer. RNG seed passed to `ComputeScores`.

- A.SL.library1, A.SL.library2:

  Character vectors. SuperLearner libraries for stage-1 and stage-2
  treatment models.

- Y.SL.library:

  Character vector. Learners for survivalSL prognostic modeling.

- A.method, Y.method:

  Optional. Scoring metrics passed to `ComputeScores`.

- param.weights.fix, param.weights.init, optim.method, maxit, penalty1,
  penalty2, param.tune:

  Tuning/optimization controls forwarded to `ComputeScores`.

- ngrid:

  Integer. Number of grid points used when integrating survival curves
  for mean survival time.

- censmod:

  Logical. Included for interface consistency; in this wrapper the calls
  to `ComputeScores` set `censmod=FALSE` to compute treatment/prognostic
  (not censoring) scores.

- doublepg:

  Logical. If TRUE, compute `pg0` and `pg1`. If FALSE, compute a single
  `pg`.

- adjustdelta1:

  Logical. If TRUE, define `deltaadj` and use `Y1.var` as the time
  variable for stage-1 scoring; otherwise use `OY.var` and `delta.var`.

- plotps:

  Logical. If TRUE, plots propensity distributions by treatment at each
  stage using [`propensityplot()`](propensityplot.md).

- model.pg:

  Character. Prognostic model family used when `superLearn=FALSE` ("cox"
  or "aft").

- standardize:

  Logical. Whether to standardize covariates for glmnet when
  `superLearn=FALSE`.

- superLearn:

  Logical. If TRUE, use SuperLearner-based estimation inside
  `ComputeScores`; otherwise use glm/glmnet.

- pslink:

  Character. Link for binomial treatment propensity model ("logit" or
  "probit").

- pglink:

  Character. AFT distribution used when `model.pg="aft"` (passed to
  `ComputeScores`).

## Value

A data.frame equal to `data` augmented with score columns. If
`useds=FALSE`, returns `data` unchanged.

**Raw score columns** (merged back by `id.var`):

- Stage 1: `prog01, prog11, prop1` (or `prog01, prop1` if
  `doublepg=FALSE`)

- Stage 2: `prog02, prog12, prop2` (or `prog02, prop2` if
  `doublepg=FALSE`; `NA` for `eta2==0`)

**Standardized columns** (z-scored; propensity on logit scale):

- Stage 1: `pg01, pg11, ps1` (or `pg01, ps1`)

- Stage 2: `pg02, pg12, ps2` (or `pg02, ps2`)

When `doublepg=TRUE`, additional convenience columns are created:

- Stage 1: `pg1ct`, `pg1tc`

- Stage 2: `pg2ct`, `pg2tc`

## Details

**What is computed.**

- A treatment propensity score `ps = P(A=1|Xtrt)` for each stage.

- A prognostic score for each stage based on survival modeling.

If `doublepg=TRUE`, prognostic scores are computed separately under each
treatment level: `pg0` (under `A=0`) and `pg1` (under `A=1`). If
`doublepg=FALSE`, a single prognostic score `pg` is computed.

**Stage 2.** Subjects are subset to `eta2==1` and
[`ComputeScores()`](ComputeScores.md) is called with `Y=Y2.var`,
`A=A2.var`, covariates `names.var2` (prognostic model) and `Xtrt2`
(treatment model). Stage-2 results are merged back to the full dataset;
non-entrants receive `NA` for stage-2 scores.

**Stage 1.** [`ComputeScores()`](ComputeScores.md) is called on the full
cohort with `A=A1.var` and either:

- `Y=OY.var, event=delta.var` if `adjustdelta1=FALSE`, or

- `Y=Y1.var, event='deltaadj'` if `adjustdelta1=TRUE`.

When `adjustdelta1=TRUE`, `deltaadj` is created by copying `delta.var`
and setting `deltaadj=0` for `eta2==1` and `delta==1`.

**Transformations and standardization.** For each stage, `ps` is
transformed using `qlogis(ps)` (logit scale) and then all score columns
are z-scored using [`scale()`](https://rdrr.io/r/base/scale.html) to
produce standardized columns (e.g., `ps1`, `pg01`, `pg11`).

**Convenience contrasts.** When `doublepg=TRUE`, the function creates:

- `pg1ct` / `pg1tc`: stage-1 “correct-treatment” and
  “treatment-contrast” prognostic scores

- `pg2ct` / `pg2tc`: analogous stage-2 versions (for `eta2==1`)

where “correct-treatment” selects `pg0` if observed `A=0` and `pg1` if
observed `A=1`, and “treatment-contrast” selects the opposite arm’s
prognostic score.

**Switch behavior.** If `useds=FALSE`, the function returns `data`
unchanged (no-op), which is useful in pipelines where score construction
is optional.

## See also

[`ComputeScores`](ComputeScores.md),
[`propensityplot`](propensityplot.md)
