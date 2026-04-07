# Compute treatment, prognostic, and optional censoring-related scores

Computes subject-level score summaries for downstream matching,
weighting, or augmentation procedures. The function always estimates a
treatment propensity score and can additionally estimate
treatment-specific prognostic scores and censoring-related scores,
depending on the supplied flags.

Estimation may be carried out using SuperLearner / survivalSL-based
models or via glm / glmnet / Cox / AFT alternatives.

## Usage

``` r
ComputeScores(
  data,
  id,
  Y,
  event,
  X,
  A,
  Xtrt = NULL,
  doublepg = TRUE,
  outer_CV = 5,
  inner_CV = NULL,
  stratifyCV = TRUE,
  cores = 5,
  tau = NULL,
  sl.seed = 100,
  A.SL.library = c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.xgboost"),
  Y.SL.library = c("LIB_COXen", "LIB_AFTggamma"),
  A.method = "method.AUC",
  Y.method = "auc",
  param.tune = NULL,
  ngrid = 2000,
  param.weights.fix = NULL,
  param.weights.init = NULL,
  optim.method = "Nelder-Mead",
  penalty = NULL,
  pgcens = FALSE,
  pscens = TRUE,
  censmod = TRUE,
  maxit = 1000,
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

  A `data.frame` containing all variables required for fitting.

- id:

  Character scalar. Subject identifier column name.

- Y:

  Character scalar. Outcome or follow-up time column name.

- event:

  Character scalar. Event indicator column name, coded `1` for event and
  `0` for censoring.

- X:

  Character vector. Covariate names used in the prognostic model(s).

- A:

  Character scalar. Binary treatment indicator column name. Values may
  be coded as `0/1` or `-1/1`; values equal to `1` are treated as the
  treated group.

- Xtrt:

  Optional character vector. Covariates used in the treatment propensity
  model. If `NULL`, `X` is used.

- doublepg:

  Logical. If `TRUE`, estimate treatment-specific prognostic scores
  `pg0` and `pg1`. If `FALSE`, these columns are returned but remain
  `NA`.

- outer_CV:

  Integer. Number of outer cross-validation folds.

- inner_CV:

  Optional integer. Number of inner cross-validation folds for nested
  SuperLearner fitting.

- stratifyCV:

  Logical. Whether to request stratified cross-validation when supported
  by the underlying fitting routine.

- cores:

  Integer. Number of cores requested for fitting.

- tau:

  Optional numeric truncation horizon used when constructing the
  prediction time grid for restricted mean calculations.

- sl.seed:

  Integer. Random seed used in SuperLearner-based fitting.

- A.SL.library:

  Character vector. SuperLearner library used for treatment propensity
  and censoring-related models.

- Y.SL.library:

  Character vector. Learners used in `survivalSL` for treatment-specific
  prognostic modeling.

- A.method:

  Character scalar. Risk or loss function passed to `CV.SuperLearner()`.

- Y.method:

  Character scalar. Metric passed to `survivalSL()`.

- param.tune:

  Optional tuning object passed to `survivalSL()`.

- ngrid:

  Integer. Number of grid points used for survival-curve prediction and
  numerical integration.

- param.weights.fix:

  Optional vector of fixed ensemble weights passed to `survivalSL()`
  when supported.

- param.weights.init:

  Optional vector of initial ensemble weights passed to `survivalSL()`
  when supported.

- optim.method:

  Character scalar. Optimization method passed to `survivalSL()`.

- penalty:

  Optional penalty value passed to `survivalSL()` or penalized
  regression routines.

- pgcens:

  Logical. If `TRUE` and `censmod = TRUE`, estimate the
  censoring-related prognostic score `pgcens`.

- pscens:

  Logical. If `TRUE` and `censmod = TRUE`, estimate the
  censoring-related propensity score `pscens`.

- censmod:

  Logical. If `TRUE`, request censoring-related scores in addition to
  treatment-based scores.

- maxit:

  Integer. Maximum number of optimization iterations passed to
  `survivalSL()`.

- model.pg:

  Character scalar. Prognostic model family used when
  `superLearn = FALSE`. Must be one of `"cox"` or `"aft"`.

- standardize:

  Logical. Whether to standardize predictors in glmnet-based fits.

- superLearn:

  Logical. If `TRUE`, use SuperLearner / survivalSL-based estimation.
  Otherwise use glm / glmnet / Cox / AFT alternatives.

- pslink:

  Character scalar. Link function for binomial propensity models. Must
  be one of `"logit"` or `"probit"`.

- pglink:

  Character scalar. Distribution used in
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  when `model.pg = "aft"`.

- sl_parallel:

  Character scalar. Parallel mode for SuperLearner-based fitting. Must
  be one of `"multicore"` or `"seq"`.

## Value

A `data.frame` with one row per subject and the following stable
columns:

- `id`: subject identifier,

- `ps`: treatment propensity score,

- `pg0`, `pg1`: treatment-specific prognostic scores,

- `pscens`: censoring-related propensity score,

- `pgcens`: censoring-related prognostic score,

- `ps_sc`, `pg0_sc`, `pg1_sc`, `pscens_sc`, `pgcens_sc`: scaled versions
  of the corresponding scores.

## Details

**Treatment propensity score.**

The function estimates `ps`, the probability of treatment conditional on
`Xtrt` when supplied, or on `X` otherwise.

**Treatment-specific prognostic scores.**

If `doublepg = TRUE`, the function estimates `pg0` and `pg1`,
corresponding to treatment-specific prognostic scores obtained by
fitting separate prognostic models within the observed treatment groups.

**Censoring-related scores.**

If `censmod = TRUE`, the function may additionally estimate `pscens` and
`pgcens`, depending on `pscens` and `pgcens`. These are constructed
using covariates `c(X, A)` and the supplied `event` indicator.

**Scaled outputs.**

The function also returns scaled versions of the raw scores: `ps_sc`,
`pg0_sc`, `pg1_sc`, `pscens_sc`, and `pgcens_sc`. Propensity-type scores
are transformed to the logit scale before standardization.

**Parallel fitting.**

When `superLearn = TRUE`, parallel behavior for SuperLearner-based
estimation is controlled by `sl_parallel`. On Windows, multicore mode is
automatically downgraded to sequential execution.
