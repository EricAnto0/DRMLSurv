# Compute propensity and prognostic scores (treatment- or censoring-based)

Computes (i) a propensity score and (ii) prognostic score(s) for use in
matching/weighting pipelines. Supports SuperLearner-based fitting or
glm/glmnet-based alternatives. Can compute treatment PS +
treatment-specific prognostic scores (pg0/pg1), or censoring PS/PG
(pscens/pgcens) when `censmod=TRUE`.

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
  stratifyCV = FALSE,
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
  pglink = "lognormal"
)
```

## Arguments

- data:

  data.frame with all variables.

- id:

  Character scalar. Subject ID column name.

- Y:

  Character scalar. Survival time column name.

- event:

  Character scalar. Event indicator column name (1=event, 0=censored).

- X:

  Character vector. Covariate column names for prognostic model(s).

- A:

  Character scalar. Treatment indicator column name (coded 0/1 or -1/1).

- Xtrt:

  Optional character vector. Covariates for treatment model if different
  from `X`.

- doublepg:

  Logical. If TRUE and `censmod=FALSE`, fit separate prognostic models
  by treatment arm.

- outer_CV:

  Integer. Outer CV folds.

- inner_CV:

  Optional integer. Inner CV folds for nested CV (SuperLearner).

- stratifyCV:

  Logical. Whether to stratify CV folds.

- cores:

  Integer. Requested cores for parallel fit (multicore only on
  non-Windows).

- tau:

  Optional numeric. Truncation horizon for mean survival.

- sl.seed:

  Integer. Seed for SuperLearner.

- A.SL.library:

  Character vector. SL learners for propensity/censoring models.

- Y.SL.library:

  Character vector. Learners for survivalSL prognostic models.

- A.method:

  Character. CV risk method for propensity SL.

- Y.method:

  Character. Metric for survivalSL.

- param.tune, ngrid, param.weights.fix, param.weights.init,
  optim.method, penalty, maxit:

  Control survivalSL fitting and prediction.

- pgcens, pscens, censmod:

  Logical flags controlling censoring scores.

- model.pg:

  Character. "cox" or "aft" for non-SL prognostic modeling.

- standardize:

  Logical. Standardize covariates for glmnet.

- superLearn:

  Logical. Use SuperLearner/survivalSL branches if TRUE.

- pslink:

  Character. "logit" or "probit".

- pglink:

  Character. AFT distribution for flexsurvreg.

## Value

A data.frame with stable columns:

- `id`: ID values

- `ps`: treatment propensity score (if computed)

- `pg0`, `pg1`: treatment-specific prognostic scores (if computed)

- `pscens`: censoring propensity score (if computed)

- `pgcens`: censoring prognostic score (if computed)
