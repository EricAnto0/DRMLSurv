# Tune and fit a random-forest DTR policy model with (optional) cross-validation

Fits a classification-based dynamic treatment regime (DTR) policy model
for a binary treatment learns a binary decision rule \\\hat d(X)\\ with
values in \\\\-1,+1\\\\ using either ranger or randomForestSRC. The
function performs grid search over random-forest tuning parameters and,
when `usecv=TRUE`, evaluates candidate policies via K-fold
cross-validation using:

- **CCR**: classification correctness rate (accuracy of predicted
  treatment labels),

- **OOB**: out-of-bag prediction error (model-reported), and

- **Score**: a user-supplied policy-value-like performance measure
  computed by
  [`my_score.Surv()`](https://ericanto0.github.io/DRMLSurv/reference/my_score.Surv.md)
  using observed and matched pseudo-outcomes.

The best tuning parameters are selected according to `metric`, after
which a final model is fit on the full observed dataset `obs`. The
output includes the final fitted model, the estimated treatment rule
`estA.obs`, and the full tuning results.

## Usage

``` r
rfdtr(
  modeltype = "ranger",
  usecv = TRUE,
  sl.seed = 123,
  obs,
  W,
  gridpar,
  metric = "ccr",
  A.obs,
  Q.obs,
  Q.match,
  score_agg = c("mean", "sum")
)
```

## Arguments

- modeltype:

  Character. Random-forest engine: `"ranger"` or `"rfsrc"`.

- usecv:

  Logical. If TRUE, uses 5-fold CV to tune hyperparameters; otherwise
  tunes on the full data.

- sl.seed:

  Integer. Seed for reproducibility. (Note: parts of the function
  currently hard-code `set.seed(123)`.)

- obs:

  data.frame. Training dataset containing a column `A` and covariates
  for predicting `A`.

- W:

  Numeric vector. Case weights of length `nrow(obs)`.

- gridpar:

  data.frame. Grid of tuning parameters. Expected columns include:
  `ntree`, `mtry`, and `nodesize`. Extra columns are ignored by the
  model fits.

- metric:

  Character. Criterion used to choose the “best” tuning row. Supported:

  - `"oob"`: minimize `OOB`

  - `"policyval"` / `"score"` / `"policy"` / `"val"`: maximize `Score`

  - otherwise: maximize `CCR`

- A.obs:

  Numeric vector. Observed treatment labels (`-1/1`) aligned with `obs`.

- Q.obs:

  Numeric vector. Observed pseudo-outcome or value component aligned
  with `obs`.

- Q.match:

  Numeric vector. Matched/pair pseudo-outcome aligned with `obs` used by
  [`my_score.Surv()`](https://ericanto0.github.io/DRMLSurv/reference/my_score.Surv.md).

- score_agg:

  Character. Aggregation for fold-level `Score`: `"sum"` or `"mean"`.

## Value

A list with elements:

- `model`: final fitted ranger or randomForestSRC object trained on all
  `obs`.

- `estA.obs`: numeric vector of predicted treatment decisions (`-1/1`)
  for each row of `obs`.

- `tune`: data.frame of tuning results for each row of `gridpar`,
  including `CCR`, `OOB`, `Score`.

- `best`: the selected row of `tune` corresponding to the chosen
  `metric`.

## Details

**Inputs and data layout.**

- `obs` is a data.frame that must contain a column named `A` and the
  covariates used to model `A`. The function coerces `obs$A` to a factor
  with levels `c(-1, 1)` to ensure a binary classification target.

- `W` is a numeric vector of case weights aligned to `obs` (e.g., IPC
  weights).

- `A.obs` is the numeric observed treatment indicator (typically `-1/1`)
  aligned to `obs`.

- `Q.obs` and `Q.match` are numeric vectors aligned to `obs` used by
  [`my_score.Surv()`](https://ericanto0.github.io/DRMLSurv/reference/my_score.Surv.md)
  to compute the policy score (e.g., observed pseudo-outcome and
  matched/pair pseudo-outcome).

**Cross-validation logic.** When `usecv=TRUE`, the function uses 5-fold
stratified folds created by `caret::createFolds(obs$A, k=5)`. For each
candidate parameter set in `gridpar`, it fits the model on the training
folds and evaluates on the held-out fold:

- *CCR*: mean(`predicted_class == obs$A[test]`).

- *OOB*: model-reported OOB error (for ranger, `mod$prediction.error`;
  for rfsrc, the last entry of `mod$err.rate[,"all"]`).

- *Score*: `my_score.Surv(pred, A.test, Q.test, Q.match.test)` where
  `pred` is the numeric treatment rule `-1/1`.

**Aggregation of score across folds.** The per-fold `Score` values are
aggregated across CV folds using:

- `score_agg="sum"`: sum of fold scores (ignoring `NA`)

- `score_agg="mean"`: mean of fold scores (ignoring `NA`)

This controls the scale used for model selection when `metric` targets
policy value.

**Parameter guards.** To avoid invalid random-forest hyperparameters, if
`gridpar` includes `mtry`, it is clamped to `[1, p]` where
`p = ncol(obs)-1` (i.e., all columns except `A`). Duplicate parameter
rows are removed via `unique(gridpar)`.

**Parallelization.** The function is written using foreach with
`%dopar%`. It assumes a parallel backend has already been registered
(e.g., via doParallel). Within each fold, `num.threads=1` is used for
ranger to avoid nested parallelism.

## See also

[`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md),
[`rfsrc`](https://www.randomforestsrc.org//reference/rfsrc.html),
[`createFolds`](https://rdrr.io/pkg/caret/man/createDataPartition.html)
