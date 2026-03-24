# Predict optimal treatment decisions from a fitted Drmatch object

Generates predicted optimal treatment assignments from a fitted
`Drmatch` object for new observations using stage-specific random forest
models fitted with the ranger package. Predictions may be obtained for
stage 1 only, stage 2 only, or both stages of a two-stage treatment
regime.

## Usage

``` r
# S3 method for class 'Drmatch'
predict(object, newdata, stage = c("both", "stage1", "stage2"), ...)
```

## Arguments

- object:

  A fitted `Drmatch` object containing the trained stage-specific ranger
  random forest models and metadata needed for deployment, including
  `stage1_model`, `stage2_model`, `names.var1`, `names.var2`, and
  `eta2.var`.

- newdata:

  A data frame containing the predictor variables required for stage-1
  and/or stage-2 prediction.

- stage:

  Character string indicating which stage predictions to return. Must be
  one of `"both"`, `"stage1"`, or `"stage2"`. The default is `"both"`.

- ...:

  Additional arguments passed through for S3 compatibility.

## Value

A data frame with one row per row of `newdata`. The returned data frame
always contains:

- `row_id`:

  Row index corresponding to the original row position in `newdata`.

Depending on `stage`, it also contains:

- `A1.opt`:

  Predicted optimal stage-1 treatment, or `NA` when stage-1 prediction
  is not available for that row.

- `A2.opt`:

  Predicted optimal stage-2 treatment, or `NA` when stage-2 prediction
  is not available for that row.

## Details

For stage 1, predictions are produced for rows in `newdata` with
complete values on the stage-1 predictor set stored in
`object$names.var1`. For stage 2, predictions are produced only for rows
that are eligible for second-stage treatment, as indicated by
`object$eta2.var == 1`, and that also have complete values on the
stage-2 predictor set stored in `object$names.var2`.

Rows that are ineligible for prediction because of missing required
covariates, or because they are not eligible for stage 2, receive `NA`
for the corresponding predicted treatment.

This function is an S3
[`predict()`](https://rdrr.io/r/stats/predict.html) method for objects
of class `Drmatch`.

When `stage = "stage1"`, the function checks that all variables listed
in `object$names.var1` are present in `newdata`. Predictions are
returned in the column `A1.opt`.

When `stage = "stage2"`, the function checks that all variables listed
in `object$names.var2` are present in `newdata`, and also verifies that
the stage-2 eligibility indicator named by `object$eta2.var` is
available in `newdata`. Predictions are returned in the column `A2.opt`.

When `stage = "both"`, both `A1.opt` and `A2.opt` are returned.

Internally, stage-specific predictions are obtained from fitted ranger
models using the `predictions` component returned by `predict.ranger()`:

    predict(model, data = newdata_subset)$predictions

For classification forests, these predicted class labels are coerced to
numeric values, so the fitted models are expected to predict treatment
classes coded as `-1` and `1`.

## Required components of `object`

The fitted `Drmatch` object is expected to contain at least the
following elements:

- `stage1_model`:

  A fitted ranger classification model for stage 1.

- `stage2_model`:

  A fitted ranger classification model for stage 2.

- `names.var1`:

  Character vector of predictor names required for stage 1.

- `names.var2`:

  Character vector of predictor names required for stage 2.

- `eta2.var`:

  Name of the stage-2 eligibility indicator in `newdata`.

## Examples

``` r
## Not run:
## Suppose `fit` is a fitted Drmatch object
## and `new_patients` is a data frame of candidate patients.
##
## Predict both stages
## pred <- predict(fit, newdata = new_patients, stage = "both")
##
## Predict stage 1 only
## pred1 <- predict(fit, newdata = new_patients, stage = "stage1")
##
## Predict stage 2 only
## pred2 <- predict(fit, newdata = new_patients, stage = "stage2")
```
