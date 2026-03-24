# Summarize a fitted `Drmatch` object on new data

Generates predicted treatment decisions from a fitted `Drmatch` object
and summarizes their performance against observed treatment assignments
in a new dataset. Depending on the requested stage, the function
evaluates stage-1 decisions, stage-2 decisions, or both jointly.

## Usage

``` r
# S3 method for class 'Drmatch'
summary(
  object,
  newdata,
  stage = c("both", "stage1", "stage2"),
  obs1_var = NULL,
  obs2_var = NULL,
  eta2_var = NULL,
  pred1_var = "A1.opt",
  pred2_var = "A2.opt",
  Tc1_var = NULL,
  Tt1_var = NULL,
  Tc2_var = NULL,
  Tt2_var = NULL,
  Tc_total_var = NULL,
  Tt_total_var = NULL,
  as.data.frame = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted object of class `"Drmatch"`.

- newdata:

  A data frame on which to compute predicted treatment assignments and
  performance summaries.

- stage:

  Character string indicating which decision stage to summarize. One of
  `"both"`, `"stage1"`, or `"stage2"`. Default is `"both"`.

- obs1_var:

  Optional character string giving the observed stage-1 treatment column
  in `newdata`. If `NULL`, the value is taken from `object$A1.var`.

- obs2_var:

  Optional character string giving the observed stage-2 treatment column
  in `newdata`. If `NULL`, the value is taken from `object$A2.var`.

- eta2_var:

  Optional character string giving the stage-2 eligibility column in
  `newdata`. If `NULL`, the value is taken from `object$eta2.var`.

- pred1_var:

  Character string giving the column name in the prediction output
  corresponding to the estimated stage-1 treatment. Default is
  `"A1.opt"`.

- pred2_var:

  Character string giving the column name in the prediction output
  corresponding to the estimated stage-2 treatment. Default is
  `"A2.opt"`.

- Tc1_var:

  Optional character string giving the stage-1 value column for
  treatment option `-1`.

- Tt1_var:

  Optional character string giving the stage-1 value column for
  treatment option `1`.

- Tc2_var:

  Optional character string giving the stage-2 value column for
  treatment option `-1`.

- Tt2_var:

  Optional character string giving the stage-2 value column for
  treatment option `1`.

- Tc_total_var:

  Optional character string giving the total value column for treatment
  option `-1`.

- Tt_total_var:

  Optional character string giving the total value column for treatment
  option `1`.

- as.data.frame:

  Logical; if `TRUE`, the output is returned as a one-row data frame.
  Otherwise, a named numeric vector is returned. Default is `TRUE`.

- ...:

  Additional arguments passed to
  [`predict.Drmatch()`](predict.Drmatch.md).

## Value

Either a one-row data frame or a named numeric vector containing the
policy performance metrics returned by
[`policy_summary_metrics()`](policy_summary_metrics.md).

## Details

Internally, predicted treatment assignments are passed to
[`policy_summary_metrics()`](policy_summary_metrics.md) to compute
stage-specific and joint-path classification metrics, with optional
value-based summaries when outcome columns are supplied.

When `stage = "stage1"`, stage-2 predictions are set to `NA` before
summary metrics are computed. When `stage = "stage2"`, the observed
stage-1 treatment from `newdata` is used as the stage-1 component of the
joint path, while stage-2 predictions are taken from
`predict(object, ...)`.

This method requires `newdata` because performance is evaluated by
comparing predicted treatment decisions with observed treatments and,
optionally, with value columns provided in the new dataset.

## See also

[`policy_summary_metrics()`](policy_summary_metrics.md),
[`print.summary.Drmatch()`](print.summary.Drmatch.md),
[`predict.Drmatch()`](predict.Drmatch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# fit <- Drmatch(...)
# summary(fit, newdata = test_dat)
# summary(fit, newdata = test_dat, stage = "stage1")
# summary(fit, newdata = test_dat, stage = "stage2", as.data.frame = FALSE)
} # }
```
