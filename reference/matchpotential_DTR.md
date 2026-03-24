# Construct matched counterfactual outcomes for a (single-stage) treatment contrast

Builds matched “potential outcomes” for each unit under both treatment
levels by matching each unit to comparable units in the opposite
treatment group, and then imputing the opposite-arm outcome using
matched donor outcomes. The function supports either:

- **full matching** (one global match structure), or

- **nearest-neighbor matching** via a two-pass strategy (ATT then ATC)
  so that every unit receives an opposite-arm imputation while
  respecting anchor definitions.

The returned dataset includes (i) the observed outcome, (ii) an imputed
opposite-arm “paired” outcome, and (iii) constructed counterfactual
outcomes for each arm for every unit.

## Usage

``` r
matchpotential_DTR(
  dat,
  txgroup,
  exact_vars,
  compY,
  vec,
  Id,
  method = c("nearest", "full"),
  k = 3,
  replace = TRUE,
  caliper = NULL,
  distance = "mahalanobis",
  compW = "ipcw.R",
  na_handling = c("drop", "zero_weight")
)
```

## Arguments

- dat:

  A data.frame containing at minimum the treatment indicator `txgroup`,
  outcome `compY`, identifier `Id`, and matching covariates in `vec` (or
  `vec[[1]]`, `vec[[2]]` for nearest).

- txgroup:

  Character scalar. Name of the treatment indicator column. Internally
  recoded to \\\\-1,+1\\\\.

- exact_vars:

  Optional. Exact matching specification passed to MatchIt. Can be:
  `NULL`, a one-sided formula (e.g., `~ sex + site`), or a character
  vector of column names.

- compY:

  Character scalar. Outcome column used to construct counterfactuals and
  paired outcomes.

- vec:

  Matching covariates. For `method="full"`, `vec` should be a character
  vector of covariate names. For `method="nearest"`, your current
  implementation expects `vec` to be a list with `vec[[1]]` (covariates
  used when treated are anchors) and `vec[[2]]` (covariates used when
  controls are anchors).

- Id:

  Character scalar. Subject identifier column name.

- method:

  Matching method: `"nearest"` (two-pass ATT/ATC) or `"full"`.

- k:

  Integer. Match ratio (number of donors per anchor) for
  nearest-neighbor matching. Ignored for full matching.

- replace:

  Logical. Whether matching is performed with replacement
  (nearest-neighbor only).

- caliper:

  Optional numeric. Caliper passed to
  [`MatchIt::matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
  (applies to both matching types when supplied).

- distance:

  Distance specification passed to
  [`MatchIt::matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
  (e.g., `"mahalanobis"` or a numeric vector).

- compW:

  Character scalar. Name of the “weight-like” column used to compute
  `diff.wt`. If missing, it is created as a copy of `compY`. Default is
  `"ipcw.R"`.

- na_handling:

  How to handle units with missing paired outcomes/weights: `"drop"` or
  `"zero_weight"`.

## Value

A data.frame containing (at least) the original columns and additional
derived columns:

- `pairedCompY`: imputed opposite-arm outcome for each unit (mean within
  matched set)

- `paired.ipcw.R`: imputed opposite-arm `compW` value (mean within
  matched set)

- `trt_cf1L`, `ctrl_cf1L`: constructed counterfactual outcomes under
  treatment/control

- `diff.wt`, `match.weight`: discrepancy in `compW` and its absolute
  value

- `newTxt`: sign-coded label based on treatment and `diff.wt`

## Details

**Treatment coding.** Internally, the treatment indicator `txgroup` is
recoded to \\\\-1,+1\\\\, where `+1` denotes treated and `-1` denotes
control. Matching formulas are constructed using `I(tx == 1)` or
`I(tx == -1)` to define anchors in the nearest-neighbor passes.

**Outcome and weight variables.**

- `compY` is the (possibly composite) observed outcome used to form
  counterfactuals.

- `compW` is an auxiliary “weight-like” variable (default `"ipcw.R"`)
  used to compute a match quality/weight proxy. If `compW` is not
  present in the data, it is initialized as `compW <- compY` (so weights
  default to outcome differences).

**Matching covariates.** The matching covariates are supplied via `vec`
and used to build:

- `fmla_bin`: full matching formula using all covariates in `vec`

- `fmla_bin_trt`: nearest-neighbor pass 1 formula (treated anchors)

- `fmla_bin_ctrl`: nearest-neighbor pass 2 formula (control anchors)

In your current implementation, `fmla_bin_trt` uses only `vec[[1]]` and
`fmla_bin_ctrl` uses only `vec[[2]]`. This implies you are allowing
different covariate sets for the two passes (which can be intentional),
but you should ensure `vec` is structured as a list with elements
`vec[[1]]` and `vec[[2]]` if you rely on this behavior.

**Full matching branch.** Runs
`MatchIt::matchit(..., method="full", estimand="ATE")` and then, within
each matched `subclass`, computes the mean outcome in treated and
control units. Each unit’s opposite-arm imputation is the mean outcome
in the opposite arm within its subclass. The same is done for the
weight-like variable `compW`.

**Nearest-neighbor branch (two-pass).** To obtain an opposite-arm
imputation for every unit (treated and control) while keeping anchors
consistent:

1.  *ATT pass (treated anchors).* Match treated units to controls
    (`estimand="ATT"`). For each treated anchor, compute the mean
    control outcome (and mean `compW`) among its matched controls;
    assign these as `pairedCompY_ATT_for_treated` and
    `paired.ipcw.R_ATT_for_treated`.

2.  *ATC pass (control anchors).* Re-run matching with control anchors
    by switching the “treated” definition in the formula to `I(tx==-1)`.
    For each control anchor, compute the mean treated outcome (and mean
    `compW`) among its matched treated units; assign these as
    `pairedCompY_ATC_for_controls` and `paired.ipcw.R_ATC_for_controls`.

3.  Merge the two summaries back to a one-row-per-ID dataset, and choose
    the appropriate opposite-arm imputation depending on the unit’s
    observed treatment.

**Constructed counterfactuals.** After matching, the function constructs
two counterfactual outcome columns for each unit:

- `trt_cf1L`: the unit’s outcome under treatment. Equals observed
  `compY` if treated, otherwise equals the matched opposite-arm
  imputation.

- `ctrl_cf1L`: the unit’s outcome under control. Equals observed `compY`
  if control, otherwise equals the matched opposite-arm imputation.

**Match-weight proxy.** A simple distance-like proxy is computed as
`diff.wt = compW - paired.ipcw.R` and `match.weight = abs(diff.wt)`.
This is not a MatchIt weight; it is a user-defined measure intended to
quantify discrepancy between the unit and its opposite-arm donor set in
`compW`.

**Handling missing opposite-arm matches.** If a unit fails to obtain a
valid opposite-arm imputation (e.g., due to empty donor sets), behavior
is controlled by `na_handling`:

- `"drop"`: remove affected units from the output.

- `"zero_weight"`: keep them but set `match.weight = 0`.

## See also

[`matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html),
[`match.data`](https://kosukeimai.github.io/MatchIt/reference/match_data.html),
[`get_matches`](https://kosukeimai.github.io/MatchIt/reference/match_data.html)
