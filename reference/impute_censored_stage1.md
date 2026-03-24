# Impute censored stage-1 composite outcomes via constrained donor matching

For each stage-1 censored subject (`death1 == 0`), constructs an
eligible donor pool of uncensored subjects (`death1 == 1`) and imputes
one or more outcome columns (default `compOY`) by matching the censored
subject to `k` donors using MatchIt. Matching can be nearest-neighbor or
optimal, with optional exact matching constraints.

Donor eligibility is enforced by filtering to donors whose primary
outcome `y_cols[1]` is at least the recipient’s observed-time boundary
`OY`. A second filter requires donors to satisfy `is.na(y1) | y1 >= OY`,
which is intended to ensure donors have follow-up long enough relative
to the censored subject.

## Usage

``` r
impute_censored_stage1(
  dat,
  Id,
  exact_vars = NULL,
  formula,
  death1,
  OY,
  y1,
  distance = "mahalanobis",
  method = c("nearest", "optimal"),
  k = 1,
  replace = TRUE,
  caliper = NULL,
  y_cols = c("compOY"),
  aggregate = c("mean", "weighted", "nearest"),
  na_handling = c("drop", "zero_weight")
)
```

## Arguments

- dat:

  A data.frame containing the stage-1 cohort and all variables required
  for donor filtering, matching, and outcome imputation.

- Id:

  Character scalar. Subject identifier column name.

- exact_vars:

  Optional exact matching specification passed to MatchIt. One of:
  `NULL`, a one-sided formula (e.g., `~ sex + site`), or a character
  vector of column names.

- formula:

  A formula giving the matching covariates. Internally
  `update(formula, tr ~ .)` is used to create a two-sided formula with
  `tr` as the matching "treatment" indicator.

- death1:

  Character scalar. Column name for stage-1 event indicator. Convention:
  `1 = observed event/outcome`, `0 = censored`.

- OY:

  Character scalar. Column name for observed-time boundary used to
  constrain donors for each censored subject.

- y1:

  Character scalar. Additional time-like column used in donor filtering:
  donors must satisfy `is.na(y1) | y1 >= OY`.

- distance:

  Character. Distance argument passed to MatchIt (e.g.,
  `"mahalanobis"`).

- method:

  Matching method: `"nearest"` or `"optimal"`.

- k:

  Integer. Donor ratio (number of matched donors per censored subject).

- replace:

  Logical. Whether donors can be reused across matches (nearest-neighbor
  only).

- caliper:

  Optional numeric. Caliper passed to nearest-neighbor matching.

- y_cols:

  Character vector. Outcome column(s) to impute. The first element
  `y_cols[1]` is also used in the donor eligibility filter
  `y_cols[1] >= OY`.

- aggregate:

  Aggregation rule for donor outcomes: `"mean"`, `"weighted"`, or
  `"nearest"`.

- na_handling:

  What to do if the primary imputed column remains missing: `"drop"`
  removes those rows from `data_merged`; `"zero_weight"` keeps them (for
  downstream handling).

## Value

A list with components:

- `imputed`: data.frame of imputed values keyed by `Id` (or `NULL` if
  none imputed),

- `data_merged`: input data with `y_cols` updated for successfully
  imputed IDs,

- `n_censored`: number of censored subjects targeted,

- `n_imputed`: number successfully imputed,

- `errors`: two-column matrix recording IDs and error messages from
  MatchIt.

## Details

The function loops over censored IDs. For each censored ID, it forms a
temporary dataset consisting of the focal censored unit plus its
eligible donors, and runs MatchIt by defining `tr = (death1 == 0)`
(censored-as-treated). Donors are taken from the same `subclass` as the
focal censored unit.

Donor outcomes are aggregated to produce the imputation:

- `"mean"`: arithmetic mean across donors,

- `"weighted"`: weighted mean using `weights` returned by MatchIt,

- `"nearest"`: take the single nearest donor by `distance`, even if
  `k > 1`.

If an ID fails to match (insufficient donors, matching error, etc.), it
is skipped and its imputation remains missing unless handled via
`na_handling`.

## See also

[`matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html),
[`get_matches`](https://kosukeimai.github.io/MatchIt/reference/match_data.html)
