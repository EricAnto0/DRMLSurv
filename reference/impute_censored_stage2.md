# Impute censored stage-2 outcomes via constrained donor matching

For each stage-2 censored subject (`delta.var == 0`), constructs an
eligible donor pool of uncensored subjects (`delta.var == 1`) and
imputes the stage-2 outcome `Y2.var`, returning the imputed/observed
outcome as `compY2`. Matching is performed using MatchIt with
nearest-neighbor or optimal matching, optional exact matching, and
configurable aggregation across `k` donors.

Donor eligibility is constrained by requiring donors to have observed
overall time `OY.var` at least as large as the recipient’s `OY.var`.
This is intended to ensure donors have follow-up long enough relative to
the censored subject.

## Usage

``` r
impute_censored_stage2(
  dat,
  id.var,
  delta.var,
  OY.var,
  Y2.var,
  formula2,
  exact.vars = NULL,
  method = c("nearest", "optimal"),
  distance = "mahalanobis",
  k = 1,
  replace = TRUE,
  caliper = NULL,
  aggregate = c("mean", "weighted", "nearest")
)
```

## Arguments

- dat:

  A data.frame containing stage-2 entrants and required variables for
  matching and filtering.

- id.var:

  Character scalar. Subject identifier column name.

- delta.var:

  Character scalar. Stage-2 event indicator column name. Convention:
  `1 = observed event/outcome`, `0 = censored`.

- OY.var:

  Character scalar. Column name for observed-time boundary used to
  constrain donors.

- Y2.var:

  Character scalar. Stage-2 outcome column to impute.

- formula2:

  A formula giving matching covariates. Internally
  `update(formula2, tr ~ .)` is used to create a two-sided formula with
  `tr` as the matching "treatment" indicator.

- exact.vars:

  Optional exact matching specification passed to MatchIt: `NULL`, a
  one-sided formula, or a character vector of names.

- method:

  Matching method: `"nearest"` or `"optimal"`.

- distance:

  Character. Distance argument passed to MatchIt (e.g.,
  `"mahalanobis"`).

- k:

  Integer. Donor ratio (number of matched donors per censored subject).

- replace:

  Logical. Whether donors can be reused across matches (nearest-neighbor
  only).

- caliper:

  Optional numeric. Caliper for nearest-neighbor matching.

- aggregate:

  Aggregation rule: `"mean"`, `"weighted"`, or `"nearest"`.

## Value

A list with components:

- `imputed`: data.frame of imputed `compY2` values keyed by `id.var` (or
  `NULL` if none),

- `data_merged`: the input `dat` augmented with `compY2`,

- `n_censored`: number of censored subjects targeted,

- `n_imputed`: number successfully imputed,

- `errors`: two-column matrix recording IDs and error messages from
  MatchIt.

## Details

The function loops over censored IDs. For each ID, it forms a temporary
dataset consisting of the focal censored unit plus eligible donors and
runs MatchIt by defining `tr = (delta.var == 0)` (censored-as-treated).
Donors are then taken from the matched `subclass` containing the focal
censored unit.

Donor outcomes `Y2.var` are aggregated according to `aggregate`:

- `"mean"`: arithmetic mean,

- `"weighted"`: weighted mean using `weights` from MatchIt,

- `"nearest"`: take the single nearest donor by `distance`.

For uncensored subjects (`delta.var == 1`), `compY2` is set to the
observed `Y2.var`.

## See also

[`matchit`](https://kosukeimai.github.io/MatchIt/reference/matchit.html),
[`get_matches`](https://kosukeimai.github.io/MatchIt/reference/match_data.html)
