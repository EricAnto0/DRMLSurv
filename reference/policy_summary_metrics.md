# Summarize policy performance metrics for a two-stage treatment regime

Computes classification and optional value-based performance metrics for
estimated treatment decisions in a two-stage dynamic treatment regime.
The function compares estimated stage-1 and stage-2 treatment
assignments against observed assignments and returns stage-specific as
well as joint treatment-path metrics.

## Usage

``` r
policy_summary_metrics(
  tmpData,
  estA1,
  estA2,
  obs1_var = "A1.var",
  obs2_var = "A2.var",
  eta2_var = "eta2",
  Tc1_var = NULL,
  Tt1_var = NULL,
  Tc2_var = NULL,
  Tt2_var = NULL,
  Tc_total_var = NULL,
  Tt_total_var = NULL,
  stage = c("both", "stage1", "stage2")
)
```

## Arguments

- tmpData:

  A data frame containing observed treatment assignments, stage-2
  eligibility indicators, and optionally outcome/value columns.

- estA1:

  A vector of estimated stage-1 treatment assignments. Must have length
  equal to `nrow(tmpData)`. Expected coding is `-1` and `1`.

- estA2:

  A vector of estimated stage-2 treatment assignments. Must have length
  equal to `nrow(tmpData)`. Expected coding is `-1` and `1`.

- obs1_var:

  Character string giving the column name in `tmpData` containing
  observed stage-1 treatment assignments. Default is `"A1.var"`.

- obs2_var:

  Character string giving the column name in `tmpData` containing
  observed stage-2 treatment assignments. Default is `"A2.var"`.

- eta2_var:

  Character string giving the column name in `tmpData` indicating
  whether an individual proceeds to stage 2. Subjects with
  `tmpData[[eta2_var]] == 1` are included in stage-2 and joint-path
  metrics. Default is `"eta2"`.

- Tc1_var:

  Optional character string giving the column name for the stage-1 value
  under control or treatment option `-1`. Default is `NULL`.

- Tt1_var:

  Optional character string giving the column name for the stage-1 value
  under treatment option `1`. Default is `NULL`.

- Tc2_var:

  Optional character string giving the column name for the stage-2 value
  under control or treatment option `-1`. Default is `NULL`.

- Tt2_var:

  Optional character string giving the column name for the stage-2 value
  under treatment option `1`. Default is `NULL`.

- Tc_total_var:

  Optional character string giving the column name for the total
  outcome/value under control or treatment option `-1`. Default is
  `NULL`.

- Tt_total_var:

  Optional character string giving the column name for the total
  outcome/value under treatment option `1`. Default is `NULL`.

- stage:

  Character string indicating which stage(s) to summarize. Must be one
  of `"both"`, `"stage1"`, or `"stage2"`. Default is `"both"`.

## Value

A named numeric vector containing some or all of the following
components:

- `Acc1L`:

  Stage-1 accuracy.

- `MCC1L`:

  Stage-1 Matthews correlation coefficient.

- `Sensitivity1L`:

  Stage-1 sensitivity.

- `Specificity1L`:

  Stage-1 specificity.

- `PPV1L`:

  Stage-1 positive predictive value.

- `NPV1L`:

  Stage-1 negative predictive value.

- `F11L`:

  Stage-1 F1 score.

- `RMST1L`:

  Optional stage-1 value summary.

- `Acc2L`:

  Stage-2 accuracy among subjects with `eta2 == 1`.

- `MCC2L`:

  Stage-2 Matthews correlation coefficient.

- `Sensitivity2L`:

  Stage-2 sensitivity.

- `Specificity2L`:

  Stage-2 specificity.

- `PPV2L`:

  Stage-2 positive predictive value.

- `NPV2L`:

  Stage-2 negative predictive value.

- `F12L`:

  Stage-2 F1 score.

- `RMST2L`:

  Optional stage-2 value summary.

- `AccTotal`:

  Joint-path accuracy.

- `MCCTotal`:

  Multiclass Matthews correlation coefficient for the joint path.

- `SensitivityTotalMacro`:

  Macro-averaged sensitivity across joint classes.

- `SpecificityTotalMacro`:

  Macro-averaged specificity across joint classes.

- `PPVTotalMacro`:

  Macro-averaged positive predictive value across joint classes.

- `NPVTotalMacro`:

  Macro-averaged negative predictive value across joint classes.

- `F1TotalMacro`:

  Macro-averaged F1 score across joint classes.

- `SensitivityTotalWtd`:

  Prevalence-weighted sensitivity across joint classes.

- `SpecificityTotalWtd`:

  Prevalence-weighted specificity across joint classes.

- `PPVTotalWtd`:

  Prevalence-weighted positive predictive value across joint classes.

- `NPVTotalWtd`:

  Prevalence-weighted negative predictive value across joint classes.

- `F1TotalWtd`:

  Prevalence-weighted F1 score across joint classes.

- `RMSTTotal`:

  Optional total value summary.

## Details

Stage-1 and stage-2 summaries include accuracy, Matthews correlation
coefficient (MCC), sensitivity, specificity, positive predictive value
(PPV), negative predictive value (NPV), and F1 score. For the joint
two-stage treatment path, the function additionally returns multiclass
MCC together with macro-averaged and prevalence-weighted classification
summaries across the four possible treatment paths.

Optional RMST-style value summaries can also be computed when
corresponding outcome columns are supplied.

The function assumes binary treatment coding with values `-1` and `1`.
Internally, estimated and observed treatments are converted to factors
with levels `c(-1, 1)`.

Stage-1 metrics are computed using all complete cases for the estimated
and observed stage-1 treatment assignments. Stage-2 and joint metrics
are computed only among subjects satisfying `eta2 == 1`.

Joint treatment paths are defined as the four-level factor
`c("-1_-1", "-1_1", "1_-1", "1_1")`, corresponding to the sequence of
stage-1 and stage-2 treatment assignments.

Macro-averaged metrics are calculated as the unweighted mean of
class-wise statistics from the joint confusion matrix. Weighted
summaries use empirical class prevalences from the observed joint
treatment paths.

Optional RMST-style summaries are computed when the corresponding column
names are supplied and found in `tmpData`. These use an internal scoring
rule that assigns subject-level values based on the estimated treatment.

## See also

[`summary.Drmatch()`](https://ericanto0.github.io/DRMLSurv/reference/summary.Drmatch.md),
[`print.summary.Drmatch()`](https://ericanto0.github.io/DRMLSurv/reference/print.summary.Drmatch.md),
[`caret::confusionMatrix()`](https://rdrr.io/pkg/caret/man/confusionMatrix.html),
[`mltools::mcc()`](https://rdrr.io/pkg/mltools/man/mcc.html)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 100
dat <- data.frame(
  A1.var = sample(c(-1, 1), n, replace = TRUE),
  A2.var = sample(c(-1, 1), n, replace = TRUE),
  eta2   = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7))
)

estA1 <- sample(c(-1, 1), n, replace = TRUE)
estA2 <- sample(c(-1, 1), n, replace = TRUE)

policy_summary_metrics(
  tmpData = dat,
  estA1 = estA1,
  estA2 = estA2
)
} # }
```
