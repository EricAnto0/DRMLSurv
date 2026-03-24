# Plot propensity (or censoring) score overlap by group

Produces a diagnostic overlap plot for a set of propensity-like scores
(e.g., treatment propensity scores or censoring propensities) stratified
by a binary group indicator `A`. The function visualizes the empirical
score distributions using semi-transparent histograms on a density
scale.

This plot is primarily intended to assess common support / overlap and
to diagnose separation or extreme predicted probabilities before
matching/weighting steps.

## Usage

``` r
propensityplot(ps, A)
```

## Arguments

- ps:

  Numeric vector of propensity-like scores, typically in \\\[0,1\]\\.

- A:

  Vector defining the grouping variable (e.g., treatment arm or event
  indicator). Will be coerced to a factor for plotting.

## Value

A ggplot2 plot object (class `"gg"` and `"ggplot"`). The plot is also
printed.

## Details

The input `A` is coerced to a factor and used for coloring and filling
the histogram. The y-axis is scaled to density (`..density..`). The
function returns a ggplot2 object (invisibly) and also prints the plot
as a side effect, which is convenient in interactive use.

**Package note.** In a package context, avoid
[`library(ggplot2)`](https://ggplot2.tidyverse.org) inside functions.
Instead, use `ggplot2::` calls and guard availability via
[`requireNamespace("ggplot2", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html).

## See also

[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html),
[`geom_histogram`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
