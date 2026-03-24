# Compute a value score under a candidate treatment rule

Evaluates the mean outcome under a candidate treatment rule by
assigning, for each subject, the outcome corresponding to the treatment
actually recommended by the rule. If the predicted treatment matches the
observed treatment, the observed outcome under the received treatment is
used; otherwise, the alternative counterfactual or model-based outcome
is used.

## Usage

``` r
my_score.Surv(pred, A, Q, pQ)
```

## Arguments

- pred:

  A vector of predicted treatment assignments under the candidate rule.

- A:

  A vector of observed treatment assignments.

- Q:

  A numeric vector giving the outcome value to use when the predicted
  treatment matches the observed treatment, i.e., for subjects with
  `pred == A`.

- pQ:

  A numeric vector giving the outcome value to use when the predicted
  treatment does not match the observed treatment, i.e., for subjects
  with `pred != A`.

## Value

A single numeric value equal to the mean outcome under the candidate
treatment rule.

## Details

This function is useful in dynamic treatment regime or policy-learning
settings where one wants to estimate the value of a learned treatment
rule using observed and predicted treatment assignments together with
outcome quantities under concordant and discordant treatment choices.

For each subject, the function constructs \$\$ Q_i^\* = \begin{cases}
Q_i, & \text{if } pred_i = A_i, \\ pQ_i, & \text{if } pred_i \neq A_i.
\end{cases} \$\$ and returns the sample mean of \\Q_i^\*\\.

The function assumes that `pred`, `A`, `Q`, and `pQ` are all aligned and
of equal length.

## Examples

``` r
pred <- c(1, -1, 1, -1)
A    <- c(1,  1, 1, -1)
Q    <- c(10, 12,  8, 15)
pQ   <- c( 9, 11,  7, 14)

my_score.Surv(pred = pred, A = A, Q = Q, pQ = pQ)
#> [1] 11
```
