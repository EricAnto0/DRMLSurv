# DRMLSurv: Doubly Robust ML and Matching for Two-Stage Survival Outcomes

Implements a two-stage survival analysis pipeline that combines (i)
imputation of censored stage-specific survival times via donor matching,
(ii) construction of matched counterfactual outcomes under static and
dynamic treatment regimes, and (iii) estimation of optimal dynamic
treatment rules using machine learning and cross-fitting. The toolkit
supports nearest and full matching, optional exact matching, balance
diagnostics, and flexible censoring/propensity score models via
DoubleScore/SuperLearner-style learners.

## See also

Useful links:

- <https://ericanto0.github.io/DRMLSurv/>

## Author

**Maintainer**: Eric Anto <eric.anto@utah.edu>
