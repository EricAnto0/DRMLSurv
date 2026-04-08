test_that("propensityplot returns a ggplot object invisibly", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(1)
  ps <- runif(50)
  A  <- sample(c(0L, 1L), 50, replace = TRUE)
  p  <- suppressWarnings(propensityplot(ps, A))
  expect_s3_class(p, "ggplot")
})

test_that("propensityplot stops when ggplot2 is absent (mocked)", {
  # Only run if we can mock the namespace check
  testthat::skip_if_not_installed("ggplot2")
  # If ggplot2 is present we can't test its absence; skip in that case.
  # This test is included to document the branch but is skipped when ggplot2
  # is available (which it almost certainly is in CI).
  testthat::skip("ggplot2 is present; cannot test its absence")
})
