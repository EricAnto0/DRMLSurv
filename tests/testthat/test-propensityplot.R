test_that("propensityplot returns a ggplot object invisibly", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(1)
  ps <- runif(50)
  A  <- sample(c(0L, 1L), 50, replace = TRUE)
  p  <- suppressWarnings(propensityplot(ps, A))
  expect_s3_class(p, "ggplot")
})
