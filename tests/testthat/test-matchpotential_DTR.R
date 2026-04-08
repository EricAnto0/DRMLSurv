make_match_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  data.frame(
    id  = seq_len(n),
    tx  = sample(c(0L, 1L), n, replace = TRUE),
    x1  = rnorm(n),
    x2  = rnorm(n),
    y   = abs(rnorm(n, mean = 10, sd = 2)),
    w   = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )
}

test_that("matchpotential_DTR stops when required columns are missing", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = NULL,
      compY = "MISSING_Y", vec = "x1", Id = "id",
      method = "full"
    ),
    "Missing required columns"
  )
})

test_that("matchpotential_DTR stops when Id column is missing", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = NULL,
      compY = "y", vec = "x1", Id = "MISSING_ID",
      method = "full"
    ),
    "Missing required columns"
  )
})

test_that("matchpotential_DTR stops when vec covariate missing (full method)", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = NULL,
      compY = "y", vec = "MISSING_COV", Id = "id",
      method = "full"
    ),
    "covariates not found"
  )
})

test_that("matchpotential_DTR stops when vec is not character for full method", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = NULL,
      compY = "y", vec = list("x1"), Id = "id",
      method = "full"
    ),
    "vec must be a character vector"
  )
})

test_that("matchpotential_DTR stops when vec is not a list for nearest method", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = NULL,
      compY = "y", vec = "x1", Id = "id",
      method = "nearest"
    ),
    "vec must be a list"
  )
})

test_that("matchpotential_DTR stops when vec list too short for nearest method", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = NULL,
      compY = "y", vec = list("x1"), Id = "id",
      method = "nearest"
    ),
    "vec must be a list"
  )
})

test_that("matchpotential_DTR stops on invalid exact_vars type", {
  dat <- make_match_data()
  expect_error(
    matchpotential_DTR(
      dat, txgroup = "tx", exact_vars = 123L,
      compY = "y", vec = "x1", Id = "id",
      method = "full"
    ),
    "exact_vars must be NULL"
  )
})

test_that("matchpotential_DTR full method returns expected columns", {
  dat <- make_match_data(n = 80, seed = 7)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full", distance = "mahalanobis"
  )
  expect_true(is.data.frame(res))
  expect_true("pairedCompY" %in% names(res))
  expect_true("paired.ipcw.R" %in% names(res))
  expect_true("trt_cf1L" %in% names(res))
  expect_true("ctrl_cf1L" %in% names(res))
  expect_true("diff.wt" %in% names(res))
  expect_true("match.weight" %in% names(res))
  expect_true("newTxt" %in% names(res))
})

test_that("matchpotential_DTR full method with character exact_vars works", {
  dat <- make_match_data(n = 80, seed = 8)
  # add a binary exact-matching variable
  dat$grp <- sample(c("A","B"), nrow(dat), replace = TRUE)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = "grp",
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full", distance = "mahalanobis"
  )
  expect_true(is.data.frame(res))
})

test_that("matchpotential_DTR full method with formula exact_vars works", {
  dat <- make_match_data(n = 80, seed = 9)
  dat$grp <- sample(c("A","B"), nrow(dat), replace = TRUE)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = ~ grp,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full", distance = "mahalanobis"
  )
  expect_true(is.data.frame(res))
})

test_that("matchpotential_DTR nearest method returns expected columns", {
  dat <- make_match_data(n = 80, seed = 10)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y",
    vec = list(c("x1", "x2"), c("x1", "x2")),
    Id = "id",
    method = "nearest", k = 2, replace = TRUE,
    distance = "mahalanobis"
  )
  expect_true(is.data.frame(res))
  expect_true("trt_cf1L"  %in% names(res))
  expect_true("ctrl_cf1L" %in% names(res))
})

test_that("matchpotential_DTR uses compW fallback when column absent", {
  dat <- make_match_data(n = 80, seed = 11)
  # ipcw.R column does not exist -> should fall back to compY
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full", compW = "ipcw.R"
  )
  expect_true(is.data.frame(res))
})

test_that("matchpotential_DTR na_handling='zero_weight' keeps bad rows", {
  # Use full-matching scenario where we can observe the behavior
  dat <- make_match_data(n = 60, seed = 12)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full", na_handling = "zero_weight"
  )
  expect_true(is.data.frame(res))
})

test_that("matchpotential_DTR match.weight is non-negative", {
  dat <- make_match_data(n = 80, seed = 13)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full"
  )
  expect_true(all(res$match.weight >= 0, na.rm = TRUE))
})

test_that("matchpotential_DTR trt_cf1L equals observed y for treated", {
  dat <- make_match_data(n = 80, seed = 14)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full"
  )
  # internally tx is recoded to {-1,+1}; find treated (original tx==1)
  orig_trt <- dat$tx[match(res$id, dat$id)] == 1
  # for treated, trt_cf1L should equal their observed outcome y
  expect_equal(
    res$trt_cf1L[orig_trt],
    res$y[orig_trt]
  )
})

test_that("matchpotential_DTR ctrl_cf1L equals observed y for controls", {
  dat <- make_match_data(n = 80, seed = 15)
  res <- matchpotential_DTR(
    dat, txgroup = "tx", exact_vars = NULL,
    compY = "y", vec = c("x1", "x2"), Id = "id",
    method = "full"
  )
  orig_ctrl <- dat$tx[match(res$id, dat$id)] == 0
  expect_equal(
    res$ctrl_cf1L[orig_ctrl],
    res$y[orig_ctrl]
  )
})
