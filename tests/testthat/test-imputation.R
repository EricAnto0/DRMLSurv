make_imputation_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  data.frame(
    id     = seq_len(n),
    delta  = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.3, 0.7)),
    OY     = sort(runif(n, 1, 24)),
    Y2     = abs(rnorm(n, 6, 2)),
    x1     = rnorm(n),
    x2     = rnorm(n),
    stringsAsFactors = FALSE
  )
}

make_stage1_data <- function(n = 60, seed = 2) {
  set.seed(seed)
  data.frame(
    id      = seq_len(n),
    death1  = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.3, 0.7)),
    compOY  = abs(rnorm(n, 12, 4)),
    OY      = sort(runif(n, 1, 24)),
    Y1      = abs(rnorm(n, 5, 2)),
    x1      = rnorm(n),
    x2      = rnorm(n),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# impute_censored_stage2 tests
# ---------------------------------------------------------------------------

test_that("impute_censored_stage2 stops when required columns are missing", {
  dat <- make_imputation_data()
  expect_error(
    impute_censored_stage2(
      dat,
      id.var    = "id",
      delta.var = "MISSING",
      OY.var    = "OY",
      Y2.var    = "Y2",
      formula2  = ~ x1 + x2
    ),
    "Missing required columns"
  )
})

test_that("impute_censored_stage2 stops when data is not a data.frame", {
  expect_error(
    impute_censored_stage2(
      list(id = 1),
      id.var    = "id",
      delta.var = "delta",
      OY.var    = "OY",
      Y2.var    = "Y2",
      formula2  = ~ x1
    )
  )
})

test_that("impute_censored_stage2 returns list with expected components", {
  dat <- make_imputation_data(n = 60, seed = 3)
  # ensure at least one censored and several donors
  dat$delta <- c(rep(0L, 10), rep(1L, 50))
  dat$OY    <- sort(dat$OY)
  res <- impute_censored_stage2(
    dat,
    id.var    = "id",
    delta.var = "delta",
    OY.var    = "OY",
    Y2.var    = "Y2",
    formula2  = ~ x1 + x2,
    method    = "nearest",
    distance  = "mahalanobis",
    k         = 2,
    replace   = TRUE
  )
  expect_type(res, "list")
  expect_true("data_merged" %in% names(res))
  expect_true("n_censored"  %in% names(res))
  expect_true("n_imputed"   %in% names(res))
  expect_true("compY2"      %in% names(res$data_merged))
})

test_that("impute_censored_stage2 compY2 equals Y2 for uncensored subjects", {
  dat <- make_imputation_data(n = 60, seed = 4)
  dat$delta <- c(rep(0L, 10), rep(1L, 50))
  dat$OY    <- sort(dat$OY)
  res <- impute_censored_stage2(
    dat,
    id.var    = "id",
    delta.var = "delta",
    OY.var    = "OY",
    Y2.var    = "Y2",
    formula2  = ~ x1 + x2,
    method    = "nearest",
    replace   = TRUE
  )
  uncens <- res$data_merged[res$data_merged$delta == 1, ]
  expect_equal(uncens$compY2, uncens$Y2)
})

test_that("impute_censored_stage2 n_censored matches number of censored", {
  dat <- make_imputation_data(n = 60, seed = 5)
  dat$delta <- c(rep(0L, 8), rep(1L, 52))
  dat$OY    <- sort(dat$OY)
  res <- impute_censored_stage2(
    dat,
    id.var    = "id",
    delta.var = "delta",
    OY.var    = "OY",
    Y2.var    = "Y2",
    formula2  = ~ x1 + x2,
    method    = "nearest",
    replace   = TRUE
  )
  expect_equal(res$n_censored, 8L)
})

test_that("impute_censored_stage2 with aggregate='weighted' runs without error", {
  dat <- make_imputation_data(n = 60, seed = 6)
  dat$delta <- c(rep(0L, 10), rep(1L, 50))
  dat$OY    <- sort(dat$OY)
  res <- impute_censored_stage2(
    dat,
    id.var    = "id",
    delta.var = "delta",
    OY.var    = "OY",
    Y2.var    = "Y2",
    formula2  = ~ x1 + x2,
    aggregate = "weighted",
    replace   = TRUE
  )
  expect_type(res, "list")
})

test_that("impute_censored_stage2 with no censored subjects returns n_imputed=0", {
  dat <- make_imputation_data(n = 40, seed = 7)
  dat$delta <- 1L  # all observed
  res <- impute_censored_stage2(
    dat,
    id.var    = "id",
    delta.var = "delta",
    OY.var    = "OY",
    Y2.var    = "Y2",
    formula2  = ~ x1 + x2,
    method    = "nearest",
    replace   = TRUE
  )
  expect_equal(res$n_imputed, 0L)
})

# ---------------------------------------------------------------------------
# impute_censored_stage1 tests
# ---------------------------------------------------------------------------

test_that("impute_censored_stage1 stops when required columns are missing", {
  dat <- make_stage1_data()
  expect_error(
    impute_censored_stage1(
      dat,
      Id     = "id",
      formula = ~ x1 + x2,
      death1 = "MISSING",
      OY     = "OY",
      y1     = "Y1"
    ),
    "Missing required columns"
  )
})

test_that("impute_censored_stage1 stops when data is not a data.frame", {
  expect_error(
    impute_censored_stage1(
      as.matrix(make_stage1_data()),
      Id     = "id",
      formula = ~ x1 + x2,
      death1 = "death1",
      OY     = "OY",
      y1     = "Y1"
    )
  )
})

test_that("impute_censored_stage1 stops when y_cols is empty", {
  dat <- make_stage1_data()
  expect_error(
    impute_censored_stage1(
      dat,
      Id      = "id",
      formula = ~ x1 + x2,
      death1  = "death1",
      OY      = "OY",
      y1      = "Y1",
      y_cols  = character(0)
    ),
    "y_cols must be a non-empty character vector"
  )
})

test_that("impute_censored_stage1 returns list with expected components", {
  dat <- make_stage1_data(n = 60, seed = 8)
  dat$death1 <- c(rep(0L, 10), rep(1L, 50))
  dat$OY     <- sort(dat$OY)
  dat$compOY <- abs(rnorm(60, 12, 4))
  res <- impute_censored_stage1(
    dat,
    Id      = "id",
    formula = ~ x1 + x2,
    death1  = "death1",
    OY      = "OY",
    y1      = "Y1",
    y_cols  = "compOY",
    method  = "nearest",
    replace = TRUE
  )
  expect_type(res, "list")
  expect_true("data_merged" %in% names(res))
  expect_true("n_censored"  %in% names(res))
  expect_true("n_imputed"   %in% names(res))
})

test_that("impute_censored_stage1 n_censored matches number of censored subjects", {
  dat <- make_stage1_data(n = 60, seed = 9)
  dat$death1 <- c(rep(0L, 6), rep(1L, 54))
  dat$OY     <- sort(dat$OY)
  dat$compOY <- abs(rnorm(60, 12, 4))
  res <- impute_censored_stage1(
    dat,
    Id      = "id",
    formula = ~ x1 + x2,
    death1  = "death1",
    OY      = "OY",
    y1      = "Y1",
    y_cols  = "compOY",
    method  = "nearest",
    replace = TRUE
  )
  expect_equal(res$n_censored, 6L)
})

test_that("impute_censored_stage1 with no censored subjects returns n_imputed=0", {
  dat <- make_stage1_data(n = 40, seed = 10)
  dat$death1 <- 1L
  res <- impute_censored_stage1(
    dat,
    Id      = "id",
    formula = ~ x1 + x2,
    death1  = "death1",
    OY      = "OY",
    y1      = "Y1",
    y_cols  = "compOY",
    method  = "nearest",
    replace = TRUE
  )
  expect_equal(res$n_imputed, 0L)
})

# ---------------------------------------------------------------------------
# impute_censored_outcomes validation tests
# ---------------------------------------------------------------------------

test_that("impute_censored_outcomes stops when data is not a data.frame", {
  expect_error(
    impute_censored_outcomes(
      data      = list(a = 1),
      id.var    = "id",
      eta2.var  = "eta2",
      Y1.var    = "Y1",
      Y2.var    = "Y2",
      delta.var = "delta",
      OY.var    = "OY",
      A1.var    = "A1",
      A2.var    = "A2",
      names.var1 = "x1",
      names.var2 = "x1",
      tau        = 24
    )
  )
})

test_that("impute_censored_outcomes stops when required columns are missing", {
  dat <- data.frame(id = 1:10, eta2 = 1, Y1 = 1, Y2 = 2,
                    delta = 1, OY = 5, A1 = 1, A2 = 1, x1 = rnorm(10))
  expect_error(
    impute_censored_outcomes(
      data       = dat,
      id.var     = "id",
      eta2.var   = "eta2",
      Y1.var     = "Y1",
      Y2.var     = "Y2",
      delta.var  = "MISSING",
      OY.var     = "OY",
      A1.var     = "A1",
      A2.var     = "A2",
      names.var1 = "x1",
      names.var2 = "x1",
      tau        = 24,
      A.SL.library1 = "SL.glm",
      A.SL.library2 = "SL.glm",
      Y.SL.library  = "LIB_COXall"
    ),
    "Missing required columns"
  )
})

test_that("impute_censored_outcomes stops when names.var1 is empty", {
  dat <- data.frame(id = 1:10, eta2 = 1, Y1 = 1, Y2 = 2,
                    delta = 1, OY = 5, A1 = 1, A2 = 1, x1 = rnorm(10))
  expect_error(
    impute_censored_outcomes(
      data       = dat,
      id.var     = "id",
      eta2.var   = "eta2",
      Y1.var     = "Y1",
      Y2.var     = "Y2",
      delta.var  = "delta",
      OY.var     = "OY",
      A1.var     = "A1",
      A2.var     = "A2",
      names.var1 = character(0),
      names.var2 = "x1",
      tau        = 24,
      A.SL.library1 = "SL.glm",
      A.SL.library2 = "SL.glm",
      Y.SL.library  = "LIB_COXall"
    ),
    "names.var1 must be a non-empty character vector"
  )
})

test_that("impute_censored_outcomes stops when names.var2 columns not found", {
  dat <- data.frame(id = 1:10, eta2 = 1, Y1 = 1, Y2 = 2,
                    delta = 1, OY = 5, A1 = 1, A2 = 1, x1 = rnorm(10))
  expect_error(
    impute_censored_outcomes(
      data       = dat,
      id.var     = "id",
      eta2.var   = "eta2",
      Y1.var     = "Y1",
      Y2.var     = "Y2",
      delta.var  = "delta",
      OY.var     = "OY",
      A1.var     = "A1",
      A2.var     = "A2",
      names.var1 = "x1",
      names.var2 = "MISSING_VAR",
      tau        = 24,
      A.SL.library1 = "SL.glm",
      A.SL.library2 = "SL.glm",
      Y.SL.library  = "LIB_COXall"
    ),
    "names.var2 columns not found"
  )
})
