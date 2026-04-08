make_score_data <- function(n = 40, seed = 1) {
  set.seed(seed)
  data.frame(
    id    = seq_len(n),
    Y     = abs(rnorm(n, 10, 3)),
    event = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.3, 0.7)),
    A     = sample(c(0L, 1L), n, replace = TRUE),
    x1    = rnorm(n),
    x2    = rnorm(n),
    stringsAsFactors = FALSE
  )
}

test_that("ComputeScores stops when data is not a data.frame", {
  expect_error(
    ComputeScores(
      data  = list(id = 1),
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A"
    )
  )
})

test_that("ComputeScores stops when required column missing", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "MISSING",
      event = "event",
      X     = "x1",
      A     = "A"
    ),
    "not found in data"
  )
})

test_that("ComputeScores stops when X is not a character vector", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = 123,
      A     = "A"
    ),
    "X must be a character vector"
  )
})

test_that("ComputeScores stops when X column not found", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "MISSING_COV",
      A     = "A"
    ),
    "Some X columns not found"
  )
})

test_that("ComputeScores stops when Xtrt is not a character vector", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A",
      Xtrt  = 999
    ),
    "Xtrt must be NULL or a character vector"
  )
})

test_that("ComputeScores stops when Xtrt column not found in data", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A",
      Xtrt  = "MISSING_TRT_COV"
    ),
    "Some Xtrt columns not found"
  )
})

test_that("ComputeScores stops when outer_CV < 2", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data     = dat,
      id       = "id",
      Y        = "Y",
      event    = "event",
      X        = "x1",
      A        = "A",
      outer_CV = 1L
    ),
    "outer_CV must be an integer >= 2"
  )
})

test_that("ComputeScores stops when inner_CV < 2", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data     = dat,
      id       = "id",
      Y        = "Y",
      event    = "event",
      X        = "x1",
      A        = "A",
      inner_CV = 1L
    ),
    "inner_CV must be NULL or an integer >= 2"
  )
})

test_that("ComputeScores stops when cores < 1", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A",
      cores = 0L
    ),
    "cores must be an integer >= 1"
  )
})

test_that("ComputeScores stops when ngrid < 2", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A",
      ngrid = 1L
    ),
    "ngrid must be an integer >= 2"
  )
})

test_that("ComputeScores stops when tau is non-positive", {
  dat <- make_score_data()
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A",
      tau   = -5
    ),
    "tau must be NULL or a positive numeric value"
  )
})

test_that("ComputeScores stops when event column not coded 0/1", {
  dat <- make_score_data()
  dat$event <- sample(c(0L, 2L), nrow(dat), replace = TRUE)
  expect_error(
    ComputeScores(
      data  = dat,
      id    = "id",
      Y     = "Y",
      event = "event",
      X     = "x1",
      A     = "A"
    ),
    "must be coded 0/1"
  )
})

# ---------------------------------------------------------------------------
# get_doublescores tests
# ---------------------------------------------------------------------------

make_ds_data <- function(n = 40, seed = 1) {
  set.seed(seed)
  data.frame(
    id    = seq_len(n),
    eta2  = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.4, 0.6)),
    Y1    = abs(rnorm(n, 5, 2)),
    Y2    = abs(rnorm(n, 6, 2)),
    delta = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.3, 0.7)),
    OY    = abs(rnorm(n, 12, 4)),
    A1    = sample(c(0L, 1L), n, replace = TRUE),
    A2    = sample(c(0L, 1L), n, replace = TRUE),
    x1    = rnorm(n),
    x2    = rnorm(n),
    stringsAsFactors = FALSE
  )
}

test_that("get_doublescores with useds=FALSE returns original data unchanged", {
  dat <- make_ds_data()
  res <- get_doublescores(
    data       = dat,
    id.var     = "id",
    eta2.var   = "eta2",
    Y1.var     = "Y1",
    Y2.var     = "Y2",
    delta.var  = "delta",
    OY.var     = "OY",
    A1.var     = "A1",
    A2.var     = "A2",
    names.var1 = c("x1", "x2"),
    names.var2 = c("x1", "x2"),
    useds      = FALSE,
    tau        = 24,
    A.SL.library1 = "SL.glm",
    A.SL.library2 = "SL.glm",
    Y.SL.library  = "LIB_COXall"
  )
  expect_identical(res, dat)
})

test_that("get_doublescores stops when required columns are missing", {
  dat <- make_ds_data()
  expect_error(
    get_doublescores(
      data       = dat,
      id.var     = "id",
      eta2.var   = "eta2",
      Y1.var     = "Y1",
      Y2.var     = "Y2",
      delta.var  = "MISSING_DELTA",
      OY.var     = "OY",
      A1.var     = "A1",
      A2.var     = "A2",
      names.var1 = c("x1", "x2"),
      names.var2 = c("x1", "x2"),
      useds      = TRUE,
      tau        = 24,
      A.SL.library1 = "SL.glm",
      A.SL.library2 = "SL.glm",
      Y.SL.library  = "LIB_COXall"
    ),
    "Missing required columns"
  )
})

test_that("get_doublescores stops when data is not a data.frame", {
  expect_error(
    get_doublescores(
      data       = list(a = 1),
      id.var     = "id",
      eta2.var   = "eta2",
      Y1.var     = "Y1",
      Y2.var     = "Y2",
      delta.var  = "delta",
      OY.var     = "OY",
      A1.var     = "A1",
      A2.var     = "A2",
      names.var1 = "x1",
      names.var2 = "x1",
      useds      = FALSE,
      tau        = 24,
      A.SL.library1 = "SL.glm",
      A.SL.library2 = "SL.glm",
      Y.SL.library  = "LIB_COXall"
    )
  )
})
