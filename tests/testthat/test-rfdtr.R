make_rfdtr_data <- function(n = 80, p = 3, seed = 1) {
  set.seed(seed)
  covs <- as.data.frame(matrix(rnorm(n * p), nrow = n,
                               dimnames = list(NULL, paste0("x", seq_len(p)))))
  covs$A <- factor(sample(c(-1L, 1L), n, replace = TRUE), levels = c(-1, 1))
  covs
}

make_gridpar <- function() {
  data.frame(ntree = 50L, mtry = 2L, nodesize = 5L)
}

test_that("rfdtr stops when obs is not a data.frame", {
  gp <- make_gridpar()
  n  <- 30
  expect_error(
    rfdtr(
      obs      = as.matrix(make_rfdtr_data(n)),
      W        = rep(1, n),
      gridpar  = gp,
      A.obs    = rep(c(-1, 1), length.out = n),
      Q.obs    = rnorm(n),
      Q.match  = rnorm(n)
    ),
    "'obs' must be a data.frame"
  )
})

test_that("rfdtr stops when 'A' column is absent", {
  gp  <- make_gridpar()
  obs <- data.frame(x1 = rnorm(30))
  expect_error(
    rfdtr(
      obs     = obs,
      W       = rep(1, 30),
      gridpar = gp,
      A.obs   = rep(c(-1, 1), 15),
      Q.obs   = rnorm(30),
      Q.match = rnorm(30)
    ),
    "column named 'A'"
  )
})

test_that("rfdtr stops when W length != nrow(obs)", {
  gp  <- make_gridpar()
  obs <- make_rfdtr_data(30)
  expect_error(
    rfdtr(
      obs     = obs,
      W       = rep(1, 10),   # wrong length
      gridpar = gp,
      A.obs   = rep(c(-1, 1), 15),
      Q.obs   = rnorm(30),
      Q.match = rnorm(30)
    ),
    "Length of W"
  )
})

test_that("rfdtr stops when A.obs length != nrow(obs)", {
  gp  <- make_gridpar()
  obs <- make_rfdtr_data(30)
  expect_error(
    rfdtr(
      obs     = obs,
      W       = rep(1, 30),
      gridpar = gp,
      A.obs   = c(-1, 1),   # wrong length
      Q.obs   = rnorm(30),
      Q.match = rnorm(30)
    ),
    "must have length nrow"
  )
})

test_that("rfdtr stops when gridpar is not a data.frame", {
  obs <- make_rfdtr_data(30)
  expect_error(
    rfdtr(
      obs     = obs,
      W       = rep(1, 30),
      gridpar = list(ntree = 50, mtry = 2, nodesize = 5),
      A.obs   = rep(c(-1, 1), 15),
      Q.obs   = rnorm(30),
      Q.match = rnorm(30)
    ),
    "non-empty data.frame"
  )
})

test_that("rfdtr stops when gridpar is an empty data.frame", {
  obs <- make_rfdtr_data(30)
  expect_error(
    rfdtr(
      obs     = obs,
      W       = rep(1, 30),
      gridpar = data.frame(ntree = integer(0), mtry = integer(0), nodesize = integer(0)),
      A.obs   = rep(c(-1, 1), 15),
      Q.obs   = rnorm(30),
      Q.match = rnorm(30)
    ),
    "non-empty data.frame"
  )
})

test_that("rfdtr with usecv=FALSE and ranger returns model and predictions", {
  foreach::registerDoSEQ()

  obs <- make_rfdtr_data(n = 60, seed = 42)
  n   <- nrow(obs)
  gp  <- data.frame(ntree = 50L, mtry = 2L, nodesize = 5L)
  A.obs   <- as.numeric(as.character(obs$A))
  Q.obs   <- rnorm(n)
  Q.match <- rnorm(n)

  res <- rfdtr(
    modeltype = "ranger",
    usecv     = FALSE,
    sl.seed   = 1L,
    obs       = obs,
    W         = rep(1, n),
    gridpar   = gp,
    metric    = "ccr",
    A.obs     = A.obs,
    Q.obs     = Q.obs,
    Q.match   = Q.match
  )

  expect_type(res, "list")
  expect_true("model"    %in% names(res))
  expect_true("estA.obs" %in% names(res))
  expect_true("tune"     %in% names(res))
  expect_true("best"     %in% names(res))
  expect_length(res$estA.obs, n)
  expect_true(all(res$estA.obs %in% c(-1, 1)))
})

test_that("rfdtr metric='oob' selects by minimum OOB", {
  foreach::registerDoSEQ()

  obs <- make_rfdtr_data(n = 60, seed = 55)
  n   <- nrow(obs)
  gp  <- data.frame(ntree = c(50L, 100L), mtry = c(2L, 2L), nodesize = c(5L, 5L))
  A.obs   <- as.numeric(as.character(obs$A))
  Q.obs   <- rnorm(n)
  Q.match <- rnorm(n)

  res <- rfdtr(
    modeltype = "ranger",
    usecv     = FALSE,
    sl.seed   = 1L,
    obs       = obs,
    W         = rep(1, n),
    gridpar   = gp,
    metric    = "oob",
    A.obs     = A.obs,
    Q.obs     = Q.obs,
    Q.match   = Q.match
  )
  expect_equal(res$best$OOB, min(res$tune$OOB))
})

test_that("rfdtr metric='policyval' selects by maximum Score", {
  foreach::registerDoSEQ()

  obs <- make_rfdtr_data(n = 60, seed = 66)
  n   <- nrow(obs)
  gp  <- data.frame(ntree = c(50L, 100L), mtry = c(2L, 2L), nodesize = c(5L, 5L))
  A.obs   <- as.numeric(as.character(obs$A))
  Q.obs   <- rnorm(n)
  Q.match <- rnorm(n)

  res <- rfdtr(
    modeltype = "ranger",
    usecv     = FALSE,
    sl.seed   = 1L,
    obs       = obs,
    W         = rep(1, n),
    gridpar   = gp,
    metric    = "policyval",
    A.obs     = A.obs,
    Q.obs     = Q.obs,
    Q.match   = Q.match
  )
  expect_equal(res$best$Score, max(res$tune$Score))
})

test_that("rfdtr mtry is clamped to [1, p] automatically", {
  foreach::registerDoSEQ()

  obs <- make_rfdtr_data(n = 60, p = 3, seed = 77)
  n   <- nrow(obs)
  gp  <- data.frame(ntree = 50L, mtry = 999L, nodesize = 5L)
  A.obs   <- as.numeric(as.character(obs$A))
  Q.obs   <- rnorm(n)
  Q.match <- rnorm(n)

  expect_no_error(
    rfdtr(
      modeltype = "ranger",
      usecv     = FALSE,
      sl.seed   = 1L,
      obs       = obs,
      W         = rep(1, n),
      gridpar   = gp,
      metric    = "ccr",
      A.obs     = A.obs,
      Q.obs     = Q.obs,
      Q.match   = Q.match
    )
  )
})
