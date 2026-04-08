make_policy_data <- function(n = 60) {
  set.seed(42)
  data.frame(
    A1.var = sample(c(-1L, 1L), n, replace = TRUE),
    A2.var = sample(c(-1L, 1L), n, replace = TRUE),
    eta2   = sample(c(0L, 1L),  n, replace = TRUE, prob = c(0.3, 0.7)),
    val_c1 = runif(n, 5, 10),
    val_t1 = runif(n, 8, 15),
    val_c2 = runif(n, 5, 10),
    val_t2 = runif(n, 8, 15),
    val_c_tot = runif(n, 5, 10),
    val_t_tot = runif(n, 8, 15)
  )
}

test_that("policy_summary_metrics returns a list for stage='both'", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2")
  expect_type(res, "list")
  expect_true("Acc1L" %in% names(res))
  expect_true("Acc2L" %in% names(res))
  expect_true("AccTotal" %in% names(res))
})

test_that("policy_summary_metrics stage='stage1' only populates stage-1 metrics", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2", stage = "stage1")
  expect_true(is.numeric(res$Acc1L))
  expect_true(is.na(res$Acc2L))
  expect_true(is.na(res$AccTotal))
})

test_that("policy_summary_metrics stage='stage2' only populates stage-2 metrics", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2", stage = "stage2")
  expect_true(is.na(res$Acc1L))
  expect_true(is.numeric(res$Acc2L))
  expect_true(is.na(res$AccTotal))
})

test_that("policy_summary_metrics computes RMST1L when value cols supplied", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2",
                                   Tc1_var  = "val_c1", Tt1_var = "val_t1",
                                   stage    = "stage1")
  expect_true(is.numeric(res$RMST1L))
  expect_false(is.na(res$RMST1L))
})

test_that("policy_summary_metrics computes RMST2L for stage2 with value cols", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2",
                                   Tc2_var  = "val_c2", Tt2_var = "val_t2",
                                   stage    = "stage2")
  expect_true(is.numeric(res$RMST2L))
  expect_false(is.na(res$RMST2L))
})

test_that("policy_summary_metrics computes RMSTTotal with total value cols", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var     = "A1.var", obs2_var = "A2.var",
                                   eta2_var     = "eta2",
                                   Tc_total_var = "val_c_tot",
                                   Tt_total_var = "val_t_tot")
  expect_true(is.numeric(res$RMSTTotal))
  expect_false(is.na(res$RMSTTotal))
})

test_that("policy_summary_metrics accuracy is in [0,1]", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2")
  expect_gte(res$Acc1L, 0)
  expect_lte(res$Acc1L, 1)
})

test_that("policy_summary_metrics: perfect predictions give Acc1L=1", {
  dat   <- make_policy_data(40)
  estA1 <- dat$A1.var
  estA2 <- dat$A2.var
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2")
  expect_equal(res$Acc1L, 1)
})

test_that("policy_summary_metrics stops when tmpData not a data.frame", {
  expect_error(
    policy_summary_metrics(list(A1.var = 1), c(1), c(1)),
    "must be a data.frame"
  )
})

test_that("policy_summary_metrics stops when obs1_var missing", {
  dat <- make_policy_data(10)
  expect_error(
    policy_summary_metrics(dat, rep(1, 10), rep(1, 10),
                            obs1_var = "MISSING_COL",
                            obs2_var = "A2.var", eta2_var = "eta2"),
    "not found"
  )
})

test_that("policy_summary_metrics stops when obs2_var missing", {
  dat <- make_policy_data(10)
  expect_error(
    policy_summary_metrics(dat, rep(1, 10), rep(1, 10),
                            obs1_var = "A1.var",
                            obs2_var = "MISSING_COL", eta2_var = "eta2"),
    "not found"
  )
})

test_that("policy_summary_metrics stops when eta2_var missing", {
  dat <- make_policy_data(10)
  expect_error(
    policy_summary_metrics(dat, rep(1, 10), rep(1, 10),
                            obs1_var = "A1.var",
                            obs2_var = "A2.var", eta2_var = "MISSING_COL"),
    "not found"
  )
})

test_that("policy_summary_metrics stops when estA1 wrong length", {
  dat <- make_policy_data(10)
  expect_error(
    policy_summary_metrics(dat, rep(1, 5), rep(1, 10),
                            obs1_var = "A1.var", obs2_var = "A2.var",
                            eta2_var = "eta2"),
    "length nrow"
  )
})

test_that("policy_summary_metrics stops when estA2 wrong length", {
  dat <- make_policy_data(10)
  expect_error(
    policy_summary_metrics(dat, rep(1, 10), rep(1, 5),
                            obs1_var = "A1.var", obs2_var = "A2.var",
                            eta2_var = "eta2"),
    "length nrow"
  )
})

test_that("policy_summary_metrics handles no stage-2 subjects gracefully", {
  dat   <- make_policy_data(20)
  dat$eta2 <- 0L
  estA1 <- rep(1L, nrow(dat))
  estA2 <- rep(1L, nrow(dat))
  res <- policy_summary_metrics(dat, estA1, estA2,
                                 obs1_var = "A1.var", obs2_var = "A2.var",
                                 eta2_var = "eta2")
  expect_true(is.na(res$Acc2L))
  expect_true(is.na(res$AccTotal))
})

test_that("policy_summary_metrics MCC is in [-1,1] or NA", {
  dat   <- make_policy_data()
  estA1 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  estA2 <- sample(c(-1L, 1L), nrow(dat), replace = TRUE)
  res   <- policy_summary_metrics(dat, estA1, estA2,
                                   obs1_var = "A1.var", obs2_var = "A2.var",
                                   eta2_var = "eta2")
  if (!is.na(res$MCC1L)) {
    expect_gte(res$MCC1L, -1)
    expect_lte(res$MCC1L, 1)
  }
})
