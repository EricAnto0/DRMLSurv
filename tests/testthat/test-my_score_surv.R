test_that("my_score.Surv returns correct mean when all predictions match", {
  pred <- c(1, 1, -1, -1)
  A    <- c(1, 1, -1, -1)
  Q    <- c(10, 20, 30, 40)
  pQ   <- c(5, 15, 25, 35)
  # all pred == A, so result = mean(Q)
  expect_equal(my_score.Surv(pred, A, Q, pQ), mean(Q))
})

test_that("my_score.Surv returns correct mean when no predictions match", {
  pred <- c(1, 1, -1, -1)
  A    <- c(-1, -1, 1, 1)
  Q    <- c(10, 20, 30, 40)
  pQ   <- c(5, 15, 25, 35)
  # all pred != A, so result = mean(pQ)
  expect_equal(my_score.Surv(pred, A, Q, pQ), mean(pQ))
})

test_that("my_score.Surv handles mixed concordance correctly", {
  pred <- c(1, -1, 1, -1)
  A    <- c(1,  1, 1, -1)
  Q    <- c(10, 12,  8, 15)
  pQ   <- c( 9, 11,  7, 14)
  # pred==A at pos 1,3,4 -> Q[1]=10, Q[3]=8, Q[4]=15
  # pred!=A at pos 2    -> pQ[2]=11
  expect_equal(my_score.Surv(pred, A, Q, pQ), mean(c(10, 11, 8, 15)))
})

test_that("my_score.Surv stops on non-numeric Q", {
  expect_error(
    my_score.Surv(c(1, -1), c(1, -1), c("a", "b"), c(1, 2)),
    "'Q' and 'pQ' must be numeric"
  )
})

test_that("my_score.Surv stops on non-numeric pQ", {
  expect_error(
    my_score.Surv(c(1, -1), c(1, -1), c(1, 2), c("a", "b")),
    "'Q' and 'pQ' must be numeric"
  )
})

test_that("my_score.Surv stops when lengths differ", {
  expect_error(
    my_score.Surv(c(1, -1), c(1, -1), c(1, 2, 3), c(1, 2, 3)),
    "must all have the same length"
  )
})

test_that("my_score.Surv handles NA values with na.rm", {
  pred <- c(1, -1, 1)
  A    <- c(1,  1, 1)
  Q    <- c(10, NA, 8)
  pQ   <- c(9, 11, 7)
  # pos 1: pred==A -> Q[1]=10; pos 2: pred!=A -> pQ[2]=11 (NA skipped); pos 3: pred==A -> Q[3]=8
  result <- my_score.Surv(pred, A, Q, pQ)
  expect_equal(result, mean(c(10, 11, 8), na.rm = TRUE))
})

test_that("my_score.Surv works with length-1 inputs", {
  expect_equal(my_score.Surv(1, 1, 5, 9), 5)
  expect_equal(my_score.Surv(1, -1, 5, 9), 9)
})
