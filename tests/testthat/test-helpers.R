test_that(".tag_error re-labels error messages with step name", {
  f <- DRMLSurv:::.tag_error
  expect_error(
    f(stop("original error"), "MyStep"),
    "\\[MyStep\\] original error"
  )
})

test_that(".tag_error passes through successful expressions", {
  f <- DRMLSurv:::.tag_error
  result <- f(1 + 1, "MyStep")
  expect_equal(result, 2)
})

test_that(".capture_step returns result on success", {
  f <- DRMLSurv:::.capture_step
  result <- f(42L, step = "test_step")
  expect_equal(result, 42L)
})

test_that(".capture_step returns NULL and shows message on error", {
  f <- DRMLSurv:::.capture_step
  expect_message(
    result <- f(stop("something went wrong"), step = "fail_step"),
    "something went wrong"
  )
  expect_null(result)
})

test_that(".capture_step saves context to debug log on error", {
  f <- DRMLSurv:::.capture_step
  tmp_dir <- tempfile("debug_test_")
  suppressMessages(
    f(stop("logged error"), step = "ctx_step",
      save_dir = tmp_dir, context = list(foo = "bar"))
  )
  saved_files <- list.files(tmp_dir, full.names = TRUE)
  expect_true(length(saved_files) >= 1L)
  info <- readRDS(saved_files[[1]])
  expect_equal(info$step, "ctx_step")
  expect_equal(info$context$foo, "bar")
  unlink(tmp_dir, recursive = TRUE)
})
