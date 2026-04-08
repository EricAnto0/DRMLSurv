#internal helper functions
#' Internal helper to relabel errors with a step name
#' @keywords internal
#' @noRd
.tag_error <- function(expr, label) {
  tryCatch(
    expr,
    error = function(e) {
      stop(sprintf("[%s] %s", label, conditionMessage(e)), call. = FALSE)
    }
  )
}

#' Internal helper to capture step failures
#' @keywords internal
#' @noRd
.capture_step <- function(expr, step, save_dir = "debug_logs", context = list()) {
  tryCatch(
    expr,
    error = function(e) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      file  <- file.path(save_dir, paste0(step, "_", stamp, ".rds"))
      
      info <- list(
        step = step,
        message = conditionMessage(e),
        class = class(e),
        call = conditionCall(e),
        sys.calls = vapply(
          sys.calls(),
          function(x) paste(deparse(x), collapse = " "),
          character(1)
        ),
        context = context,
        time = Sys.time()
      )
      
      saveRDS(info, file)
      
      message(sprintf("[%s] %s", step, conditionMessage(e)))
      message(sprintf("[%s] debug saved to %s", step, file))
      
      NULL
    }
  )
}
