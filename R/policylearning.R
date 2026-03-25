#' Tune and fit a random-forest DTR policy model with (optional) cross-validation
#'
#' @description
#' Fits a classification-based dynamic treatment regime (DTR) policy model for a binary treatment
#' learns a binary decision rule \eqn{\hat d(X)} with values in \eqn{\{-1,+1\}} using either \pkg{ranger} or \pkg{randomForestSRC}. The function
#' performs grid search over random-forest tuning parameters and, when \code{usecv=TRUE}, evaluates
#' candidate policies via K-fold cross-validation using:
#' \itemize{
#'   \item \strong{CCR}: classification correctness rate (accuracy of predicted treatment labels),
#'   \item \strong{OOB}: out-of-bag prediction error (model-reported), and
#'   \item \strong{Score}: a user-supplied policy-value-like performance measure computed by
#'   \code{my_score.Surv()} using observed and matched pseudo-outcomes.
#' }
#'
#' The best tuning parameters are selected according to \code{metric}, after which a final model is fit
#' on the full observed dataset \code{obs}. The output includes the final fitted model, the estimated
#' treatment rule \code{estA.obs}, and the full tuning results.
#'
#' @details
#' \strong{Inputs and data layout.}
#' \itemize{
#'   \item \code{obs} is a data.frame that must contain a column named \code{A} and the covariates used to
#'   model \code{A}. The function coerces \code{obs$A} to a factor with levels \code{c(-1, 1)} to ensure
#'   a binary classification target.
#'   \item \code{W} is a numeric vector of case weights aligned to \code{obs} (e.g., IPC weights).
#'   \item \code{A.obs} is the numeric observed treatment indicator (typically \code{-1/1}) aligned to \code{obs}.
#'   \item \code{Q.obs} and \code{Q.match} are numeric vectors aligned to \code{obs} used by \code{my_score.Surv()}
#'   to compute the policy score (e.g., observed pseudo-outcome and matched/pair pseudo-outcome).
#' }
#'
#' \strong{Cross-validation logic.}
#' When \code{usecv=TRUE}, the function uses 5-fold stratified folds created by \code{caret::createFolds(obs$A, k=5)}.
#' For each candidate parameter set in \code{gridpar}, it fits the model on the training folds and evaluates on the
#' held-out fold:
#' \itemize{
#'   \item \emph{CCR}: mean(\code{predicted_class == obs$A[test]}).
#'   \item \emph{OOB}: model-reported OOB error (for \pkg{ranger}, \code{mod$prediction.error}; for \pkg{rfsrc},
#'   the last entry of \code{mod$err.rate[,"all"]}).
#'   \item \emph{Score}: \code{my_score.Surv(pred, A.test, Q.test, Q.match.test)} where \code{pred} is the numeric
#'   treatment rule \code{-1/1}.
#' }
#'
#' \strong{Aggregation of score across folds.}
#' The per-fold \code{Score} values are aggregated across CV folds using:
#' \itemize{
#'   \item \code{score_agg="sum"}: sum of fold scores (ignoring \code{NA})
#'   \item \code{score_agg="mean"}: mean of fold scores (ignoring \code{NA})
#' }
#' This controls the scale used for model selection when \code{metric} targets policy value.
#'
#' \strong{Parameter guards.}
#' To avoid invalid random-forest hyperparameters, if \code{gridpar} includes \code{mtry}, it is clamped
#' to \code{[1, p]} where \code{p = ncol(obs)-1} (i.e., all columns except \code{A}). Duplicate parameter
#' rows are removed via \code{unique(gridpar)}.
#'
#' \strong{Parallelization.}
#' The function is written using \pkg{foreach} with \code{%dopar%}. It assumes a parallel backend has
#' already been registered (e.g., via \pkg{doParallel}). Within each fold, \code{num.threads=1} is used
#' for \pkg{ranger} to avoid nested parallelism.
#'
#' @param modeltype Character. Random-forest engine: \code{"ranger"} or \code{"rfsrc"}.
#'
#' @param usecv Logical. If TRUE, uses 5-fold CV to tune hyperparameters; otherwise tunes on the full data.
#'
#' @param sl.seed Integer. Seed for reproducibility. (Note: parts of the function currently hard-code \code{set.seed(123)}.)
#'
#' @param obs data.frame. Training dataset containing a column \code{A} and covariates for predicting \code{A}.
#'
#' @param W Numeric vector. Case weights of length \code{nrow(obs)}.
#'
#' @param gridpar data.frame. Grid of tuning parameters. Expected columns include:
#' \code{ntree}, \code{mtry}, and \code{nodesize}. Extra columns are ignored by the model fits.
#'
#' @param metric Character. Criterion used to choose the “best” tuning row. Supported:
#' \itemize{
#'   \item \code{"oob"}: minimize \code{OOB}
#'   \item \code{"policyval"} / \code{"score"} / \code{"policy"} / \code{"val"}: maximize \code{Score}
#'   \item otherwise: maximize \code{CCR}
#' }
#'
#' @param A.obs Numeric vector. Observed treatment labels (\code{-1/1}) aligned with \code{obs}.
#'
#' @param Q.obs Numeric vector. Observed pseudo-outcome or value component aligned with \code{obs}.
#'
#' @param Q.match Numeric vector. Matched/pair pseudo-outcome aligned with \code{obs} used by \code{my_score.Surv()}.
#'
#' @param score_agg Character. Aggregation for fold-level \code{Score}: \code{"sum"} or \code{"mean"}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{model}: final fitted \pkg{ranger} or \pkg{randomForestSRC} object trained on all \code{obs}.
#'   \item \code{estA.obs}: numeric vector of predicted treatment decisions (\code{-1/1}) for each row of \code{obs}.
#'   \item \code{tune}: data.frame of tuning results for each row of \code{gridpar}, including \code{CCR}, \code{OOB}, \code{Score}.
#'   \item \code{best}: the selected row of \code{tune} corresponding to the chosen \code{metric}.
#' }
#'
#' @seealso \code{\link[ranger]{ranger}}, \code{\link[randomForestSRC]{rfsrc}}, \code{\link[caret]{createFolds}}
#' @export
rfdtr <- function(
    modeltype = "ranger", usecv = TRUE, sl.seed = 123,
    obs, W, gridpar, metric = "ccr", A.obs, Q.obs, Q.match,
    score_agg = c("mean", "sum")
) {
  score_agg <- match.arg(score_agg)
  agg_fun <- if (score_agg == "sum") function(x) sum(x, na.rm = TRUE) else function(x) mean(x, na.rm = TRUE)

  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' is required (register a backend for %dopar%).", call. = FALSE)
  }
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required for createFolds().", call. = FALSE)
  }
  if (!exists("my_score.Surv", mode = "function")) {
    stop("Function my_score.Surv() must exist in the environment/package.", call. = FALSE)
  }
  if (!is.data.frame(obs)) stop("'obs' must be a data.frame.", call. = FALSE)
  if (!"A" %in% names(obs)) stop("'obs' must contain a column named 'A'.", call. = FALSE)
  if (length(W) != nrow(obs)) stop("Length of W must equal nrow(obs).", call. = FALSE)
  if (length(A.obs) != nrow(obs) || length(Q.obs) != nrow(obs) || length(Q.match) != nrow(obs)) {
    stop("A.obs, Q.obs, and Q.match must have length nrow(obs).", call. = FALSE)
  }
  if (!is.data.frame(gridpar) || nrow(gridpar) < 1L) stop("'gridpar' must be a non-empty data.frame.", call. = FALSE)

  # Basic guards & prep
  pp <- ncol(obs) - 1L
  if ("mtry" %in% names(gridpar)) {
    gridpar$mtry <- pmin(pmax(1L, gridpar$mtry), pp)
  }
  gridpar <- unique(gridpar)

  # ensure factor for classification; preserve {-1,1} levels
  obs$A <- factor(obs$A, levels = c(-1, 1))

  metric <- tolower(metric)

  if (usecv) {
    # 5-fold CV (stratified by A)
    set.seed(sl.seed)
    fold_ids <- caret::createFolds(obs$A, k = 5)
    n_vec    <- vapply(fold_ids, length, integer(1))

    if (modeltype == "ranger") {
      if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' is required.", call. = FALSE)

      tunedf <- foreach::foreach(
        g = seq_len(nrow(gridpar)), .combine = rbind,
        .packages = "ranger", .export = "my_score.Surv"
      ) %dopar% {
        pars <- gridpar[g, , drop = TRUE]

        fold_perf <- sapply(fold_ids, function(test_idx) {
          train_idx <- setdiff(seq_len(nrow(obs)), test_idx)

          mod <- ranger::ranger(
            A ~ .,
            data = obs[train_idx, , drop = FALSE],
            num.trees = pars$ntree,
            mtry = pars$mtry,
            min.node.size = pars$nodesize,
            case.weights = W[train_idx],
            classification = TRUE,
            importance = "none",
            num.threads = 1, # avoid nested parallelism
            respect.unordered.factors = "order",
            seed = sl.seed
          )

          R.test       <- Q.obs[test_idx]
          R.match.test <- Q.match[test_idx]
          A.test       <- A.obs[test_idx]

          pred_fac <- predict(mod, data = obs[test_idx, , drop = FALSE])$predictions
          pred_num <- as.numeric(as.character(pred_fac))

          score <- my_score.Surv(pred_num, as.numeric(as.character(A.test)), R.test, R.match.test)

          c(
            mean(pred_fac == obs$A[test_idx]),
            mod$prediction.error,
            score
          )
        })

        data.frame(
          ntree    = pars$ntree,
          mtry     = pars$mtry,
          nodesize = pars$nodesize,
          CCR      = mean(fold_perf[1, ], na.rm = TRUE),
          OOB      = mean(fold_perf[2, ], na.rm = TRUE),
          Score    = agg_fun(fold_perf[3, ]),
          stringsAsFactors = FALSE
        )
      }

    } else if (modeltype == "rfsrc") {
      if (!requireNamespace("randomForestSRC", quietly = TRUE)) stop("Package 'randomForestSRC' is required.", call. = FALSE)

      tunedf <- foreach::foreach(
        g = seq_len(nrow(gridpar)), .combine = rbind,
        .packages = "randomForestSRC", .export = "my_score.Surv"
      ) %dopar% {
        pars <- gridpar[g, , drop = TRUE]

        fold_perf <- sapply(fold_ids, function(test_idx) {
          train_idx <- setdiff(seq_len(nrow(obs)), test_idx)

          mod <- randomForestSRC::rfsrc(
            A ~ .,
            data = obs[train_idx, , drop = FALSE],
            case.wt = W[train_idx],
            ntree = pars$ntree,
            mtry = pars$mtry,
            nodesize = pars$nodesize,
            importance = TRUE,
            samptype = "swr",
            do.trace = FALSE,
            seed = sl.seed
          )

          R.test       <- Q.obs[test_idx]
          R.match.test <- Q.match[test_idx]
          A.test       <- A.obs[test_idx]

          pred_num <- as.numeric(as.character(predict(mod, obs[test_idx, , drop = FALSE])$class))
          score <- my_score.Surv(pred_num, as.numeric(as.character(A.test)), R.test, R.match.test)

          c(
            mean(pred_num == as.numeric(as.character(obs$A[test_idx]))),
            tail(mod$err.rate[, "all"], 1),
            score
          )
        })

        data.frame(
          ntree    = pars$ntree,
          mtry     = pars$mtry,
          nodesize = pars$nodesize,
          CCR      = mean(fold_perf[1, ], na.rm = TRUE),
          OOB      = mean(fold_perf[2, ], na.rm = TRUE),
          Score    = agg_fun(fold_perf[3, ]),
          stringsAsFactors = FALSE
        )
      }

    } else {
      stop("Unsupported modeltype: use 'ranger' or 'rfsrc'.", call. = FALSE)
    }

  } else {
    # No CV: evaluate on full data
    if (modeltype == "ranger" && !requireNamespace("ranger", quietly = TRUE)) {
      stop("Package 'ranger' is required.", call. = FALSE)
    }
    if (modeltype == "rfsrc" && !requireNamespace("randomForestSRC", quietly = TRUE)) {
      stop("Package 'randomForestSRC' is required.", call. = FALSE)
    }

    tunedf <- foreach::foreach(
      i = seq_len(nrow(gridpar)), .combine = rbind,
      .packages = c("ranger", "randomForestSRC"), .export = "my_score.Surv"
    ) %dopar% {
      pars <- gridpar[i, , drop = TRUE]

      if (modeltype == "ranger") {
        mod <- ranger::ranger(
          A ~ ., data = obs,
          num.trees = pars$ntree,
          mtry = pars$mtry,
          min.node.size = pars$nodesize,
          case.weights = W,
          classification = TRUE,
          importance = "none",
          num.threads = 1,
          respect.unordered.factors = "order",
          seed = sl.seed
        )
        pred_fac <- predict(mod, data = obs)$predictions
        pred_num <- as.numeric(as.character(pred_fac))
        ccr     <- mean(pred_fac == obs$A, na.rm = TRUE)
        oob_err <- mod$prediction.error
        score   <- my_score.Surv(pred_num, as.numeric(as.character(A.obs)), Q.obs, Q.match)

      } else {
        mod <- randomForestSRC::rfsrc(
          A ~ ., data = obs, case.wt = W,
          ntree = pars$ntree, mtry = pars$mtry, nodesize = pars$nodesize,
          importance = TRUE, samptype = "swr", do.trace = FALSE, seed = sl.seed
        )
        pred_num <- as.numeric(as.character(predict(mod, obs)$class))
        ccr     <- mean(pred_num == as.numeric(as.character(obs$A)), na.rm = TRUE)
        oob_err <- tail(mod$err.rate[, "all"], 1)
        score   <- my_score.Surv(pred_num, as.numeric(as.character(A.obs)), Q.obs, Q.match)
      }

      data.frame(
        ntree = pars$ntree,
        mtry = pars$mtry,
        nodesize = pars$nodesize,
        CCR = ccr,
        OOB = oob_err,
        Score = score,
        stringsAsFactors = FALSE
      )
    }
  }

  # Pick best by requested metric
  best <- if (metric == "oob") {
    tunedf[which.min(tunedf$OOB), , drop = FALSE]
  } else if (metric %in% c("policyval","score","policy","val")) {
    tunedf[which.max(tunedf$Score), , drop = FALSE]
  } else {
    tunedf[which.max(tunedf$CCR), , drop = FALSE]
  }

  cat("Best set of tuning parameters and metrics overall\n\n",
      paste0(names(best), " = ", unlist(best), collapse = "\n"), "\n")

  # Fit final model and predictions on obs
  if (modeltype == "ranger") {
    final <- ranger::ranger(
      A ~ ., data = obs,
      num.trees = best$ntree, mtry = best$mtry, min.node.size = best$nodesize,
      case.weights = W, classification = TRUE, importance = "none",
      respect.unordered.factors = "order", seed = sl.seed
    )
    estA.obs <- as.numeric(as.character(predict(final, data = obs)$predictions))

  } else if (modeltype == "rfsrc") {
    final <- randomForestSRC::rfsrc(
      A ~ ., data = obs, case.wt = W,
      ntree = best$ntree, mtry = best$mtry, nodesize = best$nodesize,
      importance = TRUE, samptype = "swr", do.trace = FALSE, seed = sl.seed
    )
    estA.obs <- as.numeric(as.character(predict(final, obs)$class))

  } else {
    stop("Unsupported modeltype.", call. = FALSE)
  }

  list(model = final, estA.obs = estA.obs, tune = tunedf, best = best)
}



#' Compute a value score under a candidate treatment rule
#'
#' Evaluates the mean outcome under a candidate treatment rule by assigning,
#' for each subject, the outcome corresponding to the treatment actually
#' recommended by the rule. If the predicted treatment matches the observed
#' treatment, the observed outcome under the received treatment is used;
#' otherwise, the alternative counterfactual or model-based outcome is used.
#'
#' This function is useful in dynamic treatment regime or policy-learning
#' settings where one wants to estimate the value of a learned treatment rule
#' using observed and predicted treatment assignments together with outcome
#' quantities under concordant and discordant treatment choices.
#'
#' @param pred A vector of predicted treatment assignments under the candidate
#'   rule.
#' @param A A vector of observed treatment assignments.
#' @param Q A numeric vector giving the outcome value to use when the predicted
#'   treatment matches the observed treatment, i.e., for subjects with
#'   `pred == A`.
#' @param pQ A numeric vector giving the outcome value to use when the predicted
#'   treatment does not match the observed treatment, i.e., for subjects with
#'   `pred != A`.
#'
#' @details
#' For each subject, the function constructs
#' \deqn{
#' Q_i^* =
#' \begin{cases}
#' Q_i, & \text{if } pred_i = A_i, \\
#' pQ_i, & \text{if } pred_i \neq A_i.
#' \end{cases}
#' }
#' and returns the sample mean of \eqn{Q_i^*}.
#'
#' The function assumes that `pred`, `A`, `Q`, and `pQ` are all aligned and
#' of equal length.
#'
#' @return
#' A single numeric value equal to the mean outcome under the candidate
#' treatment rule.
#'
#' @examples
#' pred <- c(1, -1, 1, -1)
#' A    <- c(1,  1, 1, -1)
#' Q    <- c(10, 12,  8, 15)
#' pQ   <- c( 9, 11,  7, 14)
#'
#' my_score.Surv(pred = pred, A = A, Q = Q, pQ = pQ)
#'
#' @export
my_score.Surv <- function(pred, A, Q, pQ) {
  if (!is.numeric(Q) || !is.numeric(pQ)) {
    stop("'Q' and 'pQ' must be numeric.")
  }
  if (!(length(pred) == length(A) &&
        length(A) == length(Q) &&
        length(Q) == length(pQ))) {
    stop("'pred', 'A', 'Q', and 'pQ' must all have the same length.")
  }

  Q0 <- ifelse(pred == A, Q, pQ)
  mean(Q0, na.rm = TRUE)
}
