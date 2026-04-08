
#' Compute treatment, prognostic, and optional censoring-related scores
#'
#' @description
#' Computes subject-level score summaries for downstream matching, weighting, or
#' augmentation procedures. The function always estimates a treatment propensity
#' score and can additionally estimate treatment-specific prognostic scores and
#' censoring-related scores, depending on the supplied flags.
#'
#' Estimation may be carried out using SuperLearner / survivalSL-based models or
#' via glm / glmnet / Cox / AFT alternatives.
#'
#' @details
#' \strong{Treatment propensity score.}
#'
#' The function estimates \code{ps}, the probability of treatment conditional on
#' \code{Xtrt} when supplied, or on \code{X} otherwise.
#'
#' \strong{Treatment-specific prognostic scores.}
#'
#' If \code{doublepg = TRUE}, the function estimates \code{pg0} and \code{pg1},
#' corresponding to treatment-specific prognostic scores obtained by fitting
#' separate prognostic models within the observed treatment groups.
#'
#' \strong{Censoring-related scores.}
#'
#' If \code{censmod = TRUE}, the function may additionally estimate
#' \code{pscens} and \code{pgcens}, depending on \code{pscens} and
#' \code{pgcens}. These are constructed using covariates \code{c(X, A)} and the
#' supplied \code{event} indicator.
#'
#' \strong{Scaled outputs.}
#'
#' The function also returns scaled versions of the raw scores:
#' \code{ps_sc}, \code{pg0_sc}, \code{pg1_sc}, \code{pscens_sc}, and
#' \code{pgcens_sc}. Propensity-type scores are transformed to the logit scale
#' before standardization.
#'
#' \strong{Parallel fitting.}
#'
#' When \code{superLearn = TRUE}, parallel behavior for SuperLearner-based
#' estimation is controlled by \code{sl_parallel}. On Windows, multicore mode is
#' automatically downgraded to sequential execution.
#'
#' @param data A \code{data.frame} containing all variables required for fitting.
#' @param id Character scalar. Subject identifier column name.
#' @param Y Character scalar. Outcome or follow-up time column name.
#' @param event Character scalar. Event indicator column name, coded \code{1} for
#' event and \code{0} for censoring.
#' @param X Character vector. Covariate names used in the prognostic model(s).
#' @param A Character scalar. Binary treatment indicator column name. Values may
#' be coded as \code{0/1} or \code{-1/1}; values equal to \code{1} are treated as
#' the treated group.
#' @param Xtrt Optional character vector. Covariates used in the treatment
#' propensity model. If \code{NULL}, \code{X} is used.
#' @param doublepg Logical. If \code{TRUE}, estimate treatment-specific
#' prognostic scores \code{pg0} and \code{pg1}. If \code{FALSE}, these columns are
#' returned but remain \code{NA}.
#' @param outer_CV Integer. Number of outer cross-validation folds.
#' @param inner_CV Optional integer. Number of inner cross-validation folds for
#' nested SuperLearner fitting.
#' @param stratifyCV Logical. Whether to request stratified cross-validation when
#' supported by the underlying fitting routine.
#' @param cores Integer. Number of cores requested for fitting.
#' @param tau Optional numeric truncation horizon used when constructing the
#' prediction time grid for restricted mean calculations.
#' @param sl.seed Integer. Random seed used in SuperLearner-based fitting.
#' @param A.SL.library Character vector. SuperLearner library used for treatment
#' propensity and censoring-related models.
#' @param Y.SL.library Character vector. Learners used in \code{survivalSL} for
#' treatment-specific prognostic modeling.
#' @param A.method Character scalar. Risk or loss function passed to
#' \code{CV.SuperLearner()}.
#' @param Y.method Character scalar. Metric passed to \code{survivalSL()}.
#' @param param.tune Optional tuning object passed to \code{survivalSL()}.
#' @param ngrid Integer. Number of grid points used for survival-curve prediction
#' and numerical integration.
#' @param param.weights.fix Optional vector of fixed ensemble weights passed to
#' \code{survivalSL()} when supported.
#' @param param.weights.init Optional vector of initial ensemble weights passed to
#' \code{survivalSL()} when supported.
#' @param optim.method Character scalar. Optimization method passed to
#' \code{survivalSL()}.
#' @param penalty Optional penalty value passed to \code{survivalSL()} or
#' penalized regression routines.
#' @param pgcens Logical. If \code{TRUE} and \code{censmod = TRUE}, estimate the
#' censoring-related prognostic score \code{pgcens}.
#' @param pscens Logical. If \code{TRUE} and \code{censmod = TRUE}, estimate the
#' censoring-related propensity score \code{pscens}.
#' @param censmod Logical. If \code{TRUE}, request censoring-related scores in
#' addition to treatment-based scores.
#' @param maxit Integer. Maximum number of optimization iterations passed to
#' \code{survivalSL()}.
#' @param model.pg Character scalar. Prognostic model family used when
#' \code{superLearn = FALSE}. Must be one of \code{"cox"} or \code{"aft"}.
#' @param standardize Logical. Whether to standardize predictors in glmnet-based
#' fits.
#' @param superLearn Logical. If \code{TRUE}, use SuperLearner / survivalSL-based
#' estimation. Otherwise use glm / glmnet / Cox / AFT alternatives.
#' @param pslink Character scalar. Link function for binomial propensity models.
#' Must be one of \code{"logit"} or \code{"probit"}.
#' @param pglink Character scalar. Distribution used in
#' \code{flexsurv::flexsurvreg()} when \code{model.pg = "aft"}.
#' @param sl_parallel Character scalar. Parallel mode for SuperLearner-based
#' fitting. Must be one of \code{"multicore"} or \code{"seq"}.
#'
#' @return
#' A \code{data.frame} with one row per subject and the following stable columns:
#' \itemize{
#'   \item \code{id}: subject identifier,
#'   \item \code{ps}: treatment propensity score,
#'   \item \code{pg0}, \code{pg1}: treatment-specific prognostic scores,
#'   \item \code{pscens}: censoring-related propensity score,
#'   \item \code{pgcens}: censoring-related prognostic score,
#'   \item \code{ps_sc}, \code{pg0_sc}, \code{pg1_sc}, \code{pscens_sc},
#'   \code{pgcens_sc}: scaled versions of the corresponding scores.
#' }
#'
#'
#' @export
ComputeScores <- function(data, id, Y, event, X, A,
                          Xtrt = NULL,
                          doublepg = TRUE,
                          outer_CV = 5,
                          inner_CV = NULL,
                          stratifyCV = TRUE,
                          cores = 5,
                          tau = NULL,
                          sl.seed = 100,
                          A.SL.library = c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.xgboost"),
                          Y.SL.library = c("LIB_COXen", "LIB_AFTggamma"),
                          A.method = "method.AUC",
                          Y.method = "auc",
                          param.tune = NULL,
                          ngrid = 2000,
                          param.weights.fix = NULL,
                          param.weights.init = NULL,
                          optim.method = "Nelder-Mead",
                          penalty = NULL,
                          pgcens = FALSE,
                          pscens = TRUE,
                          censmod = TRUE,
                          maxit = 1000,
                          model.pg = "cox",
                          standardize = FALSE,
                          superLearn = TRUE,
                          pslink = "logit",
                          pglink = "lognormal",
                          sl_parallel = c("multicore", "seq")) {

  sl_parallel <- match.arg(sl_parallel)
  model.pg <- match.arg(model.pg, c("cox", "aft"))
  pslink <- match.arg(pslink, c("logit", "probit"))

  # ------------------------------------------------------------
  # input validation
  # ------------------------------------------------------------
  stopifnot(is.data.frame(data))

  for (nm in c(id, Y, event, A)) {
    if (!nm %in% names(data)) {
      stop("Column '", nm, "' not found in data.", call. = FALSE)
    }
  }

  if (!is.character(X)) {
    stop("X must be a character vector.", call. = FALSE)
  }
  if (length(X) > 0L && !all(X %in% names(data))) {
    stop("Some X columns not found in data.", call. = FALSE)
  }

  if (!is.null(Xtrt)) {
    if (!is.character(Xtrt)) {
      stop("Xtrt must be NULL or a character vector.", call. = FALSE)
    }
    if (length(Xtrt) > 0L && !all(Xtrt %in% names(data))) {
      stop("Some Xtrt columns not found in data.", call. = FALSE)
    }
  }

  outer_CV <- as.integer(outer_CV)[1]
  if (!is.finite(outer_CV) || is.na(outer_CV) || outer_CV < 2L) {
    stop("outer_CV must be an integer >= 2.", call. = FALSE)
  }

  if (!is.null(inner_CV)) {
    inner_CV <- as.integer(inner_CV)[1]
    if (!is.finite(inner_CV) || is.na(inner_CV) || inner_CV < 2L) {
      stop("inner_CV must be NULL or an integer >= 2.", call. = FALSE)
    }
  }

  cores <- as.integer(cores)[1]
  if (!is.finite(cores) || is.na(cores) || cores < 1L) {
    stop("cores must be an integer >= 1.", call. = FALSE)
  }

  ngrid <- as.integer(ngrid)[1]
  if (!is.finite(ngrid) || is.na(ngrid) || ngrid < 2L) {
    stop("ngrid must be an integer >= 2.", call. = FALSE)
  }

  if (!is.null(tau)) {
    tau <- as.numeric(tau)[1]
    if (!is.finite(tau) || is.na(tau) || tau <= 0) {
      stop("tau must be NULL or a positive numeric value.", call. = FALSE)
    }
  }

  Event <- data[[event]]
  if (is.logical(Event)) Event <- as.integer(Event)
  if (!all(is.na(Event) | Event %in% c(0L, 1L))) {
    stop("'", event, "' must be coded 0/1 with 1=event and 0=censored.", call. = FALSE)
  }

  AA <- data[[A]]
  A_bin <- ifelse(is.na(AA), NA_integer_, ifelse(AA == 1, 1L, 0L))
  if (!all(is.na(A_bin) | A_bin %in% c(0L, 1L))) {
    stop("Treatment variable must be binary, with treated values coded as 1.", call. = FALSE)
  }

  # ------------------------------------------------------------
  # dependency guards
  # ------------------------------------------------------------
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.", call. = FALSE)
  }

  if (superLearn) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required when superLearn = TRUE.", call. = FALSE)
    }
    if (isTRUE(doublepg) && !requireNamespace("survivalSL", quietly = TRUE)) {
      stop("Package 'survivalSL' is required when superLearn = TRUE and doublepg = TRUE.", call. = FALSE)
    }
  } else {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package 'glmnet' is required when superLearn = FALSE.", call. = FALSE)
    }
    if (model.pg == "aft" && isTRUE(doublepg) && !requireNamespace("flexsurv", quietly = TRUE)) {
      stop("Package 'flexsurv' is required for model.pg = 'aft'.", call. = FALSE)
    }
  }

  # ------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------
  tick <- function(...) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic(...)
  }
  tock <- function(...) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()
  }

  prep_df <- function(df) {
    df <- as.data.frame(df)
    df[] <- lapply(df, function(z) if (is.character(z)) factor(z) else z)
    df
  }

  mm_full <- function(df) {
    df <- prep_df(df)
    if (ncol(df) == 0L) {
      return(matrix(numeric(0), nrow = nrow(df), ncol = 0L))
    }
    stats::model.matrix(~ . - 1, data = df)
  }

  safe_scale <- function(z) {
    z <- as.numeric(z)
    ok <- is.finite(z)
    if (sum(ok) <= 1L) return(rep(NA_real_, length(z)))
    sdz <- stats::sd(z[ok])
    if (is.na(sdz) || sdz == 0) return(rep(NA_real_, length(z)))
    outz <- rep(NA_real_, length(z))
    outz[ok] <- as.numeric(scale(z[ok]))
    outz
  }

  safe_logit <- function(p, eps = 1e-6) {
    p <- as.numeric(p)
    p[p <= 0] <- eps
    p[p >= 1] <- 1 - eps
    stats::qlogis(p)
  }

  const_prob <- function(y, n) {
    p <- mean(y, na.rm = TRUE)
    rep(p, n)
  }

  const_value <- function(x, n) {
    rep(mean(x, na.rm = TRUE), n)
  }

  make_time_grid <- function(y, tau, ngrid) {
    ymax <- if (is.null(tau)) max(y, na.rm = TRUE) else tau
    if (!is.finite(ymax) || is.na(ymax) || ymax <= 0) ymax <- 1
    seq(0, ymax, length.out = ngrid)
  }

  normalize_surv_pred <- function(Smat, time_grid, tol = 1e-10) {
    time_grid <- as.numeric(time_grid)

    if (!length(time_grid)) {
      stop("Prediction time grid is empty.", call. = FALSE)
    }
    if (anyNA(time_grid) || !all(is.finite(time_grid))) {
      stop("Prediction time grid must be finite and non-missing.", call. = FALSE)
    }

    Smat <- if (is.null(dim(Smat))) {
      matrix(as.numeric(Smat), nrow = 1L)
    } else {
      as.matrix(Smat)
    }
    storage.mode(Smat) <- "double"

    if (ncol(Smat) != length(time_grid) && nrow(Smat) == length(time_grid)) {
      Smat <- t(Smat)
    }

    if (ncol(Smat) + 1L == length(time_grid) &&
        isTRUE(all.equal(time_grid[1], 0, tolerance = tol))) {
      Smat <- cbind(1, Smat)
    }

    if (ncol(Smat) == length(time_grid) + 1L &&
        all(abs(Smat[, 1] - 1) < 1e-8, na.rm = TRUE)) {
      time_grid <- c(0, time_grid)
    }

    if (ncol(Smat) != length(time_grid)) {
      k <- min(ncol(Smat), length(time_grid))
      warning(
        sprintf(
          "Adjusted survival prediction grid: %d columns in prediction vs %d time points; using first %d aligned points.",
          ncol(Smat), length(time_grid), k
        ),
        call. = FALSE
      )
      Smat <- Smat[, seq_len(k), drop = FALSE]
      time_grid <- time_grid[seq_len(k)]
    }

    ord <- order(time_grid)
    time_grid <- time_grid[ord]
    Smat <- Smat[, ord, drop = FALSE]

    keep <- c(TRUE, diff(time_grid) > tol)
    time_grid <- time_grid[keep]
    Smat <- Smat[, keep, drop = FALSE]

    if (time_grid[1] > tol) {
      time_grid <- c(0, time_grid)
      Smat <- cbind(1, Smat)
    } else {
      time_grid[1] <- 0
      Smat[, 1] <- 1
    }

    Smat[!is.finite(Smat)] <- NA_real_
    Smat <- pmin(pmax(Smat, 0), 1)

    Smat <- t(apply(Smat, 1, function(z) {
      if (all(is.na(z))) return(rep(NA_real_, length(z)))
      idx <- which(!is.na(z))
      z <- stats::approx(
        x = idx, y = z[idx], xout = seq_along(z),
        method = "linear", rule = 2
      )$y
      cummin(z)
    }))

    list(Smat = Smat, time_grid = time_grid)
  }

  pred_mean <- function(Smat, time_grid) {
    obj <- normalize_surv_pred(Smat, time_grid)
    Smat <- obj$Smat
    time_grid <- obj$time_grid

    if (length(time_grid) < 2L) {
      return(rep(0, nrow(Smat)))
    }

    dt <- diff(time_grid)
    mids <- (Smat[, -1, drop = FALSE] + Smat[, -ncol(Smat), drop = FALSE]) / 2
    rowSums(mids * matrix(dt, nrow = nrow(Smat), ncol = length(dt), byrow = TRUE))
  }

  predict_survival_mean <- function(fit, newdata, newtimes) {
    pred <- .tag_error(
      predict(fit, newdata = newdata, newtimes = newtimes),
      "ComputeScores: predict survivalSL"
    )

    if (is.null(pred$predictions) || is.null(pred$predictions$sl) || is.null(pred$times)) {
      stop(
        "survivalSL::predict() did not return both 'predictions$sl' and 'times'.",
        call. = FALSE
      )
    }

    pred_mean(pred$predictions$sl, pred$times)
  }

  sl_parallel_resolved <- if (sl_parallel == "multicore" && .Platform$OS.type == "windows") {
    "seq"
  } else {
    sl_parallel
  }

  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  options(mc.cores = max(1L, cores))

  # ------------------------------------------------------------
  # construct data objects
  # ------------------------------------------------------------
  X <- unique(setdiff(X, A))
  XX <- prep_df(data[, X, drop = FALSE])
  XXtrt <- if (!is.null(Xtrt)) prep_df(data[, Xtrt, drop = FALSE]) else XX
  Xcens <- unique(c(X, A))
  XXcens <- prep_df(data[, Xcens, drop = FALSE])

  y <- A_bin
  YY <- data[[Y]]
  Id <- data[[id]]

  loc1 <- which(A_bin == 1L)
  loc0 <- which(A_bin == 0L)

  out <- data.frame(
    tmp_id     = Id,
    ps         = NA_real_,
    pg0        = NA_real_,
    pg1        = NA_real_,
    pscens     = NA_real_,
    pgcens     = NA_real_,
    ps_sc      = NA_real_,
    pg0_sc     = NA_real_,
    pg1_sc     = NA_real_,
    pscens_sc  = NA_real_,
    pgcens_sc  = NA_real_
  )
  names(out)[1] <- id

  innerCvControl_value <- if (!is.null(inner_CV)) {
    rep(list(list(V = inner_CV, stratifyCV = stratifyCV)), outer_CV)
  } else {
    NULL
  }

  # ------------------------------------------------------------
  # 1) Treatment propensity score: ps
  # ------------------------------------------------------------
  if (superLearn) {
    if (ncol(XXtrt) == 0L) {
      out$ps <- const_prob(y, nrow(data))
    } else {
      tick("Treatment PS: CV.SuperLearner")
      set.seed(sl.seed, "L'Ecuyer-CMRG")
      sl_ps <- .tag_error(
        SuperLearner::CV.SuperLearner(
          Y = y,
          X = XXtrt,
          family = stats::binomial(link = pslink),
          method = A.method,
          SL.library = A.SL.library,
          cvControl = list(V = outer_CV, stratifyCV = stratifyCV),
          innerCvControl = innerCvControl_value,
          parallel = sl_parallel_resolved,
          verbose = TRUE,
          env = getNamespace("SuperLearner")
        ),
        "ComputeScores: treatment PS CV.SuperLearner"
      )
      out$ps <- as.vector(sl_ps$SL.predict)
      tock()
      rm(sl_ps)
      gc()
    }
  } else {
    if (pslink == "logit") {
      xtrt_mm <- mm_full(XXtrt)

      if (ncol(xtrt_mm) == 0L) {
        out$ps <- const_prob(y, nrow(data))
      } else if (ncol(xtrt_mm) > 1L) {
        cv.fit.ps <- .tag_error(
          glmnet::cv.glmnet(
            x = xtrt_mm,
            y = y,
            family = "binomial",
            alpha = 1,
            nfolds = outer_CV,
            standardize = standardize
          ),
          "ComputeScores: treatment PS cv.glmnet"
        )
        out$ps <- as.vector(
          stats::predict(cv.fit.ps, newx = xtrt_mm, s = cv.fit.ps$lambda.min, type = "response")
        )
        rm(cv.fit.ps)
        gc()
      } else {
        df_ps <- data.frame(y = y, x1 = xtrt_mm[, 1])
        fit_ps <- .tag_error(
          stats::glm(y ~ x1, data = df_ps, family = stats::binomial(link = "logit")),
          "ComputeScores: treatment PS glm"
        )
        out$ps <- as.vector(stats::predict(fit_ps, newdata = df_ps, type = "response"))
        rm(fit_ps, df_ps)
        gc()
      }
      rm(xtrt_mm)
      gc()
    } else {
      if (ncol(XXtrt) == 0L) {
        out$ps <- const_prob(y, nrow(data))
      } else {
        df_ps <- data.frame(y = y, XXtrt)
        fit_ps <- .tag_error(
          stats::glm(y ~ ., data = df_ps, family = stats::binomial(link = pslink)),
          "ComputeScores: treatment PS glm nonlogit"
        )
        out$ps <- as.vector(stats::predict(fit_ps, newdata = XXtrt, type = "response"))
        rm(fit_ps, df_ps)
        gc()
      }
    }
  }

  # ------------------------------------------------------------
  # 2) Treatment prognostic scores: pg0, pg1
  # ------------------------------------------------------------
  if (isTRUE(doublepg)) {
    if (superLearn) {
      surv_formula <- if (length(X) > 0L) {
        stats::as.formula(
          paste0("survival::Surv(", Y, ", ", event, ") ~ ", paste(X, collapse = " + "))
        )
      } else {
        stats::as.formula(paste0("survival::Surv(", Y, ", ", event, ") ~ 1"))
      }

      if (length(loc1) > 0L) {
        survdata1 <- prep_df(data[loc1, c(Y, event, X), drop = FALSE])
        tick("Treatment PG1: survivalSL")
        slres1 <- .tag_error(
          survivalSL::survivalSL(
            formula = surv_formula,
            methods = Y.SL.library,
            metric = Y.method,
            data = survdata1,
            cv = outer_CV,
            param.tune = param.tune,
            seed = sl.seed,
            param.weights.fix = param.weights.fix,
            param.weights.init = param.weights.init,
            maxit = maxit,
            penalty = penalty,
            show_progress = TRUE
          ),
          "ComputeScores: treatment PG1 survivalSL"
        )
        grid1 <- make_time_grid(YY, tau, ngrid)
        out$pg1 <- predict_survival_mean(slres1, newdata = XX, newtimes = grid1)
        tock()
        rm(slres1, survdata1, grid1)
        gc()
      }

      if (length(loc0) > 0L) {
        survdata0 <- prep_df(data[loc0, c(Y, event, X), drop = FALSE])
        tick("Treatment PG0: survivalSL")
        slres0 <- .tag_error(
          survivalSL::survivalSL(
            formula = surv_formula,
            methods = Y.SL.library,
            metric = Y.method,
            data = survdata0,
            cv = outer_CV,
            param.tune = param.tune,
            seed = sl.seed,
            param.weights.fix = param.weights.fix,
            param.weights.init = param.weights.init,
            maxit = maxit,
            penalty = penalty,
            show_progress = TRUE
          ),
          "ComputeScores: treatment PG0 survivalSL"
        )
        grid0 <- make_time_grid(YY, tau, ngrid)
        out$pg0 <- predict_survival_mean(slres0, newdata = XX, newtimes = grid0)
        tock()
        rm(slres0, survdata0, grid0)
        gc()
      }

    } else {
      if (model.pg == "cox") {
        x_mm <- mm_full(XX)

        if (length(loc1) > 0L) {
          if (ncol(x_mm) == 0L) {
            out$pg1 <- rep(0, nrow(data))
          } else {
            Y1 <- survival::Surv(YY[loc1], Event[loc1])
            X1 <- x_mm[loc1, , drop = FALSE]

            if (ncol(X1) > 1L) {
              cv.fit1 <- .tag_error(
                glmnet::cv.glmnet(
                  x = X1,
                  y = Y1,
                  family = "cox",
                  alpha = 1,
                  nfolds = outer_CV,
                  standardize = standardize
                ),
                "ComputeScores: treatment PG1 cv.glmnet cox"
              )
              out$pg1 <- as.vector(
                stats::predict(cv.fit1, newx = x_mm, s = cv.fit1$lambda.min, type = "link")
              )
              rm(cv.fit1)
              gc()
            } else {
              df1 <- data.frame(Y = YY[loc1], event = Event[loc1], x1 = X1[, 1])
              fit1 <- .tag_error(
                survival::coxph(survival::Surv(Y, event) ~ x1, data = df1),
                "ComputeScores: treatment PG1 coxph"
              )
              out$pg1 <- as.vector(
                stats::predict(fit1, newdata = data.frame(x1 = x_mm[, 1]), type = "lp")
              )
              rm(df1, fit1)
              gc()
            }
          }
        }

        if (length(loc0) > 0L) {
          if (ncol(x_mm) == 0L) {
            out$pg0 <- rep(0, nrow(data))
          } else {
            Y0 <- survival::Surv(YY[loc0], Event[loc0])
            X0 <- x_mm[loc0, , drop = FALSE]

            if (ncol(X0) > 1L) {
              cv.fit0 <- .tag_error(
                glmnet::cv.glmnet(
                  x = X0,
                  y = Y0,
                  family = "cox",
                  alpha = 1,
                  nfolds = outer_CV,
                  standardize = standardize
                ),
                "ComputeScores: treatment PG0 cv.glmnet cox"
              )
              out$pg0 <- as.vector(
                stats::predict(cv.fit0, newx = x_mm, s = cv.fit0$lambda.min, type = "link")
              )
              rm(cv.fit0)
              gc()
            } else {
              df0 <- data.frame(Y = YY[loc0], event = Event[loc0], x1 = X0[, 1])
              fit0 <- .tag_error(
                survival::coxph(survival::Surv(Y, event) ~ x1, data = df0),
                "ComputeScores: treatment PG0 coxph"
              )
              out$pg0 <- as.vector(
                stats::predict(fit0, newdata = data.frame(x1 = x_mm[, 1]), type = "lp")
              )
              rm(df0, fit0)
              gc()
            }
          }
        }

        rm(x_mm)
        gc()

      } else if (model.pg == "aft") {
        form_aft <- if (length(X) > 0L) {
          survival::Surv(Y, event) ~ .
        } else {
          survival::Surv(Y, event) ~ 1
        }

        if (length(loc1) > 0L) {
          data1 <- data.frame(Y = YY[loc1], event = Event[loc1], XX[loc1, , drop = FALSE])
          fit1 <- .tag_error(
            flexsurv::flexsurvreg(
              form_aft,
              data = data1,
              dist = pglink
            ),
            "ComputeScores: treatment PG1 flexsurvreg"
          )
          out$pg1 <- as.numeric(stats::predict(fit1, newdata = XX, type = "mean"))
          rm(data1, fit1)
          gc()
        }

        if (length(loc0) > 0L) {
          data0 <- data.frame(Y = YY[loc0], event = Event[loc0], XX[loc0, , drop = FALSE])
          fit0 <- .tag_error(
            flexsurv::flexsurvreg(
              form_aft,
              data = data0,
              dist = pglink
            ),
            "ComputeScores: treatment PG0 flexsurvreg"
          )
          out$pg0 <- as.numeric(stats::predict(fit0, newdata = XX, type = "mean"))
          rm(data0, fit0)
          gc()
        }
      }
    }
  }

  # ------------------------------------------------------------
  # 3) Censoring-related scores using c(X, A)
  # ------------------------------------------------------------
  if (isTRUE(censmod) && (isTRUE(pscens) || isTRUE(pgcens))) {

    if (isTRUE(pscens)) {
      if (superLearn) {
        if (ncol(XXcens) == 0L) {
          out$pscens <- const_prob(Event, nrow(data))
        } else {
          tick("Censoring PS: CV.SuperLearner")
          set.seed(sl.seed, "L'Ecuyer-CMRG")
          sl_cens <- .tag_error(
            SuperLearner::CV.SuperLearner(
              Y = Event,
              X = XXcens,
              family = stats::binomial(link = pslink),
              method = A.method,
              SL.library = A.SL.library,
              cvControl = list(V = outer_CV, stratifyCV = stratifyCV),
              innerCvControl = innerCvControl_value,
              parallel = sl_parallel_resolved,
              verbose = TRUE,
              env = getNamespace("SuperLearner")
            ),
            "ComputeScores: censoring PS CV.SuperLearner"
          )
          out$pscens <- as.vector(sl_cens$SL.predict)
          tock()
          rm(sl_cens)
          gc()
        }
      } else {
        if (pslink == "logit") {
          x_mm_cens <- mm_full(XXcens)

          if (ncol(x_mm_cens) == 0L) {
            out$pscens <- const_prob(Event, nrow(data))
          } else if (ncol(x_mm_cens) > 1L) {
            cv.fit.cens <- .tag_error(
              glmnet::cv.glmnet(
                x = x_mm_cens,
                y = Event,
                family = "binomial",
                alpha = 1,
                nfolds = outer_CV,
                standardize = standardize
              ),
              "ComputeScores: censoring PS cv.glmnet"
            )
            out$pscens <- as.vector(
              stats::predict(cv.fit.cens, newx = x_mm_cens, s = cv.fit.cens$lambda.min, type = "response")
            )
            rm(cv.fit.cens)
            gc()
          } else {
            df_cens <- data.frame(Event = Event, x1 = x_mm_cens[, 1])
            fit_cens <- .tag_error(
              stats::glm(Event ~ x1, data = df_cens, family = stats::binomial(link = "logit")),
              "ComputeScores: censoring PS glm"
            )
            out$pscens <- as.vector(stats::predict(fit_cens, newdata = df_cens, type = "response"))
            rm(fit_cens, df_cens)
            gc()
          }
          rm(x_mm_cens)
          gc()
        } else {
          if (ncol(XXcens) == 0L) {
            out$pscens <- const_prob(Event, nrow(data))
          } else {
            df_cens <- data.frame(Event = Event, XXcens)
            fit_cens <- .tag_error(
              stats::glm(Event ~ ., data = df_cens, family = stats::binomial(link = pslink)),
              "ComputeScores: censoring PS glm nonlogit"
            )
            out$pscens <- as.vector(stats::predict(fit_cens, newdata = XXcens, type = "response"))
            rm(fit_cens, df_cens)
            gc()
          }
        }
      }
    }

    if (isTRUE(pgcens)) {
      unc_idx <- which(Event == 1L)

      if (length(unc_idx) > 0L) {
        if (superLearn) {
          if (ncol(XXcens) == 0L) {
            out$pgcens <- const_value(YY[unc_idx], nrow(data))
          } else {
            tick("Censoring PG: SuperLearner")
            set.seed(sl.seed, "L'Ecuyer-CMRG")
            sl_pgcens <- .tag_error(
              SuperLearner::SuperLearner(
                Y = YY[unc_idx],
                X = XXcens[unc_idx, , drop = FALSE],
                family = stats::gaussian(),
                SL.library = A.SL.library,
                cvControl = list(V = outer_CV, stratifyCV = FALSE),
                verbose = TRUE,
                env = getNamespace("SuperLearner")
              ),
              "ComputeScores: censoring PG SuperLearner"
            )
            out$pgcens <- as.vector(
              stats::predict(sl_pgcens, newdata = XXcens)$pred
            )
            tock()
            rm(sl_pgcens)
            gc()
          }
        } else {
          XA_mm <- mm_full(XXcens)

          if (ncol(XA_mm) == 0L) {
            out$pgcens <- const_value(YY[unc_idx], nrow(data))
          } else if (ncol(XA_mm) > 1L) {
            cv.fit.pg <- .tag_error(
              glmnet::cv.glmnet(
                x = XA_mm[unc_idx, , drop = FALSE],
                y = YY[unc_idx],
                family = "gaussian",
                alpha = 1,
                nfolds = outer_CV,
                standardize = standardize
              ),
              "ComputeScores: censoring PG cv.glmnet"
            )
            out$pgcens <- as.vector(
              stats::predict(cv.fit.pg, newx = XA_mm, s = cv.fit.pg$lambda.min, type = "response")
            )
            rm(cv.fit.pg)
            gc()
          } else {
            df_pg <- data.frame(Y = YY[unc_idx], x1 = XA_mm[unc_idx, 1])
            fit_pg <- .tag_error(
              stats::lm(Y ~ x1, data = df_pg),
              "ComputeScores: censoring PG lm"
            )
            out$pgcens <- as.vector(
              stats::predict(fit_pg, newdata = data.frame(x1 = XA_mm[, 1]))
            )
            rm(fit_pg, df_pg)
            gc()
          }
          rm(XA_mm)
          gc()
        }
      }
    }
  }

  # ------------------------------------------------------------
  # 4) Scaled versions
  # ------------------------------------------------------------
  out$ps_sc     <- safe_scale(safe_logit(out$ps))
  out$pg0_sc    <- safe_scale(out$pg0)
  out$pg1_sc    <- safe_scale(out$pg1)
  out$pscens_sc <- safe_scale(safe_logit(out$pscens))
  out$pgcens_sc <- safe_scale(out$pgcens)

  out
}

#' Compute stage-specific treatment, prognostic, and optional censoring scores
#'
#' @description
#' Computes and attaches stage-specific score summaries for a two-stage treatment setting
#' by calling \code{\link{ComputeScores}} separately at stage 2 and stage 1.
#'
#' The function is a wrapper that:
#' \enumerate{
#'   \item restricts to subjects with \code{eta2 == 1} and computes stage-2 scores using
#'   \code{Y2.var}, \code{A2.var}, \code{names.var2}, and \code{Xtrt2};
#'   \item computes stage-1 scores on the full cohort using either \code{OY.var} or
#'   \code{Y1.var} depending on \code{adjustdelta1}, together with \code{A1.var},
#'   \code{names.var1}, and \code{Xtrt1};
#'   \item renames the outputs from \code{ComputeScores()} into stage-specific columns
#'   and merges them back into the original dataset by subject ID.
#' }
#'
#' @details
#' \strong{Scores returned by \code{ComputeScores()}.}
#'
#' For each call, \code{ComputeScores()} returns a fixed set of score columns:
#' \itemize{
#'   \item \code{ps}: treatment propensity score,
#'   \item \code{pg0}, \code{pg1}: treatment-specific prognostic scores,
#'   \item \code{pscens}: censoring propensity score,
#'   \item \code{pgcens}: censoring prognostic score,
#'   \item \code{ps_sc}, \code{pg0_sc}, \code{pg1_sc}, \code{pscens_sc}, \code{pgcens_sc}:
#'   scaled versions of the corresponding raw scores.
#' }
#'
#' This wrapper renames those outputs to stage-specific names and attaches them to
#' \code{data}. For stage 1, the suffix \code{1} is used; for stage 2, the suffix
#' \code{2} is used.
#'
#' \strong{Stage 2.}
#'
#' Stage-2 scores are computed only among subjects satisfying \code{data[[eta2.var]] == 1}.
#' These scores are then merged back into the full dataset. Subjects who do not enter
#' stage 2 receive \code{NA} for all stage-2 score columns.
#'
#' \strong{Stage 1.}
#'
#' Stage-1 scores are computed on the full dataset. By default, the stage-1 scoring call
#' uses \code{Y = OY.var} and \code{event = delta.var}. If \code{adjustdelta1 = TRUE},
#' the function instead uses \code{Y = Y1.var} and a modified event indicator
#' \code{deltaadj}, where \code{deltaadj} is initialized as \code{delta.var} and then set
#' to 0 for subjects with \code{eta2 == 1} and \code{delta == 1}.
#'
#' \strong{Treatment and censoring covariates.}
#'
#' The prognostic model covariates are supplied through \code{names.var1} and
#' \code{names.var2}. The treatment propensity model covariates are supplied separately
#' through \code{Xtrt1} and \code{Xtrt2}. If \code{Xtrt1} or \code{Xtrt2} is \code{NULL},
#' then \code{ComputeScores()} uses the corresponding prognostic covariates.
#'
#' \strong{Censoring-related scores.}
#'
#' If \code{censmod = TRUE}, the wrapper also requests censoring-related scores from
#' \code{ComputeScores()}. The arguments \code{pscens} and \code{pgcens} determine whether
#' censoring propensity and censoring prognostic scores are actively estimated. Since
#' \code{ComputeScores()} returns a fixed output structure, the corresponding columns are
#' still present in the returned data even when those components are not estimated; in
#' such cases they are typically \code{NA}.
#'
#' \strong{Treatment-specific prognostic scores.}
#'
#' If \code{doublepg = TRUE}, the wrapper additionally creates convenience variables
#' comparing the observed-treatment and opposite-treatment prognostic scores:
#' \itemize{
#'   \item \code{pg1ct}, \code{pg1tc} for stage 1,
#'   \item \code{pg2ct}, \code{pg2tc} for stage 2.
#' }
#' These are constructed from the standardized prognostic scores:
#' \code{pg01}, \code{pg11}, \code{pg02}, and \code{pg12}.
#'
#' \strong{No-op behavior.}
#'
#' If \code{useds = FALSE}, the function returns \code{data} unchanged.
#'
#' @param data A \code{data.frame} containing subject identifiers, stage indicators,
#' outcomes, treatment variables, and covariates required for stage-1 and stage-2 score
#' estimation.
#'
#' @param id.var Character scalar. Name of the subject identifier column.
#'
#' @param eta2.var Character scalar. Name of the stage-2 entry indicator column, where
#' \code{1} denotes entry into stage 2 and \code{0} denotes no entry.
#'
#' @param Y1.var Character scalar. Name of the stage-1 outcome or time variable used when
#' \code{adjustdelta1 = TRUE}.
#'
#' @param Y2.var Character scalar. Name of the stage-2 outcome or time variable.
#'
#' @param delta.var Character scalar. Name of the event indicator variable used in the
#' survival or censoring models.
#'
#' @param OY.var Character scalar. Name of the overall outcome or time variable used for
#' stage-1 scoring when \code{adjustdelta1 = FALSE}.
#'
#' @param A1.var Character scalar. Name of the stage-1 treatment indicator variable.
#'
#' @param A2.var Character scalar. Name of the stage-2 treatment indicator variable.
#'
#' @param names.var1 Character vector. Covariate names used in the stage-1 prognostic
#' score model.
#'
#' @param names.var2 Character vector. Covariate names used in the stage-2 prognostic
#' score model.
#'
#' @param Xtrt1 Character vector or \code{NULL}. Covariate names used in the stage-1
#' treatment propensity model. If \code{NULL}, \code{ComputeScores()} uses
#' \code{names.var1}.
#'
#' @param Xtrt2 Character vector or \code{NULL}. Covariate names used in the stage-2
#' treatment propensity model. If \code{NULL}, \code{ComputeScores()} uses
#' \code{names.var2}.
#'
#' @param useds Logical. If \code{TRUE}, compute and merge the stage-specific scores.
#' If \code{FALSE}, return \code{data} unchanged.
#'
#' @param cores Integer. Number of cores passed to \code{ComputeScores()} for model fitting.
#'
#' @param tau Optional numeric truncation horizon passed to \code{ComputeScores()} for
#' restricted mean prediction or time-grid construction.
#'
#' @param sl.seed Integer. Random seed passed to \code{ComputeScores()}.
#'
#' @param A.SL.library1 Character vector. SuperLearner library for the stage-1 treatment
#' propensity model.
#'
#' @param A.SL.library2 Character vector. SuperLearner library for the stage-2 treatment
#' propensity model.
#'
#' @param Y.SL.library Character vector. Learners used for prognostic survival modeling
#' inside \code{ComputeScores()}.
#'
#' @param A.method Optional character scalar. Performance metric passed to
#' \code{ComputeScores()} for treatment propensity estimation.
#'
#' @param Y.method Optional character scalar. Performance metric passed to
#' \code{ComputeScores()} for prognostic survival estimation.
#'
#' @param param.weights.fix Optional numeric vector. Fixed ensemble weights passed to
#' \code{ComputeScores()} when supported by the underlying learner.
#'
#' @param param.weights.init Optional numeric vector. Initial ensemble weights passed to
#' \code{ComputeScores()} when supported by the underlying learner.
#'
#' @param optim.method Character scalar or \code{NULL}. Optimization method forwarded to
#' \code{ComputeScores()}.
#'
#' @param stratifyCV Logical. Passed to \code{ComputeScores()}. If \code{TRUE},
#' cross-validation folds are stratified when supported by the underlying fitting
#' procedure.
#'
#' @param maxit Integer. Maximum number of optimization iterations passed to
#' \code{ComputeScores()}.
#'
#' @param penalty1 Optional tuning parameter or penalty value passed to
#' \code{ComputeScores()} for stage-1 prognostic estimation.
#'
#' @param penalty2 Optional tuning parameter or penalty value passed to
#' \code{ComputeScores()} for stage-2 prognostic estimation.
#'
#' @param ngrid Integer. Number of grid points used by \code{ComputeScores()} when
#' approximating restricted means or evaluating predicted survival curves.
#'
#' @param censmod Logical. If \code{TRUE}, request censoring-related scores from
#' \code{ComputeScores()} in addition to treatment propensity and treatment prognostic
#' scores.
#'
#' @param pscens Logical. If \code{TRUE} and \code{censmod = TRUE}, estimate censoring
#' propensity scores within \code{ComputeScores()}.
#'
#' @param pgcens Logical. If \code{TRUE} and \code{censmod = TRUE}, estimate censoring
#' prognostic scores within \code{ComputeScores()}.
#'
#' @param doublepg Logical. Passed to \code{ComputeScores()}. If \code{TRUE}, estimate
#' treatment-specific prognostic scores separately by treatment arm. If \code{FALSE},
#' the returned prognostic components may be partially unestimated and therefore remain
#' \code{NA}.
#'
#' @param param.tune Optional list or tuning object passed to \code{ComputeScores()} for
#' learner-specific tuning.
#'
#' @param adjustdelta1 Logical. If \code{TRUE}, construct an adjusted stage-1 event
#' indicator \code{deltaadj} and use \code{Y1.var} instead of \code{OY.var} in the
#' stage-1 scoring call.
#'
#' @param plotps Logical. If \code{TRUE}, plot the raw treatment propensity score
#' distribution at each stage using \code{\link{propensityplot}}, when available.
#'
#' @param model.pg Character scalar. Prognostic model type passed to \code{ComputeScores()}.
#' Currently intended values are \code{"cox"} and \code{"aft"}.
#'
#' @param standardize Logical. Passed to \code{ComputeScores()}. If \code{TRUE},
#' standardize covariates for penalized regression fits when applicable.
#'
#' @param superLearn Logical. Passed to \code{ComputeScores()}. If \code{TRUE}, use
#' SuperLearner-based fitting; otherwise use the parametric or penalized alternatives
#' implemented there.
#'
#' @param pslink Character scalar. Link function for binomial propensity models passed to
#' \code{ComputeScores()}, typically \code{"logit"} or \code{"probit"}.
#'
#' @param pglink Character scalar. Distribution used when \code{model.pg = "aft"} inside
#' \code{ComputeScores()}, for example \code{"exponential"}, \code{"weibull"},
#' \code{"lognormal"}, or \code{"loglogistic"}.
#'
#' @param sl_parallel Character scalar. Parallel mode passed to \code{ComputeScores()}
#' for SuperLearner-based fitting. Must be one of \code{"multicore"} or \code{"seq"}.
#'
#' @return
#' A \code{data.frame} equal to \code{data} augmented with stage-specific score columns.
#' If \code{useds = FALSE}, the original \code{data} is returned unchanged.
#'
#' The following columns are attached for stage 1:
#' \itemize{
#'   \item \code{probps1}: raw treatment propensity score,
#'   \item \code{prog01}, \code{prog11}: raw treatment-specific prognostic scores,
#'   \item \code{probcens1}: raw censoring propensity score,
#'   \item \code{progcens1}: raw censoring prognostic score,
#'   \item \code{ps1}, \code{pg01}, \code{pg11}, \code{pscens1}, \code{pgcens1}:
#'   scaled versions of the above scores.
#' }
#'
#' The analogous columns
#' \code{probps2}, \code{prog02}, \code{prog12}, \code{probcens2}, \code{progcens2},
#' \code{ps2}, \code{pg02}, \code{pg12}, \code{pscens2}, and \code{pgcens2} are attached
#' for stage 2. Subjects with \code{eta2 == 0} receive \code{NA} for all stage-2 score
#' columns.
#'
#' If \code{doublepg = TRUE}, the convenience columns \code{pg1ct}, \code{pg1tc},
#' \code{pg2ct}, and \code{pg2tc} are also added.
#'
#' @seealso \code{\link{ComputeScores}}, \code{\link{propensityplot}}
#' @export


get_doublescores <- function(
    data,
    id.var, eta2.var,
    Y1.var, Y2.var,
    delta.var, OY.var,
    A1.var, A2.var,
    names.var1, names.var2,
    Xtrt1 = NULL,
    Xtrt2 = NULL,
    useds         = FALSE,
    cores         = 1,
    tau,
    sl.seed       = 123,
    A.SL.library1,
    A.SL.library2,
    Y.SL.library,
    A.method      = NULL,
    Y.method      = NULL,
    param.weights.fix   = NULL,
    param.weights.init  = NULL,
    optim.method        = NULL,
    stratifyCV    = TRUE,
    maxit         = 1000,
    penalty1      = NULL,
    penalty2      = NULL,
    ngrid         = 50,
    censmod       = TRUE,
    pscens        = TRUE,
    pgcens        = TRUE,
    doublepg      = TRUE,
    param.tune    = NULL,
    adjustdelta1  = FALSE,
    plotps        = FALSE,
    model.pg      = "cox",
    standardize   = FALSE,
    superLearn    = TRUE,
    pslink        = "logit",
    pglink        = "lognormal",
    sl_parallel   = c("multicore", "seq")
) {

  sl_parallel <- match.arg(sl_parallel)

  stopifnot(is.data.frame(data))

  if (!isTRUE(useds)) return(data)

  if (!exists("ComputeScores", mode = "function")) {
    stop("ComputeScores() not found in the current environment.", call. = FALSE)
  }

  req_cols <- unique(c(
    id.var, eta2.var, Y1.var, Y2.var, delta.var, OY.var, A1.var, A2.var,
    names.var1, names.var2, Xtrt1, Xtrt2
  ))
  req_cols <- req_cols[!is.na(req_cols) & nzchar(req_cols)]

  miss <- setdiff(req_cols, names(data))
  if (length(miss) > 0L) {
    stop("Missing required columns in data: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # ------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------
  tick <- function(...) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic(...)
  }
  tock <- function(...) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()
  }

  expected_cs_names <- function(id.var) {
    c(
      id.var,
      "ps", "pg0", "pg1", "pscens", "pgcens",
      "ps_sc", "pg0_sc", "pg1_sc", "pscens_sc", "pgcens_sc"
    )
  }

  rename_stage_scores <- function(ds, id.var, stage = c("1", "2")) {
    stage <- match.arg(stage)
    ds <- as.data.frame(ds)

    if (!id.var %in% names(ds)) {
      names(ds)[1] <- id.var
    }

    want <- expected_cs_names(id.var)
    if (!identical(names(ds), want)) {
      if (ncol(ds) != length(want)) {
        stop(
          "ComputeScores output has ", ncol(ds), " columns, but ",
          length(want), " were expected.",
          call. = FALSE
        )
      }
      names(ds) <- want
    }

    names(ds) <- c(
      id.var,
      paste0("probps",  stage),
      paste0("prog0",   stage),
      paste0("prog1",   stage),
      paste0("probcens",stage),
      paste0("progcens",stage),
      paste0("ps",      stage),
      paste0("pg0",     stage),
      paste0("pg1",     stage),
      paste0("pscens",  stage),
      paste0("pgcens",  stage)
    )

    ds
  }

  attach_by_id <- function(df, add, id.var) {
    add <- as.data.frame(add)

    if (anyDuplicated(add[[id.var]]) > 0L) {
      stop("Duplicated IDs found in score output for ", id.var, ".", call. = FALSE)
    }

    idx <- match(df[[id.var]], add[[id.var]])
    add_cols <- setdiff(names(add), id.var)

    for (nm in add_cols) {
      df[[nm]] <- add[[nm]][idx]
    }

    df
  }

  add_empty_stage_cols <- function(df, stage = c("1", "2")) {
    stage <- match.arg(stage)
    nm <- c(
      paste0("probps",  stage),
      paste0("prog0",   stage),
      paste0("prog1",   stage),
      paste0("probcens",stage),
      paste0("progcens",stage),
      paste0("ps",      stage),
      paste0("pg0",     stage),
      paste0("pg1",     stage),
      paste0("pscens",  stage),
      paste0("pgcens",  stage)
    )
    for (x in nm) df[[x]] <- NA_real_
    df
  }

  maybe_plot_ps <- function(dat, ps_col, A_col) {
    if (isTRUE(plotps) && exists("propensityplot", mode = "function")) {
      print(propensityplot(ps = dat[[ps_col]], A = dat[[A_col]]))
    }
  }

  # ------------------------------------------------------------
  # working copy
  # ------------------------------------------------------------
  df <- data

  # optional delta adjustment for stage 1
  if (isTRUE(adjustdelta1)) {
    df$deltaadj <- df[[delta.var]]
    df$deltaadj[df[[eta2.var]] == 1 & df[[delta.var]] == 1] <- 0
  }

  # ------------------------------------------------------------
  # Stage 2
  # ------------------------------------------------------------
  df2 <- df[df[[eta2.var]] == 1, , drop = FALSE]

  if (nrow(df2) > 0L) {
    tick("Time to compute stage-2 scores")

    ds2 <- ComputeScores(
      data         = df2,
      id           = id.var,
      Y            = Y2.var,
      event        = delta.var,
      X            = names.var2,
      A            = A2.var,
      Xtrt         = Xtrt2,
      doublepg     = doublepg,
      outer_CV     = 5,
      inner_CV     = 5,
      stratifyCV   = stratifyCV,
      cores        = cores,
      tau          = tau,
      sl.seed      = sl.seed,
      A.SL.library = A.SL.library2,
      Y.SL.library = Y.SL.library,
      A.method     = A.method,
      Y.method     = Y.method,
      param.tune   = param.tune,
      ngrid        = ngrid,
      param.weights.fix  = param.weights.fix,
      param.weights.init = param.weights.init,
      optim.method = optim.method,
      penalty      = penalty2,
      pgcens       = pgcens,
      pscens       = pscens,
      censmod      = censmod,
      maxit        = maxit,
      model.pg     = model.pg,
      standardize  = standardize,
      superLearn   = superLearn,
      pslink       = pslink,
      pglink       = pglink,
      sl_parallel  = sl_parallel
    )

    tock()

    ds2 <- rename_stage_scores(ds2, id.var = id.var, stage = "2")
    df  <- attach_by_id(df, ds2, id.var = id.var)

    df2_plot <- df[df[[eta2.var]] == 1, , drop = FALSE]
    maybe_plot_ps(df2_plot, "probps2", A2.var)

  } else {
    message("No stage-2 entrants (", eta2.var, " == 1); stage-2 score columns set to NA.")
    df <- add_empty_stage_cols(df, stage = "2")
  }

  # ------------------------------------------------------------
  # Stage 1
  # ------------------------------------------------------------
  tick("Time to compute stage-1 scores")

  ds1 <- ComputeScores(
    data         = df,
    id           = id.var,
    Y            = if (isTRUE(adjustdelta1)) Y1.var else OY.var,
    event        = if (isTRUE(adjustdelta1)) "deltaadj" else delta.var,
    X            = names.var1,
    A            = A1.var,
    Xtrt         = Xtrt1,
    doublepg     = doublepg,
    outer_CV     = 5,
    inner_CV     = 5,
    stratifyCV   = stratifyCV,
    cores        = cores,
    tau          = tau,
    sl.seed      = sl.seed,
    A.SL.library = A.SL.library1,
    Y.SL.library = Y.SL.library,
    A.method     = A.method,
    Y.method     = Y.method,
    param.tune   = param.tune,
    ngrid        = ngrid,
    param.weights.fix  = param.weights.fix,
    param.weights.init = param.weights.init,
    optim.method = optim.method,
    penalty      = penalty1,
    pgcens       = pgcens,
    pscens       = pscens,
    censmod      = censmod,
    maxit        = maxit,
    model.pg     = model.pg,
    standardize  = standardize,
    superLearn   = superLearn,
    pslink       = pslink,
    pglink       = pglink,
    sl_parallel  = sl_parallel
  )

  tock()

  ds1 <- rename_stage_scores(ds1, id.var = id.var, stage = "1")
  df  <- attach_by_id(df, ds1, id.var = id.var)

  maybe_plot_ps(df, "probps1", A1.var)

  # ------------------------------------------------------------
  # convenience contrasts from standardized prognostic scores
  # ------------------------------------------------------------
  if (isTRUE(doublepg)) {
    df$pg1ct <- ifelse(df[[A1.var]] == 1, df$pg11, df$pg01)
    df$pg1tc <- ifelse(df[[A1.var]] == 1, df$pg01, df$pg11)

    df$pg2ct <- ifelse(
      df[[eta2.var]] == 1,
      ifelse(df[[A2.var]] == 1, df$pg12, df$pg02),
      NA_real_
    )
    df$pg2tc <- ifelse(
      df[[eta2.var]] == 1,
      ifelse(df[[A2.var]] == 1, df$pg02, df$pg12),
      NA_real_
    )
  }

  df
}



#' Plot propensity (or censoring) score overlap by group
#'
#' @description
#' Produces a diagnostic overlap plot for a set of propensity-like scores (e.g., treatment propensity
#' scores or censoring propensities) stratified by a binary group indicator \code{A}. The function
#' visualizes the empirical score distributions using semi-transparent histograms on a density scale.
#'
#' This plot is primarily intended to assess common support / overlap and to diagnose separation or
#' extreme predicted probabilities before matching/weighting steps.
#'
#' @details
#' The input \code{A} is coerced to a factor and used for coloring and filling the histogram.
#' The y-axis is scaled to density (\code{..density..}). The function returns a \pkg{ggplot2} object
#' (invisibly) and also prints the plot as a side effect, which is convenient in interactive use.
#'
#' \strong{Package note.} In a package context, avoid \code{library(ggplot2)} inside functions.
#' Instead, use \code{ggplot2::} calls and guard availability via \code{requireNamespace("ggplot2", quietly = TRUE)}.
#'
#' @param ps Numeric vector of propensity-like scores, typically in \eqn{[0,1]}.
#' @param A Vector defining the grouping variable (e.g., treatment arm or event indicator). Will be
#' coerced to a factor for plotting.
#'
#' @return A \pkg{ggplot2} plot object (class \code{"gg"} and \code{"ggplot"}). The plot is also printed.
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_histogram}}
#' @export
propensityplot <- function(ps, A) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for propensityplot().", call. = FALSE)
  }

  ps_dfctrl <- data.frame(
    A  = factor(A),
    ps = ps
  )

  p <- ggplot2::ggplot(ps_dfctrl, ggplot2::aes(x = ps, fill = A, color = A)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 30,
      alpha = 0.3,
      position = "identity"
    ) +
    ggplot2::labs(
      title = "Density Plot of Propensity Scores",
      x = "Propensity Score",
      y = "Density",
      fill = "Group (A)",
      color = "Group (A)"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}

