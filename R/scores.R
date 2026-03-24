#' Compute propensity and prognostic scores (treatment- or censoring-based)
#'
#' @description
#' Computes (i) a propensity score and (ii) prognostic score(s) for use in
#' matching/weighting pipelines. Supports SuperLearner-based fitting or
#' glm/glmnet-based alternatives. Can compute treatment PS + treatment-specific
#' prognostic scores (pg0/pg1), or censoring PS/PG (pscens/pgcens) when
#' \code{censmod=TRUE}.
#'
#' @param data data.frame with all variables.
#' @param id Character scalar. Subject ID column name.
#' @param Y Character scalar. Survival time column name.
#' @param event Character scalar. Event indicator column name (1=event, 0=censored).
#' @param X Character vector. Covariate column names for prognostic model(s).
#' @param A Character scalar. Treatment indicator column name (coded 0/1 or -1/1).
#' @param Xtrt Optional character vector. Covariates for treatment model if different from \code{X}.
#' @param doublepg Logical. If TRUE and \code{censmod=FALSE}, fit separate prognostic models by treatment arm.
#' @param outer_CV Integer. Outer CV folds.
#' @param inner_CV Optional integer. Inner CV folds for nested CV (SuperLearner).
#' @param stratifyCV Logical. Whether to stratify CV folds.
#' @param cores Integer. Requested cores for parallel fit (multicore only on non-Windows).
#' @param tau Optional numeric. Truncation horizon for mean survival.
#' @param sl.seed Integer. Seed for SuperLearner.
#' @param A.SL.library Character vector. SL learners for propensity/censoring models.
#' @param Y.SL.library Character vector. Learners for survivalSL prognostic models.
#' @param A.method Character. CV risk method for propensity SL.
#' @param Y.method Character. Metric for survivalSL.
#' @param param.tune,ngrid,param.weights.fix,param.weights.init,optim.method,penalty,maxit
#' Control survivalSL fitting and prediction.
#' @param pgcens,pscens,censmod Logical flags controlling censoring scores.
#' @param model.pg Character. "cox" or "aft" for non-SL prognostic modeling.
#' @param standardize Logical. Standardize covariates for glmnet.
#' @param superLearn Logical. Use SuperLearner/survivalSL branches if TRUE.
#' @param pslink Character. "logit" or "probit".
#' @param pglink Character. AFT distribution for flexsurvreg.
#'
#' @return A data.frame with stable columns:
#' \itemize{
#'   \item \code{id}: ID values
#'   \item \code{ps}: treatment propensity score (if computed)
#'   \item \code{pg0}, \code{pg1}: treatment-specific prognostic scores (if computed)
#'   \item \code{pscens}: censoring propensity score (if computed)
#'   \item \code{pgcens}: censoring prognostic score (if computed)
#' }
#'
#' @export
ComputeScores <- function(data, id, Y, event, X, A,
                          Xtrt = NULL,
                          doublepg = TRUE,
                          outer_CV = 5,
                          inner_CV = NULL,
                          stratifyCV = FALSE,
                          cores = 5, tau = NULL,
                          sl.seed = 100,
                          A.SL.library = c("SL.mean","SL.glm","SL.glmnet","SL.ranger","SL.xgboost"),
                          Y.SL.library = c("LIB_COXen","LIB_AFTggamma"),
                          A.method = "method.AUC", Y.method = "auc",
                          param.tune = NULL, ngrid = 2000,
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
                          pglink = "lognormal") {

  # ---------- input validation (fast fail) ----------
  stopifnot(is.data.frame(data))
  for (nm in c(id, Y, event, A)) {
    if (!nm %in% names(data)) stop("Column '", nm, "' not found in data.", call. = FALSE)
  }
  if (!is.character(X) || length(X) < 1L) stop("X must be a non-empty character vector.", call. = FALSE)
  if (!all(X %in% names(data))) stop("Some X columns not found in data.", call. = FALSE)
  if (!is.null(Xtrt) && !all(Xtrt %in% names(data))) stop("Some Xtrt columns not found in data.", call. = FALSE)
  outer_CV <- as.integer(outer_CV)
  if (outer_CV < 2L) stop("outer_CV must be >= 2.", call. = FALSE)
  cores <- as.integer(cores)
  if (cores < 1L) stop("cores must be >= 1.", call. = FALSE)

  # Event indicator must be 0/1 (or logical)
  Event <- data[[event]]
  if (is.logical(Event)) Event <- as.integer(Event)
  if (!all(Event %in% c(0L, 1L, NA_integer_))) {
    stop("'", event, "' must be coded 0/1 (1=event, 0=censored).", call. = FALSE)
  }

  # Treatment indicator: accept -1/1 or 0/1
  AA <- data[[A]]
  A_bin <- ifelse(AA == 1, 1L, 0L)
  if (!all(A_bin %in% c(0L, 1L, NA_integer_))) stop("Treatment variable must be binary or -1/1.", call. = FALSE)

  Id <- data[[id]]
  YY <- data[[Y]]
  XX <- data[, X, drop = FALSE]
  x  <- data.matrix(XX)

  xtrt <- if (!is.null(Xtrt)) data.matrix(data[, Xtrt, drop = FALSE]) else x

  loc1 <- which(A_bin == 1L)
  loc0 <- which(A_bin == 0L)

  # stable return contract
  out <- data.frame(
    id     = Id,
    ps     = NA_real_,
    pg0    = NA_real_,
    pg1    = NA_real_,
    pscens = NA_real_,
    pgcens = NA_real_
  )
  names(out)[1] <- id

  # ---------- helpers ----------
  pred_mean <- function(Smat, time_grid) {
    if (length(time_grid) != ncol(Smat)) stop("time_grid length mismatch with Smat.", call. = FALSE)
    if (any(diff(time_grid) <= 0)) stop("time_grid must be strictly increasing.", call. = FALSE)
    apply(Smat, 1, function(surv_i) {
      dt   <- diff(time_grid)
      mids <- (surv_i[-1] + surv_i[-length(surv_i)]) / 2
      sum(mids * dt)
    })
  }

  parallel_mode <- if (cores > 1L && .Platform$OS.type != "windows") "multicore" else "seq"

  # dependency guards
  if (superLearn) {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required when superLearn=TRUE.", call. = FALSE)
    }
  }
  # survival is effectively required for Surv/coxph branches
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.", call. = FALSE)
  }
  if (!superLearn && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required when superLearn=FALSE.", call. = FALSE)
  }

  if (!superLearn && !requireNamespace("survivalSL", quietly = TRUE)) {
    stop("Package 'survivalSL' is required when superLearn=TRUE", call. = FALSE)
  }
  if (!superLearn && model.pg == "aft" && !requireNamespace("flexsurv", quietly = TRUE)) {
    stop("Package 'flexsurv' is required for model.pg='aft'.", call. = FALSE)
  }

  # ---------- main branches ----------
  if (censmod && superLearn) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic("Scores: censoring scores (SL)")

    # censoring PS: P(Event=1 | X) if that's what you intend
    if (pscens) {
      innerCvControl_value <- if (!is.null(inner_CV)) {
        rep(list(list(V = inner_CV, stratifyCV = stratifyCV)), outer_CV)
      } else NULL

      set.seed(sl.seed, "L'Ecuyer-CMRG")
      sl_out <- SuperLearner::CV.SuperLearner(
        Y = Event,
        X = as.data.frame(x),
        family = stats::binomial(link = pslink),
        method = A.method,
        SL.library = A.SL.library,
        cvControl = list(V = outer_CV, stratifyCV = stratifyCV),
        innerCvControl = innerCvControl_value,
        parallel = parallel_mode
      )
      out[["pscens"]] <- as.numeric(sl_out$SL.predict)
    }

    if (pgcens) {
      # simple PG on uncensored only; you may want to revisit this estimand
      uncensored <- data[Event == 1L, c(Y, X, A), drop = FALSE]
      if (nrow(uncensored) > 0) {
        set.seed(sl.seed, "L'Ecuyer-CMRG")
        sl_fit <- SuperLearner::SuperLearner(
          Y = uncensored[[Y]],
          X = uncensored[, c(X, A), drop = FALSE],
          family = stats::gaussian(),
          SL.library = A.SL.library,
          cvControl = list(V = outer_CV, stratifyCV = stratifyCV)
        )
        out[["pgcens"]] <- as.numeric(SuperLearner:::predict.SuperLearner(sl_fit, newdata = data[, c(X, A), drop = FALSE])$pred)
      }
    }

    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()

  } else if (!censmod && doublepg && superLearn) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic("Scores: treatment PS + pg0/pg1 (SL)")

    # treatment PS
    innerCvControl_value <- if (!is.null(inner_CV)) {
      rep(list(list(V = inner_CV, stratifyCV = stratifyCV)), outer_CV)
    } else NULL

    set.seed(sl.seed, "L'Ecuyer-CMRG")
    sl_out <- SuperLearner::CV.SuperLearner(
      Y = A_bin,
      X = as.data.frame(xtrt),
      family = stats::binomial(link = pslink),
      method = A.method,
      SL.library = A.SL.library,
      cvControl = list(V = outer_CV, stratifyCV = stratifyCV),
      innerCvControl = innerCvControl_value,
      parallel = parallel_mode
    )
    out[["ps"]] <- as.numeric(sl_out$SL.predict)

    # prognostic via survivalSL (NOTE: survivalSL must exist in your package or dependency)
    # if (!exists("survivalSL", mode = "function")) {
    #   stop("Function 'survivalSL' not found. Include it in your package or import from its source.", call. = FALSE)
    # }

    X_df <- data[, X, drop = FALSE]
    tmax <- if (is.null(tau)) max(YY, na.rm = TRUE) else tau

    if (length(loc1) > 0) {
      survdata1 <- data[loc1, c(Y, event, X), drop = FALSE]
      f1 <- stats::as.formula(paste0("survival::Surv(", Y, ",", event, ") ~ ", paste(X, collapse = "+")))
      slres1 <- survivalSL::survivalSL(
        formula = f1, methods = Y.SL.library, metric = Y.method, data = survdata1,
        cv = outer_CV, param.tune = param.tune, seed = sl.seed,
        param.weights.fix = param.weights.fix, param.weights.init = param.weights.init,
        maxit = maxit, penalty = penalty, show_progress = TRUE
      )
      grid1 <- seq(0, tmax, length.out = ngrid)
      pred1 <- survivalSL:::predict.sltime(slres1, newdata = X_df, newtimes = grid1)
      out[["pg1"]] <- pred_mean(pred1$predictions$sl, pred1$times)
    }

    if (length(loc0) > 0) {
      survdata0 <- data[loc0, c(Y, event, X), drop = FALSE]
      f0 <- stats::as.formula(paste0("survival::Surv(", Y, ",", event, ") ~ ", paste(X, collapse = "+")))
      slres0 <- survivalSL::survivalSL(
        formula = f0, methods = Y.SL.library, metric = Y.method, data = survdata0,
        cv = outer_CV, param.tune = param.tune, seed = sl.seed,
        param.weights.fix = param.weights.fix, param.weights.init = param.weights.init,
        maxit = maxit, penalty = penalty, show_progress = TRUE
      )
      grid0 <- seq(0, tmax, length.out = ngrid)
      pred0 <- survivalSL:::predict.sltime(slres0, newdata = X_df, newtimes = grid0)
      out[["pg0"]] <- pred_mean(pred0$predictions$sl, pred0$times)
    }

    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()

  } else if (censmod && !superLearn) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic("Scores: censoring scores (glm/glmnet)")

    # censoring PS
    if (pscens) {
      if (ncol(x) > 1L) {
        cvfit <- glmnet::cv.glmnet(x, Event, family = stats::binomial(link = pslink),
                                   nfolds = outer_CV, standardize = standardize)
        out[["pscens"]] <- as.numeric(predict(cvfit, newx = x, s = "lambda.min", type = "response"))
      } else {
        df1 <- data.frame(Event = Event, x = x[, 1])
        fit <- stats::glm(Event ~ x, data = df1, family = stats::binomial(link = pslink))
        out[["pscens"]] <- as.numeric(stats::predict(fit, newdata = df1, type = "response"))
      }
    }

    # censoring PG (you used gaussian on Event previously; that seems off; keep as-is but consistent)
    if (pgcens) {
      Xpg <- data.matrix(data[, c(X, A), drop = FALSE])
      cvfit <- glmnet::cv.glmnet(Xpg, Event, family = "gaussian",
                                 nfolds = outer_CV, standardize = standardize)
      out[["pgcens"]] <- as.numeric(predict(cvfit, newx = Xpg, s = "lambda.min", type = "response"))
    }

    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()

  } else if (!censmod && doublepg && !superLearn) {
    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic("Scores: treatment PS + pg0/pg1 (glm/glmnet)")

    # treatment PS
    if (ncol(xtrt) > 1L) {
      cvfit <- glmnet::cv.glmnet(xtrt, A_bin, family = stats::binomial(link = pslink),
                                 nfolds = outer_CV, standardize = standardize)
      out[["ps"]] <- as.numeric(predict(cvfit, newx = xtrt, s = "lambda.min", type = "response"))
    } else {
      df1 <- data.frame(A_bin = A_bin, xtrt = xtrt[, 1])
      fit <- stats::glm(A_bin ~ xtrt, data = df1, family = stats::binomial(link = pslink))
      out[["ps"]] <- as.numeric(stats::predict(fit, newdata = df1, type = "response"))
    }

    # prognostic models
    if (model.pg == "cox") {
      if (length(loc1) > 0) {
        Y1 <- survival::Surv(YY[loc1], Event[loc1])
        X1 <- data.matrix(XX[loc1, , drop = FALSE])
        if (ncol(X1) > 1L) {
          cv1 <- glmnet::cv.glmnet(X1, Y1, family = "cox", nfolds = outer_CV, standardize = standardize)
          out[["pg1"]] <- as.numeric(predict(cv1, newx = x, s = "lambda.min", type = "link"))
        } else {
          dfc <- data.frame(time = YY[loc1], status = Event[loc1], x1 = X1[, 1])
          fit <- survival::coxph(survival::Surv(time, status) ~ x1, data = dfc)
          out[["pg1"]] <- as.numeric(stats::predict(fit, newdata = data.frame(x1 = x[, 1]), type = "lp"))
        }
      }
      if (length(loc0) > 0) {
        Y0 <- survival::Surv(YY[loc0], Event[loc0])
        X0 <- data.matrix(XX[loc0, , drop = FALSE])
        if (ncol(X0) > 1L) {
          cv0 <- glmnet::cv.glmnet(X0, Y0, family = "cox", nfolds = outer_CV, standardize = standardize)
          out[["pg0"]] <- as.numeric(predict(cv0, newx = x, s = "lambda.min", type = "link"))
        } else {
          dfc <- data.frame(time = YY[loc0], status = Event[loc0], x0 = X0[, 1])
          fit <- survival::coxph(survival::Surv(time, status) ~ x0, data = dfc)
          out[["pg0"]] <- as.numeric(stats::predict(fit, newdata = data.frame(x0 = x[, 1]), type = "lp"))
        }
      }

    } else if (model.pg == "aft") {
      # AFT pg as mean survival time from flexsurvreg
      if (length(loc1) > 0) {
        d1 <- data.frame(time = YY[loc1], status = Event[loc1], XX[loc1, , drop = FALSE])
        fit1 <- flexsurv::flexsurvreg(survival::Surv(time, status) ~ ., data = d1, dist = pglink)
        out[["pg1"]] <- as.numeric(predict(fit1, newdata = as.data.frame(XX), type = "mean"))
      }
      if (length(loc0) > 0) {
        d0 <- data.frame(time = YY[loc0], status = Event[loc0], XX[loc0, , drop = FALSE])
        fit0 <- flexsurv::flexsurvreg(survival::Surv(time, status) ~ ., data = d0, dist = pglink)
        out[["pg0"]] <- as.numeric(predict(fit0, newdata = as.data.frame(XX), type = "mean"))
      }
    }

    if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()
  } else {
    stop("Unsupported combination of flags (censmod/doublepg/superLearn).", call. = FALSE)
  }

  out
}




#' Compute stage-specific treatment propensity and prognostic “double scores”
#'
#' @description
#' Computes and attaches stage-specific \emph{double scores} for a two-stage treatment setting.
#' The function is a thin orchestrator around \code{\link{ComputeScores}} that:
#' (i) restricts to stage-2 entrants (\code{eta2==1}) to compute stage-2 scores using \code{A2} and
#' \code{Y2}, then (ii) computes stage-1 scores on the full cohort using \code{A1} and either \code{OY}
#' (overall outcome) or \code{Y1} (stage-1 time) depending on \code{adjustdelta1}.
#'
#' The output is the original dataset augmented with both \emph{raw} score columns (propensities and
#' prognostic scores) and \emph{standardized} score columns intended for distance-based matching or
#' downstream modeling.
#'
#' @details
#' \strong{What is computed.}
#' \itemize{
#'   \item A treatment propensity score \code{ps = P(A=1|Xtrt)} for each stage.
#'   \item A prognostic score for each stage based on survival modeling.
#' }
#'
#' If \code{doublepg=TRUE}, prognostic scores are computed separately under each treatment level:
#' \code{pg0} (under \code{A=0}) and \code{pg1} (under \code{A=1}). If \code{doublepg=FALSE}, a single
#' prognostic score \code{pg} is computed.
#'
#' \strong{Stage 2.} Subjects are subset to \code{eta2==1} and \code{ComputeScores()} is called with
#' \code{Y=Y2.var}, \code{A=A2.var}, covariates \code{names.var2} (prognostic model) and \code{Xtrt2}
#' (treatment model). Stage-2 results are merged back to the full dataset; non-entrants receive \code{NA}
#' for stage-2 scores.
#'
#' \strong{Stage 1.} \code{ComputeScores()} is called on the full cohort with \code{A=A1.var} and either:
#' \itemize{
#'   \item \code{Y=OY.var, event=delta.var} if \code{adjustdelta1=FALSE}, or
#'   \item \code{Y=Y1.var, event='deltaadj'} if \code{adjustdelta1=TRUE}.
#' }
#' When \code{adjustdelta1=TRUE}, \code{deltaadj} is created by copying \code{delta.var} and setting
#' \code{deltaadj=0} for \code{eta2==1} and \code{delta==1}.
#'
#' \strong{Transformations and standardization.} For each stage, \code{ps} is transformed using
#' \code{qlogis(ps)} (logit scale) and then all score columns are z-scored using \code{scale()} to produce
#' standardized columns (e.g., \code{ps1}, \code{pg01}, \code{pg11}).
#'
#' \strong{Convenience contrasts.} When \code{doublepg=TRUE}, the function creates:
#' \itemize{
#'   \item \code{pg1ct} / \code{pg1tc}: stage-1 “correct-treatment” and “treatment-contrast” prognostic scores
#'   \item \code{pg2ct} / \code{pg2tc}: analogous stage-2 versions (for \code{eta2==1})
#' }
#' where “correct-treatment” selects \code{pg0} if observed \code{A=0} and \code{pg1} if observed \code{A=1},
#' and “treatment-contrast” selects the opposite arm’s prognostic score.
#'
#' \strong{Switch behavior.} If \code{useds=FALSE}, the function returns \code{data} unchanged (no-op),
#' which is useful in pipelines where score construction is optional.
#'
#' @param data A data.frame containing all required stage-1 and stage-2 variables.
#' @param id.var Character scalar. Subject identifier column name.
#' @param eta2.var Character scalar. Stage-2 entry indicator column name (1=entered stage 2, 0=did not).
#'
#' @param Y1.var Character scalar. Stage-1 time/outcome component (used only when \code{adjustdelta1=TRUE}).
#' @param Y2.var Character scalar. Stage-2 outcome/time column used for stage-2 prognostic scoring.
#' @param delta.var Character scalar. Event indicator column name used for survival modeling.
#' @param OY.var Character scalar. Overall outcome/time column used for stage-1 prognostic scoring when
#' \code{adjustdelta1=FALSE}.
#'
#' @param A1.var Character scalar. Stage-1 treatment indicator column name.
#' @param A2.var Character scalar. Stage-2 treatment indicator column name.
#'
#' @param names.var1 Character vector. Covariate names for the stage-1 prognostic model.
#' @param names.var2 Character vector. Covariate names for the stage-2 prognostic model (stage-2 entrants only).
#'
#' @param Xtrt1 Character vector. Covariate names for the stage-1 treatment propensity model (if different from
#' \code{names.var1}).
#' @param Xtrt2 Character vector. Covariate names for the stage-2 treatment propensity model (if different from
#' \code{names.var2}).
#'
#' @param useds Logical. If TRUE, compute and merge scores. If FALSE, return \code{data} unchanged.
#'
#' @param cores Integer. Number of cores passed to \code{ComputeScores} (if supported by the backend).
#' @param tau Optional numeric. Truncation horizon used in prognostic mean calculations inside \code{ComputeScores}.
#' @param sl.seed Integer. RNG seed passed to \code{ComputeScores}.
#'
#' @param A.SL.library1,A.SL.library2 Character vectors. SuperLearner libraries for stage-1 and stage-2 treatment models.
#' @param Y.SL.library Character vector. Learners for survivalSL prognostic modeling.
#' @param A.method,Y.method Optional. Scoring metrics passed to \code{ComputeScores}.
#'
#' @param param.weights.fix,param.weights.init,optim.method,maxit,penalty1,penalty2,param.tune
#' Tuning/optimization controls forwarded to \code{ComputeScores}.
#'
#' @param ngrid Integer. Number of grid points used when integrating survival curves for mean survival time.
#'
#' @param censmod Logical. Included for interface consistency; in this wrapper the calls to \code{ComputeScores}
#' set \code{censmod=FALSE} to compute treatment/prognostic (not censoring) scores.
#'
#' @param doublepg Logical. If TRUE, compute \code{pg0} and \code{pg1}. If FALSE, compute a single \code{pg}.
#'
#' @param adjustdelta1 Logical. If TRUE, define \code{deltaadj} and use \code{Y1.var} as the time variable
#' for stage-1 scoring; otherwise use \code{OY.var} and \code{delta.var}.
#'
#' @param plotps Logical. If TRUE, plots propensity distributions by treatment at each stage using \code{propensityplot()}.
#'
#' @param model.pg Character. Prognostic model family used when \code{superLearn=FALSE} ("cox" or "aft").
#' @param standardize Logical. Whether to standardize covariates for glmnet when \code{superLearn=FALSE}.
#' @param superLearn Logical. If TRUE, use SuperLearner-based estimation inside \code{ComputeScores}; otherwise use glm/glmnet.
#' @param pslink Character. Link for binomial treatment propensity model ("logit" or "probit").
#' @param pglink Character. AFT distribution used when \code{model.pg="aft"} (passed to \code{ComputeScores}).
#'
#' @return A data.frame equal to \code{data} augmented with score columns. If \code{useds=FALSE},
#' returns \code{data} unchanged.
#'
#' \strong{Raw score columns} (merged back by \code{id.var}):
#' \itemize{
#'   \item Stage 1: \code{prog01, prog11, prop1} (or \code{prog01, prop1} if \code{doublepg=FALSE})
#'   \item Stage 2: \code{prog02, prog12, prop2} (or \code{prog02, prop2} if \code{doublepg=FALSE}; \code{NA} for \code{eta2==0})
#' }
#'
#' \strong{Standardized columns} (z-scored; propensity on logit scale):
#' \itemize{
#'   \item Stage 1: \code{pg01, pg11, ps1} (or \code{pg01, ps1})
#'   \item Stage 2: \code{pg02, pg12, ps2} (or \code{pg02, ps2})
#' }
#'
#' When \code{doublepg=TRUE}, additional convenience columns are created:
#' \itemize{
#'   \item Stage 1: \code{pg1ct}, \code{pg1tc}
#'   \item Stage 2: \code{pg2ct}, \code{pg2tc}
#' }
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
    Xtrt1,
    Xtrt2,
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
    maxit         = 1000,
    penalty1      = NULL,
    penalty2      = NULL,
    ngrid         = 50,
    censmod       = TRUE,
    doublepg      = TRUE,
    param.tune    = NULL,
    adjustdelta1  = FALSE,
    plotps        = FALSE,
    model.pg      = "cox",      # "cox" or "aft"
    standardize   = FALSE,      # Standardize covariates for glmnet
    superLearn    = TRUE,       # Whether to use SuperLearner or glmnet
    pslink        = "logit",    # "logit" or "probit"
    pglink        = NULL        # e.g., "lognormal" when model.pg == "aft"
) {
  stopifnot(is.data.frame(data))

  # basic column checks
  req_cols <- c(id.var, eta2.var, Y1.var, Y2.var, delta.var, OY.var, A1.var, A2.var)
  miss <- setdiff(req_cols, names(data))
  if (length(miss) > 0L) {
    stop("Missing required columns in data: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!all(names.var1 %in% names(data))) stop("Some names.var1 not found in data.", call. = FALSE)
  if (!all(names.var2 %in% names(data))) stop("Some names.var2 not found in data.", call. = FALSE)
  if (!all(Xtrt1 %in% names(data)))      stop("Some Xtrt1 not found in data.", call. = FALSE)
  if (!all(Xtrt2 %in% names(data)))      stop("Some Xtrt2 not found in data.", call. = FALSE)

  if (!useds) return(data)

  if (!exists("ComputeScores", mode = "function")) {
    stop("ComputeScores() not found. It must be available in your package/environment.", call. = FALSE)
  }

  # optional timing (no hard dependency)
  tick <- function(...) if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic(...)
  tock <- function(...) if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()

  # safe logit to avoid +/-Inf when ps=0/1
  safe_qlogis <- function(p, eps = 1e-6) {
    p <- as.numeric(p)
    p <- pmax(pmin(p, 1 - eps), eps)
    stats::qlogis(p)
  }

  df <- data

  # optional delta adjustment for stage-1 scoring
  if (adjustdelta1) {
    df$deltaadj <- df[[delta.var]]
    df$deltaadj[df[[eta2.var]] == 1 & df[[delta.var]] == 1] <- 0
  }

  # -------------------------
  # Stage 2: stage-2 entrants
  # -------------------------
  df2 <- df[df[[eta2.var]] == 1, , drop = FALSE]

  if (nrow(df2) > 0L) {
    tick("Stage 2 DoubleScore")

    ds2 <- ComputeScores(
      data         = df2,
      id           = id.var,
      Y            = Y2.var,
      event        = delta.var,
      X            = names.var2,
      Xtrt         = Xtrt2,
      A            = A2.var,
      doublepg     = doublepg,
      outer_CV     = 5,
      inner_CV     = 5,
      stratifyCV   = FALSE,
      cores        = cores,
      tau          = tau,
      sl.seed      = sl.seed,
      A.SL.library = A.SL.library2,
      Y.SL.library = Y.SL.library,
      A.method     = A.method,
      Y.method     = Y.method,
      param.weights.fix  = param.weights.fix,
      param.weights.init = param.weights.init,
      optim.method = optim.method,
      maxit        = maxit,
      penalty      = penalty2,
      ngrid        = ngrid,
      pscens       = FALSE,
      pgcens       = FALSE,
      censmod      = FALSE,      # treatment score + prognostic score (not censoring)
      param.tune   = param.tune,
      model.pg     = model.pg,
      standardize  = standardize,
      superLearn   = superLearn,
      pslink       = pslink,
      pglink       = pglink
    )

    tock()

    ds2 <- as.data.frame(ds2)
    if (!id.var %in% names(ds2)) names(ds2)[1] <- id.var

    if (doublepg) {
      # expected columns: id, pg0, pg1, ps
      if (!all(c("pg0", "pg1", "ps") %in% names(ds2))) {
        stop("Stage-2 ComputeScores output must contain pg0, pg1, ps when doublepg=TRUE.", call. = FALSE)
      }
      ds2$pg0 <- as.numeric(ds2$pg0)
      ds2$pg1 <- as.numeric(ds2$pg1)
      ds2$ps  <- as.numeric(ds2$ps)

      if (plotps && exists("propensityplot", mode = "function")) {
        propensityplot(ps = ds2[["ps"]], A = df2[[A2.var]])
      }

      # logit transform + standardize for matching
      ds2$ps_logit <- safe_qlogis(ds2$ps)
      dsp2 <- as.data.frame(scale(ds2[, c("pg0", "pg1", "ps_logit"), drop = FALSE]))
      colnames(dsp2) <- c("pg02", "pg12", "ps2")
      dsp2[[id.var]] <- ds2[[id.var]]

      # keep raw (renamed) columns for interpretability
      ds2_out <- ds2[, c(id.var, "pg0", "pg1", "ps"), drop = FALSE]
      colnames(ds2_out) <- c(id.var, "prog02", "prog12", "prop2")

      # merge
      df <- merge(df, ds2_out, by = id.var, all.x = TRUE)
      df <- merge(df, dsp2,   by = id.var, all.x = TRUE)

      # convenience: correct-treatment vs treatment-contrast prognostic score
      df$pg2ct <- ifelse(df[[A2.var]] == 1, df$pg12, df$pg02)
      df$pg2tc <- ifelse(df[[A2.var]] == 1, df$pg02, df$pg12)

    } else {
      # expected columns: id, pg, ps
      if (!all(c("pg", "ps") %in% names(ds2))) {
        stop("Stage-2 ComputeScores output must contain pg, ps when doublepg=FALSE.", call. = FALSE)
      }
      ds2$pg <- as.numeric(ds2$pg)
      ds2$ps <- as.numeric(ds2$ps)

      if (plotps && exists("propensityplot", mode = "function")) {
        propensityplot(ps = ds2[["ps"]], A = df2[[A2.var]])
      }

      ds2$ps_logit <- safe_qlogis(ds2$ps)
      dsp2 <- as.data.frame(scale(ds2[, c("pg", "ps_logit"), drop = FALSE]))
      colnames(dsp2) <- c("pg02", "ps2")
      dsp2[[id.var]] <- ds2[[id.var]]

      ds2_out <- ds2[, c(id.var, "pg", "ps"), drop = FALSE]
      colnames(ds2_out) <- c(id.var, "prog02", "prop2")

      df <- merge(df, ds2_out, by = id.var, all.x = TRUE)
      df <- merge(df, dsp2,   by = id.var, all.x = TRUE)
    }
  } else {
    # no stage-2 entrants; still create empty columns? (leave as-is)
    message("No stage-2 entrants (eta2==1); skipping stage-2 score computation.")
  }

  # -------------------------
  # Stage 1: full cohort
  # -------------------------
  tick("Stage 1 DoubleScore")

  ds1 <- ComputeScores(
    data         = df,
    id           = id.var,
    Y            = if (adjustdelta1) Y1.var else OY.var,
    event        = if (adjustdelta1) "deltaadj" else delta.var,
    X            = names.var1,
    Xtrt         = Xtrt1,
    A            = A1.var,
    doublepg     = doublepg,
    outer_CV     = 5,
    inner_CV     = 5,
    stratifyCV   = FALSE,
    cores        = cores,
    tau          = tau,
    sl.seed      = sl.seed,
    A.SL.library = A.SL.library1,
    Y.SL.library = Y.SL.library,
    A.method     = A.method,
    Y.method     = Y.method,
    param.weights.fix  = param.weights.fix,
    param.weights.init = param.weights.init,
    optim.method = optim.method,
    maxit        = maxit,
    penalty      = penalty1,
    ngrid        = ngrid,
    pscens       = FALSE,
    pgcens       = FALSE,
    censmod      = FALSE,
    param.tune   = param.tune,
    model.pg     = model.pg,
    standardize  = standardize,
    superLearn   = superLearn,
    pslink       = pslink,
    pglink       = pglink
  )

  tock()

  ds1 <- as.data.frame(ds1)
  if (!id.var %in% names(ds1)) names(ds1)[1] <- id.var

  if (doublepg) {
    if (!all(c("pg0", "pg1", "ps") %in% names(ds1))) {
      stop("Stage-1 ComputeScores output must contain pg0, pg1, ps when doublepg=TRUE.", call. = FALSE)
    }
    ds1$pg0 <- as.numeric(ds1$pg0)
    ds1$pg1 <- as.numeric(ds1$pg1)
    ds1$ps  <- as.numeric(ds1$ps)

    if (plotps && exists("propensityplot", mode = "function")) {
      propensityplot(ps = ds1[["ps"]], A = df[[A1.var]])
    }

    ds1$ps_logit <- safe_qlogis(ds1$ps)
    dsp1 <- as.data.frame(scale(ds1[, c("pg0", "pg1", "ps_logit"), drop = FALSE]))
    colnames(dsp1) <- c("pg01", "pg11", "ps1")
    dsp1[[id.var]] <- ds1[[id.var]]

    ds1_out <- ds1[, c(id.var, "pg0", "pg1", "ps"), drop = FALSE]
    colnames(ds1_out) <- c(id.var, "prog01", "prog11", "prop1")

    df <- merge(df, ds1_out, by = id.var, all.x = TRUE)
    df <- merge(df, dsp1,   by = id.var, all.x = TRUE)

    df$pg1ct <- ifelse(df[[A1.var]] == 1, df$pg11, df$pg01)
    df$pg1tc <- ifelse(df[[A1.var]] == 1, df$pg01, df$pg11)

  } else {
    if (!all(c("pg", "ps") %in% names(ds1))) {
      stop("Stage-1 ComputeScores output must contain pg, ps when doublepg=FALSE.", call. = FALSE)
    }
    ds1$pg <- as.numeric(ds1$pg)
    ds1$ps <- as.numeric(ds1$ps)

    if (plotps && exists("propensityplot", mode = "function")) {
      propensityplot(ps = ds1[["ps"]], A = df[[A1.var]])
    }

    ds1$ps_logit <- safe_qlogis(ds1$ps)
    dsp1 <- as.data.frame(scale(ds1[, c("pg", "ps_logit"), drop = FALSE]))
    colnames(dsp1) <- c("pg01", "ps1")
    dsp1[[id.var]] <- ds1[[id.var]]

    ds1_out <- ds1[, c(id.var, "pg", "ps"), drop = FALSE]
    colnames(ds1_out) <- c(id.var, "prog01", "prop1")

    df <- merge(df, ds1_out, by = id.var, all.x = TRUE)
    df <- merge(df, dsp1,   by = id.var, all.x = TRUE)
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

