#' Impute censored stage-specific outcomes via matching with optional learned censoring scores
#'
#' @description
#' Performs two-stage imputation of censored outcomes in a sequential treatment / DTR setting.
#' The function (i) constructs a stage-2 composite outcome among subjects who enter stage 2
#' (\code{eta2==1}) by matching censored subjects (\code{delta2==0}) to eligible uncensored donors,
#' and (ii) constructs a stage-1 composite outcome for all subjects by matching stage-1 censored
#' subjects (\code{delta1==0}) to eligible donors, optionally leveraging learned censoring propensity
#' and prognostic scores in the matching distance.
#'
#' Matching is performed by helper routines \code{impute_censored_stage2()} and
#' \code{impute_censored_stage1()}, which implement nearest-neighbor or optimal matching using
#' Mahalanobis or other distances, with optional exact matching constraints. When \code{usecov=FALSE},
#' matching covariates are replaced by learned censoring scores computed by \code{ComputeScores()}:
#' \code{pscens} (a censoring propensity score) and/or \code{pgcens} (a censoring prognostic score).
#'
#' @details
#' \strong{Stage-2 imputation.} The function subsets to \code{eta2==1} and imputes \code{Y2} for
#' censored subjects using donors with observed stage-2 outcomes, subject to donor eligibility
#' constraints encoded in \code{impute_censored_stage2()} (e.g., donor observed time exceeding the
#' recipient’s observed time). The result is stored as \code{compY2}. For uncensored subjects,
#' \code{compY2} is set to the observed \code{Y2}.
#'
#' \strong{Stage-1 composite outcome.} A composite outcome \code{compOY} is constructed as
#' \code{compOY = Y1 + compY2} for stage-2 entrants; for non-entrants (\code{eta2==0}) with observed
#' final outcome, \code{compOY} is set to the observed overall outcome \code{OY}. Stage-1 censoring
#' indicators \code{delta1} and \code{delta2} are defined from whether \code{compOY} is observed.
#' If \code{adjustdelta1=TRUE}, the stage-1 event indicator is modified to treat certain stage-2
#' events as censored (see the \code{adjustdelta1} section below).
#'
#' \strong{Censoring-score matching (optional).} If \code{usecov=FALSE}, censoring scores are computed
#' via \code{ComputeScores(censmod=TRUE, doublepg=FALSE)} separately for stage 2 and stage 1.
#' The raw probability scores are optionally logit-transformed and standardized (z-scored) prior
#' to matching. These scores are then used as matching covariates in place of the original covariates.
#'
#' \strong{Important assumption.} This function assumes \code{delta==1} denotes an observed event/outcome
#' (i.e., not censored) and \code{delta==0} denotes censoring, consistent with \code{Surv(time, delta)}.
#' If your data use the opposite convention, you must recode before calling this function.
#'
#' @param data A data.frame containing all variables needed for stage-1 and stage-2 processing.
#'
#' @param id.var Character scalar. Subject identifier column name.
#' @param eta2.var Character scalar. Stage-2 entry indicator column name (1=entered stage 2, 0=did not).
#'
#' @param Y1.var Character scalar. Stage-1 time/outcome component used in the composite outcome.
#' @param Y2.var Character scalar. Stage-2 time/outcome component to be imputed for censored stage-2 subjects.
#' @param delta.var Character scalar. Stage-2 event indicator column name (1=observed, 0=censored).
#' @param OY.var Character scalar. Overall outcome column name (used when \code{eta2==0} and observed).
#'
#' @param A1.var Character scalar. Stage-1 treatment column name.
#' @param A2.var Character scalar. Stage-2 treatment column name.
#'
#' @param names.var1 Character vector. Covariate names available at stage 1 (used when \code{usecov=TRUE}).
#' @param names.var2 Character vector. Covariate names available at stage 2 (used when \code{usecov=TRUE}).
#'
#' @param exact1.vars,exact2.vars Character vectors. Variables used for exact matching at stage 1 / stage 2.
#' Default is none.
#'
#' @param usecov Logical. If TRUE, matching uses the covariates in \code{names.var1} / \code{names.var2}.
#' If FALSE, matching uses learned censoring scores (\code{pscens} and/or \code{pgcens}).
#'
#' @param useds Logical. Reserved for future use (currently not used in the provided implementation).
#'
#' @param adjustdelta1 Logical. If TRUE, modifies the stage-1 event indicator to treat certain stage-2
#' events as censored by creating \code{deltaadj}. (Current code sets \code{deltaadj=0} for
#' \code{eta2==1 & delta==1}.)
#'
#' @param cores Integer. Requested number of cores for downstream scoring routines (if supported).
#' @param tau Optional numeric. Truncation horizon passed to scoring routines for mean survival calculations.
#' @param sl.seed Integer. RNG seed passed to SuperLearner-based scoring.
#'
#' @param A.SL.library1,A.SL.library2 Character vectors. SuperLearner libraries for stage-1 and stage-2
#' censoring propensity models (passed to \code{ComputeScores}).
#' @param Y.SL.library Character vector. Learners for survivalSL prognostic modeling (passed through).
#' @param A.method,Y.method Optional characters. Risk/metric identifiers for SuperLearner / survivalSL.
#'
#' @param param.weights.fix,param.weights.init,optim.method,maxit,penalty1,penalty2,param.tune
#' Tuning and optimization controls passed to \code{ComputeScores} / survivalSL scoring.
#'
#' @param ngrid Integer. Number of grid points for survival-curve integration in prognostic scoring.
#'
#' @param pscens,pgcens Logical. Whether to compute and use censoring propensity (\code{pscens}) and/or
#' censoring prognostic (\code{pgcens}) scores when \code{usecov=FALSE}.
#'
#' @param plotps Logical. If TRUE, produces diagnostic propensity plots for censoring scores.
#'
#' @param model.pg Character. Prognostic modeling family when \code{superLearn=FALSE} ("cox" or "aft").
#' @param standardize Logical. Whether to standardize covariates for glmnet when \code{superLearn=FALSE}.
#' @param superLearn Logical. If TRUE, uses SuperLearner-based scoring; otherwise uses glm/glmnet-based scoring.
#' @param pslink Character. Link for binomial models ("logit" or "probit").
#'
#' @param distance Character. Distance type used by matching routines (e.g., "mahalanobis").
#' @param method Character. Matching method ("nearest" or "optimal").
#' @param K Integer. Donor ratio (number of matched donors per censored subject).
#' @param replacement Logical. Whether donors can be reused across matches.
#'
#' @return A data.frame containing the original data augmented with imputed/composite outcomes.
#' At minimum, the output includes:
#' \itemize{
#'   \item \code{compY2}: imputed/observed stage-2 component
#'   \item \code{compOY}: composite overall outcome used for stage-1 imputation
#'   \item \code{delta1}, \code{delta2}: derived indicators of composite outcome observability
#' }
#' If \code{usecov=FALSE}, additional columns containing raw and standardized censoring scores may also be present,
#' depending on \code{pscens} and \code{pgcens}.
#'
#' @seealso \code{\link{ComputeScores}}, \code{impute_censored_stage2}, \code{impute_censored_stage1}
#' @export


impute_censored_outcomes <- function(
    data,
    id.var, eta2.var,
    Y1.var, Y2.var,
    delta.var, OY.var,
    A1.var, A2.var,
    names.var1, names.var2,
    exact1.vars = character(0),
    exact2.vars = character(0),
    usecov        = TRUE,
    useds         = FALSE,
    adjustdelta1  = FALSE,
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
    pscens        = TRUE,
    pgcens        = TRUE,
    param.tune    = NULL,
    plotps        = FALSE,
    model.pg      = "cox",      # "cox" or "aft"
    standardize   = FALSE,      # Standardize covariates for glmnet
    superLearn    = TRUE,       # Whether to use SuperLearner or glmnet
    pslink        = "logit",    # "logit" or "probit"
    distance      = "mahalanobis",
    method        = "nearest",  # "nearest" or "optimal"
    K             = 3,          # donor ratio
    replacement   = TRUE
) {
  # ---- fast validation ----
  stopifnot(is.data.frame(data))
  req_cols <- c(id.var, eta2.var, Y1.var, Y2.var, delta.var, OY.var, A1.var, A2.var)
  miss_cols <- setdiff(req_cols, names(data))
  if (length(miss_cols) > 0L) {
    stop("Missing required columns in data: ", paste(miss_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.character(names.var1) || length(names.var1) < 1L) stop("names.var1 must be a non-empty character vector.", call. = FALSE)
  if (!is.character(names.var2) || length(names.var2) < 1L) stop("names.var2 must be a non-empty character vector.", call. = FALSE)
  if (!all(names.var1 %in% names(data))) stop("Some names.var1 columns not found in data.", call. = FALSE)
  if (!all(names.var2 %in% names(data))) stop("Some names.var2 columns not found in data.", call. = FALSE)

  # helpers must exist (they can be internal package functions)
  if (!exists("impute_censored_stage2", mode = "function")) {
    stop("Function 'impute_censored_stage2()' not found. It must be included in your package.", call. = FALSE)
  }
  if (!exists("impute_censored_stage1", mode = "function")) {
    stop("Function 'impute_censored_stage1()' not found. It must be included in your package.", call. = FALSE)
  }

  # timing helpers (no hard dependency)
  tick <- function(...) if (requireNamespace("tictoc", quietly = TRUE)) tictoc::tic(...)
  tock <- function(...) if (requireNamespace("tictoc", quietly = TRUE)) tictoc::toc()
  # if (requireNamespace("cobalt", quietly = TRUE)) {
  #   cobalt::love.plot(...)
  # }

  # safe logit with clipping
  safe_qlogis <- function(p, eps = 1e-6) {
    p <- as.numeric(p)
    p <- pmax(pmin(p, 1 - eps), eps)
    stats::qlogis(p)
  }

  # stable message helper
  say <- function(...) message(...)

  df <- data

  # ---- optional delta adjustment ----
  if (adjustdelta1) {
    df$deltaadj <- df[[delta.var]]
    df$deltaadj[df[[eta2.var]] == 1 & df[[delta.var]] == 1] <- 0
  }

  # ---- build matching formulas ----
  if (usecov) {
    f1 <- paste(names.var1, collapse = " + ")
    f2 <- paste(names.var2, collapse = " + ")
    # NOTE: formula variables delta1 and (delta.var) are created/used later; stage-2 formula uses delta.var directly.
    formula1 <- stats::as.formula(paste0("(delta1 == 0) ~ ", f1))
    formula2 <- stats::as.formula(paste0("(", delta.var, " == 0) ~ ", f2))
  } else {
    # score-based matching: rely on pscens/pgcens-derived columns we will create below
    if (pscens && pgcens) {
      formula1 <- stats::as.formula("(delta1 == 0) ~ pscens1 + pgcens1")
      formula2 <- stats::as.formula(paste0("(", delta.var, " == 0) ~ pscens2 + pgcens2"))
    } else if (pgcens) {
      formula1 <- stats::as.formula("(delta1 == 0) ~ pgcens1")
      formula2 <- stats::as.formula(paste0("(", delta.var, " == 0) ~ pgcens2"))
    } else if (pscens) {
      formula1 <- stats::as.formula("(delta1 == 0) ~ pscens1")
      formula2 <- stats::as.formula(paste0("(", delta.var, " == 0) ~ pscens2"))
    } else {
      stop("When usecov=FALSE, at least one of pscens or pgcens must be TRUE.", call. = FALSE)
    }
  }

  # ---- stage 2 subset ----
  df2 <- df[df[[eta2.var]] == 1, , drop = FALSE]

  # ---- compute censoring scores for stage 2 (optional) ----
  if (!usecov) {
    if (!exists("ComputeScores", mode = "function")) {
      stop("ComputeScores() not found. It must be included in your package when usecov=FALSE.", call. = FALSE)
    }

    ds2 <- ComputeScores(
      data         = df2,
      id           = id.var,
      Y            = Y2.var,
      event        = delta.var,
      X            = names.var2,
      A            = A2.var,
      censmod      = TRUE,
      doublepg     = FALSE,
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
      pscens       = pscens,
      pgcens       = pgcens,
      param.tune   = param.tune,
      model.pg     = model.pg,
      standardize  = standardize,
      superLearn   = superLearn,
      pslink       = pslink
    )

    ds2 <- as.data.frame(ds2)
    if (!id.var %in% names(ds2)) names(ds2)[1] <- id.var

    # Standardized score frame for matching covariates (pscens2 / pgcens2)
    dsp2 <- ds2[, id.var, drop = FALSE]

    if (pscens) {
      if (!"pscens" %in% names(ds2)) stop("ComputeScores did not return 'pscens' for stage 2.", call. = FALSE)
      ds2$pscens <- as.numeric(ds2$pscens)

      if (plotps) {
        if (!exists("propensityplot", mode = "function")) {
          stop("plotps=TRUE requires a function 'propensityplot()' available in the package.", call. = FALSE)
        }
        p <- propensityplot(ps = ds2[["pscens"]], A = df2[[delta.var]])
        if (requireNamespace("ggplot2", quietly = TRUE)) {
          p <- p + ggplot2::labs(fill = "Event", colour = "Event")
        }
        print(p)
      }

      # logit-transform then z-score for matching distance
      dsp2$pscens2 <- as.numeric(scale(safe_qlogis(ds2$pscens)))
      # keep raw probability too, but do not overwrite pscens
      ds2$psprobcens2 <- ds2$pscens
    }

    if (pgcens) {
      if (!"pgcens" %in% names(ds2)) stop("ComputeScores did not return 'pgcens' for stage 2.", call. = FALSE)
      ds2$pgcens <- as.numeric(ds2$pgcens)
      dsp2$pgcens2 <- as.numeric(scale(ds2$pgcens))
      ds2$pgprobcens2 <- ds2$pgcens
    }

    # keep only ID + raw score columns we created, to avoid clutter/duplication
    keep_raw2 <- c(id.var,
                   if (pscens) "psprobcens2",
                   if (pgcens) "pgprobcens2")
    ds2_out <- ds2[, keep_raw2, drop = FALSE]

    # merge back into df (all subjects; stage-2 non-entrants get NAs)
    df <- merge(df, ds2_out, by = id.var, all.x = TRUE)
    df <- merge(df, dsp2,   by = id.var, all.x = TRUE)
  }

  # refresh df2 after merges
  df2 <- df[df[[eta2.var]] == 1, , drop = FALSE]

  # ---- stage 2 matching / imputation ----
  tick("Stage 2 matching")

  res2 <- impute_censored_stage2(
    dat        = df2,
    id.var     = id.var,
    delta.var  = delta.var,
    OY.var     = OY.var,
    Y2.var     = Y2.var,
    formula2   = formula2,
    exact.vars = exact2.vars,
    method     = method,
    distance   = distance,
    k          = K,
    replace    = replacement,
    caliper    = NULL,
    aggregate  = "mean"
  )

  say("With ", nrow(df2), " subjects who entered stage 2, ",
      res2$n_imputed, " censored subjects were matched/imputed out of ",
      res2$n_censored, ".")

  comp2 <- res2$data_merged
  if (!all(c(id.var, "compY2") %in% names(comp2))) {
    stop("impute_censored_stage2() must return data_merged with columns ", id.var, " and compY2.", call. = FALSE)
  }

  df <- merge(df, comp2[, c(id.var, "compY2")], by = id.var, all.x = TRUE)

  # For uncensored stage-2 subjects, compY2 is observed Y2
  df$compY2[df[[delta.var]] == 1] <- df[[Y2.var]][df[[delta.var]] == 1]

  # Drop stage-2 entrants that still have missing compY2 after imputation (structural failure)
  df <- df[!(df[[eta2.var]] == 1 & is.na(df$compY2)), , drop = FALSE]

  # Composite outcome for stage 1
  df$compOY <- df$compY2 + df[[Y1.var]]
  df$compOY[df[[eta2.var]] == 0 & df[[delta.var]] == 1] <- df[[OY.var]][df[[eta2.var]] == 0 & df[[delta.var]] == 1]

  # Derived censoring indicators for stage-1 matching
  df$delta1 <- as.numeric(!is.na(df$compOY))
  df$delta2 <- as.numeric(!is.na(df$compOY))

  tock()

  # ---- censoring scores for stage 1 (optional) ----
  if (!usecov) {
    ds1 <- ComputeScores(
      data         = df,
      id           = id.var,
      Y            = if (adjustdelta1) Y1.var else OY.var,
      event        = if (adjustdelta1) "deltaadj" else delta.var,
      X            = names.var1,
      A            = A1.var,
      doublepg     = FALSE,
      censmod      = TRUE,
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
      pscens       = pscens,
      pgcens       = pgcens,
      ngrid        = ngrid,
      param.tune   = param.tune,
      model.pg     = model.pg,
      standardize  = standardize,
      superLearn   = superLearn,
      pslink       = pslink
    )

    ds1 <- as.data.frame(ds1)
    if (!id.var %in% names(ds1)) names(ds1)[1] <- id.var

    dsp1 <- ds1[, id.var, drop = FALSE]

    if (pscens) {
      if (!"pscens" %in% names(ds1)) stop("ComputeScores did not return 'pscens' for stage 1.", call. = FALSE)
      ds1$pscens <- as.numeric(ds1$pscens)

      if (plotps) {
        if (!exists("propensityplot", mode = "function")) {
          stop("plotps=TRUE requires a function 'propensityplot()' available in the package.", call. = FALSE)
        }
        p <- propensityplot(ps = ds1[["pscens"]], A = df[[delta.var]])
        if (requireNamespace("ggplot2", quietly = TRUE)) {
          p <- p + ggplot2::labs(fill = "Event", colour = "Event")
        }
        print(p)
      }

      dsp1$pscens1 <- as.numeric(scale(safe_qlogis(ds1$pscens)))
      ds1$psprobcens1 <- ds1$pscens
    }

    if (pgcens) {
      if (!"pgcens" %in% names(ds1)) stop("ComputeScores did not return 'pgcens' for stage 1.", call. = FALSE)
      ds1$pgcens <- as.numeric(ds1$pgcens)
      dsp1$pgcens1 <- as.numeric(scale(ds1$pgcens))
      ds1$pgprobcens1 <- ds1$pgcens
    }

    keep_raw1 <- c(id.var,
                   if (pscens) "psprobcens1",
                   if (pgcens) "pgprobcens1")
    ds1_out <- ds1[, keep_raw1, drop = FALSE]

    df <- merge(df, ds1_out, by = id.var, all.x = TRUE)
    df <- merge(df, dsp1,   by = id.var, all.x = TRUE)
  }

  # ---- stage 1 matching / imputation ----
  tick("Stage 1 matching")

  res1 <- impute_censored_stage1(
    dat        = df,
    Id         = id.var,
    exact_vars = exact1.vars,
    formula    = formula1,
    death1     = "delta1",
    OY         = OY.var,
    y1         = Y1.var,
    distance   = distance,
    method     = method,
    k          = K,
    replace    = replacement,
    y_cols     = c("compOY"),
    aggregate  = "mean"
  )

  say("With ", nrow(df), " stage-1 entrants, ",
      res1$n_imputed, " censored were matched/imputed out of ",
      res1$n_censored, ".")

  tock()

  finaldf <- res1$data_merged
  finaldf
}



#' Impute censored stage-1 composite outcomes via constrained donor matching
#'
#' @description
#' For each stage-1 censored subject (\code{death1 == 0}), constructs an eligible donor pool of
#' uncensored subjects (\code{death1 == 1}) and imputes one or more outcome columns (default
#' \code{compOY}) by matching the censored subject to \code{k} donors using \pkg{MatchIt}.
#' Matching can be nearest-neighbor or optimal, with optional exact matching constraints.
#'
#' Donor eligibility is enforced by filtering to donors whose primary outcome \code{y_cols[1]} is at
#' least the recipient’s observed-time boundary \code{OY}. A second filter requires donors to satisfy
#' \code{is.na(y1) | y1 >= OY}, which is intended to ensure donors have follow-up long enough relative to
#' the censored subject.
#'
#' @details
#' The function loops over censored IDs. For each censored ID, it forms a temporary dataset consisting
#' of the focal censored unit plus its eligible donors, and runs \pkg{MatchIt} by defining
#' \code{tr = (death1 == 0)} (censored-as-treated). Donors are taken from the same \code{subclass} as the
#' focal censored unit.
#'
#' Donor outcomes are aggregated to produce the imputation:
#' \itemize{
#'   \item \code{"mean"}: arithmetic mean across donors,
#'   \item \code{"weighted"}: weighted mean using \code{weights} returned by \pkg{MatchIt},
#'   \item \code{"nearest"}: take the single nearest donor by \code{distance}, even if \code{k > 1}.
#' }
#'
#' If an ID fails to match (insufficient donors, matching error, etc.), it is skipped and its
#' imputation remains missing unless handled via \code{na_handling}.
#'
#' @param dat A data.frame containing the stage-1 cohort and all variables required for donor filtering,
#' matching, and outcome imputation.
#' @param Id Character scalar. Subject identifier column name.
#' @param exact_vars Optional exact matching specification passed to \pkg{MatchIt}. One of:
#' \code{NULL}, a one-sided formula (e.g., \code{~ sex + site}), or a character vector of column names.
#' @param formula A formula giving the matching covariates. Internally \code{update(formula, tr ~ .)}
#' is used to create a two-sided formula with \code{tr} as the matching "treatment" indicator.
#' @param death1 Character scalar. Column name for stage-1 event indicator. Convention:
#' \code{1 = observed event/outcome}, \code{0 = censored}.
#' @param OY Character scalar. Column name for observed-time boundary used to constrain donors for each censored subject.
#' @param y1 Character scalar. Additional time-like column used in donor filtering: donors must satisfy
#' \code{is.na(y1) | y1 >= OY}.
#' @param distance Character. Distance argument passed to \pkg{MatchIt} (e.g., \code{"mahalanobis"}).
#' @param method Matching method: \code{"nearest"} or \code{"optimal"}.
#' @param k Integer. Donor ratio (number of matched donors per censored subject).
#' @param replace Logical. Whether donors can be reused across matches (nearest-neighbor only).
#' @param caliper Optional numeric. Caliper passed to nearest-neighbor matching.
#' @param y_cols Character vector. Outcome column(s) to impute. The first element \code{y_cols[1]} is
#' also used in the donor eligibility filter \code{y_cols[1] >= OY}.
#' @param aggregate Aggregation rule for donor outcomes: \code{"mean"}, \code{"weighted"}, or \code{"nearest"}.
#' @param na_handling What to do if the primary imputed column remains missing: \code{"drop"} removes
#' those rows from \code{data_merged}; \code{"zero_weight"} keeps them (for downstream handling).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{imputed}: data.frame of imputed values keyed by \code{Id} (or \code{NULL} if none imputed),
#'   \item \code{data_merged}: input data with \code{y_cols} updated for successfully imputed IDs,
#'   \item \code{n_censored}: number of censored subjects targeted,
#'   \item \code{n_imputed}: number successfully imputed,
#'   \item \code{errors}: two-column matrix recording IDs and error messages from \pkg{MatchIt}.
#' }
#'
#' @seealso \code{\link[MatchIt]{matchit}}, \code{\link[MatchIt]{get_matches}}
#' @export



impute_censored_stage1 <- function(dat, Id,
                                   exact_vars = NULL,
                                   formula,                 # single RHS formula; we use update(formula, tr ~ .)
                                   death1,                  # column name for stage-1 event indicator (1=event, 0=censored)
                                   OY,                      # column name for observed-time boundary used in donor filter
                                   y1,                      # extra time variable used in donor filter
                                   distance = "mahalanobis",
                                   method = c("nearest","optimal"),
                                   k = 1, replace = TRUE, caliper = NULL,
                                   y_cols = c("compOY"),
                                   aggregate = c("mean","weighted","nearest"),
                                   na_handling = c("drop","zero_weight")) {

  method      <- match.arg(method)
  aggregate   <- match.arg(aggregate)
  na_handling <- match.arg(na_handling)

  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is required for impute_censored_stage1().", call. = FALSE)
  }

  .make_exact <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "formula")) return(x)
    if (is.character(x)) return(stats::reformulate(x))
    stop("exact_vars must be NULL, a one-sided formula, or a character vector", call. = FALSE)
  }

  .agg <- function(df, vars) {
    if (nrow(df) == 0L) return(stats::setNames(rep(NA_real_, length(vars)), vars))

    if (aggregate == "weighted" && "weights" %in% names(df) && all(is.finite(df$weights))) {
      w <- df$weights
      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0) w[] <- 1 / nrow(df) else w <- w / sw
      return(vapply(vars, function(v) stats::weighted.mean(df[[v]], w, na.rm = TRUE), numeric(1)))
    }

    if (aggregate == "nearest" && "distance" %in% names(df)) {
      df <- df[order(df$distance), , drop = FALSE]
      return(vapply(vars, function(v) as.numeric(df[[v]][1]), numeric(1)))
    }

    vapply(vars, function(v) mean(df[[v]], na.rm = TRUE), numeric(1))
  }

  stopifnot(is.data.frame(dat))
  needed <- unique(c(Id, death1, OY, y1, y_cols))
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0L) {
    stop("Missing required columns in 'dat': ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!is.character(y_cols) || length(y_cols) < 1L) {
    stop("y_cols must be a non-empty character vector.", call. = FALSE)
  }

  subHs <- dat

  dataComp1L <- NULL
  errorData  <- NULL

  cens_ids <- subHs[[Id]][subHs[[death1]] == 0]
  cens_ids <- cens_ids[!is.na(cens_ids)]

  for (idv in cens_ids) {
    tmp  <- subHs[subHs[[Id]] == idv, , drop = FALSE]
    tmp1 <- subHs[subHs[[death1]] == 1 & subHs[[y_cols[1]]] >= tmp[[OY]], , drop = FALSE]
    tmp1 <- tmp1[is.na(tmp1[[y1]]) | tmp1[[y1]] >= tmp[[OY]], , drop = FALSE]
    tmp  <- rbind(tmp, tmp1)

    tmp$tr <- tmp[[death1]] == 0
    tmp    <- tmp[!is.na(tmp$tr), , drop = FALSE]
    if (length(unique(tmp$tr)) < 2L) next

    mobj <- tryCatch({
      if (identical(method, "optimal")) {
        MatchIt::matchit(
          update(formula, tr ~ .),
          data = tmp, distance = distance,
          method = "optimal",
          exact  = .make_exact(exact_vars),
          ratio  = k
        )
      } else {
        MatchIt::matchit(
          update(formula, tr ~ .),
          data = tmp, distance = distance,
          method = "nearest",
          replace = replace,
          caliper = caliper,
          exact  = .make_exact(exact_vars),
          ratio  = k
        )
      }
    }, error = function(e) e)

    if (inherits(mobj, "error")) {
      message("ERROR: ", conditionMessage(mobj))
      errorData <- rbind(errorData, c(idv, paste("ERROR:", conditionMessage(mobj))))
      next
    }

    mm <- as.data.frame(MatchIt::get_matches(mobj, data = tmp))

    sc_focal <- unique(mm$subclass[mm[[Id]] == idv & mm$tr])
    if (length(sc_focal) == 0L) next

    donors <- mm[mm$subclass %in% sc_focal & mm[[death1]] == 1, , drop = FALSE]
    if (nrow(donors) == 0L) next

    agg_vals <- .agg(donors, y_cols)
    dataComp1L <- rbind(dataComp1L, c(idv, agg_vals))
  }

  out <- subHs

  if (!is.null(dataComp1L) && nrow(dataComp1L) > 0) {
    dataComp1L <- as.data.frame(dataComp1L, stringsAsFactors = FALSE)
    colnames(dataComp1L) <- c(Id, "compOY") #, "compOS_2LChemo", "compOS_2LCheMon"
    dataComp1L$compOY <- as.numeric(dataComp1L$compOY)

    idx_match <- match(dataComp1L[[Id]], out[[Id]])
    out$compOY[idx_match] <- dataComp1L$compOY
    # out$CompOS_2LChemo1Lobs[idx_match] <- dataComp1L$compOS_2LChemo
    # out$CompOS_2LCheMon1Lobs[idx_match] <- dataComp1L$compOS_2LCheMon
  }

  out$compOY <- as.numeric(out$compOY)
  # out$CompOS_2LChemo1Lobs <- as.numeric(out$CompOS_2LChemo1Lobs)
  # out$CompOS_2LCheMon1Lobs <- as.numeric(out$CompOS_2LCheMon1Lobs)
  #
  bad <- is.na(out$compOY)
  if (any(bad)) {
    if (na_handling == "drop") {
      out <- out[!bad, , drop = FALSE]
      message(sprintf("Dropped %d row(s) with missing compOY", sum(bad)))
    } else {
      # keep but set a zero weight column if you track match weights elsewhere
      # here we simply keep the rows; user can filter downstream
      message(sprintf("Kept %d row(s) with missing compOY", sum(bad)))
    }
  }

  list(
    imputed     = if (is.null(dataComp1L)) NULL else dataComp1L,
    data_merged = out,
    n_censored  = length(cens_ids),
    n_imputed   = if (is.null(dataComp1L)) 0L else nrow(dataComp1L),
    errors      = errorData
  )
}


#' Impute censored stage-2 outcomes via constrained donor matching
#'
#' @description
#' For each stage-2 censored subject (\code{delta.var == 0}), constructs an eligible donor pool of
#' uncensored subjects (\code{delta.var == 1}) and imputes the stage-2 outcome \code{Y2.var},
#' returning the imputed/observed outcome as \code{compY2}. Matching is performed using \pkg{MatchIt}
#' with nearest-neighbor or optimal matching, optional exact matching, and configurable aggregation
#' across \code{k} donors.
#'
#' Donor eligibility is constrained by requiring donors to have observed overall time \code{OY.var}
#' at least as large as the recipient’s \code{OY.var}. This is intended to ensure donors have
#' follow-up long enough relative to the censored subject.
#'
#' @details
#' The function loops over censored IDs. For each ID, it forms a temporary dataset consisting of the
#' focal censored unit plus eligible donors and runs \pkg{MatchIt} by defining
#' \code{tr = (delta.var == 0)} (censored-as-treated). Donors are then taken from the matched
#' \code{subclass} containing the focal censored unit.
#'
#' Donor outcomes \code{Y2.var} are aggregated according to \code{aggregate}:
#' \itemize{
#'   \item \code{"mean"}: arithmetic mean,
#'   \item \code{"weighted"}: weighted mean using \code{weights} from \pkg{MatchIt},
#'   \item \code{"nearest"}: take the single nearest donor by \code{distance}.
#' }
#'
#' For uncensored subjects (\code{delta.var == 1}), \code{compY2} is set to the observed \code{Y2.var}.
#'
#' @param dat A data.frame containing stage-2 entrants and required variables for matching and filtering.
#' @param id.var Character scalar. Subject identifier column name.
#' @param delta.var Character scalar. Stage-2 event indicator column name. Convention:
#' \code{1 = observed event/outcome}, \code{0 = censored}.
#' @param OY.var Character scalar. Column name for observed-time boundary used to constrain donors.
#' @param Y2.var Character scalar. Stage-2 outcome column to impute.
#' @param formula2 A formula giving matching covariates. Internally \code{update(formula2, tr ~ .)}
#' is used to create a two-sided formula with \code{tr} as the matching "treatment" indicator.
#' @param exact.vars Optional exact matching specification passed to \pkg{MatchIt}: \code{NULL}, a
#' one-sided formula, or a character vector of names.
#' @param method Matching method: \code{"nearest"} or \code{"optimal"}.
#' @param distance Character. Distance argument passed to \pkg{MatchIt} (e.g., \code{"mahalanobis"}).
#' @param k Integer. Donor ratio (number of matched donors per censored subject).
#' @param replace Logical. Whether donors can be reused across matches (nearest-neighbor only).
#' @param caliper Optional numeric. Caliper for nearest-neighbor matching.
#' @param aggregate Aggregation rule: \code{"mean"}, \code{"weighted"}, or \code{"nearest"}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{imputed}: data.frame of imputed \code{compY2} values keyed by \code{id.var} (or \code{NULL} if none),
#'   \item \code{data_merged}: the input \code{dat} augmented with \code{compY2},
#'   \item \code{n_censored}: number of censored subjects targeted,
#'   \item \code{n_imputed}: number successfully imputed,
#'   \item \code{errors}: two-column matrix recording IDs and error messages from \pkg{MatchIt}.
#' }
#'
#' @seealso \code{\link[MatchIt]{matchit}}, \code{\link[MatchIt]{get_matches}}
#' @export
impute_censored_stage2 <- function(dat,
                                   id.var,
                                   delta.var,          # 1 = event, 0 = censored
                                   OY.var,
                                   Y2.var,
                                   formula2,
                                   exact.vars = NULL,
                                   method = c("nearest","optimal"),
                                   distance = "mahalanobis",
                                   k = 1,
                                   replace = TRUE,
                                   caliper = NULL,
                                   aggregate = c("mean","weighted","nearest")) {

  method    <- match.arg(method)
  aggregate <- match.arg(aggregate)


  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is required for impute_censored_stage2().", call. = FALSE)
  }

  .make_exact <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "formula")) return(x)
    if (is.character(x)) return(stats::reformulate(x))
    stop("exact.vars must be NULL, a one-sided formula, or a character vector", call. = FALSE)
  }

  .agg <- function(df, var) {
    if (nrow(df) == 0L) return(NA_real_)
    if (aggregate == "weighted" && "weights" %in% names(df) && all(is.finite(df$weights))) {
      w <- df$weights
      sw <- sum(w)
      if (is.finite(sw) && sw > 0) w <- w / sw else w[] <- 1 / nrow(df)
      return(stats::weighted.mean(df[[var]], w, na.rm = TRUE))
    }
    if (aggregate == "nearest" && "distance" %in% names(df)) {
      df <- df[order(df$distance), , drop = FALSE]
      return(as.numeric(df[[var]][1]))
    }
    mean(df[[var]], na.rm = TRUE)
  }

  stopifnot(is.data.frame(dat))
  needed <- unique(c(id.var, delta.var, OY.var, Y2.var))
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0L) {
    stop("Missing required columns in 'dat': ", paste(miss, collapse = ", "), call. = FALSE)
  }

  dataComp2L <- NULL
  errorData  <- NULL

  cens_ids <- dat[[id.var]][dat[[delta.var]] == 0]
  cens_ids <- cens_ids[!is.na(cens_ids)]

  for (idv in cens_ids) {
    tmp  <- dat[dat[[id.var]] == idv, ]
    tmp1 <- dat[dat[[delta.var]] == 1 & dat[[OY.var]] >= tmp[[OY.var]], ]
    tmp  <- rbind(tmp, tmp1)

    # censored-as-treated
    tmp$tr <- tmp[[delta.var]] == 0
    tmp    <- tidyr::drop_na(tmp, tr)
    if (dplyr::n_distinct(tmp$tr) < 2L) next

    mobj <- tryCatch({
      if (identical(method, "optimal")) {
        MatchIt::matchit(
          update(formula2, tr ~ .),
          data = tmp, distance = distance,
          method = "optimal",
          exact  = .make_exact(exact.vars),
          ratio  = k
        )
      } else {
        MatchIt::matchit(
          update(formula2, tr ~ .),
          data = tmp, distance = distance,
          method = "nearest",
          replace = replace,
          caliper = caliper,
          exact  = .make_exact(exact.vars),
          ratio  = k
        )
      }
    }, error = function(e) e)

    if (inherits(mobj, "error")) {
      message("ERROR: ", conditionMessage(mobj))
      errorData <- rbind(errorData, c(idv, paste("ERROR:", conditionMessage(mobj))))
      next
    }

    mm <- as.data.frame(MatchIt::get_matches(mobj, data = tmp))

    # same subclass as focal censored unit
    sc_focal <- unique(mm$subclass[mm[[id.var]] == idv & mm$tr])
    if (length(sc_focal) == 0L) next

    donors <- mm[mm$subclass %in% sc_focal & mm[[delta.var]] == 1, , drop = FALSE]
    if (nrow(donors) == 0L) next

    imputed_val <- .agg(donors, Y2.var)
    dataComp2L  <- rbind(dataComp2L, c(idv, imputed_val))
  }

  dataComp2L <- as.data.frame(dataComp2L, stringsAsFactors = FALSE)
  if (!is.null(dataComp2L) && nrow(dataComp2L) > 0) {
    colnames(dataComp2L) <- c(id.var, "compY2")
    dataComp2L$compY2  <- as.numeric(dataComp2L$compY2)

    outDat <- merge(dat, dataComp2L, by = id.var, all.x = TRUE)

    is_event <- outDat[[delta.var]] == 1
    outDat$compY2[is_event] <- outDat[[Y2.var]][is_event]

    outDat <- outDat[!is.na(outDat$compY2), ]
  } else {
    outDat <- dat
    outDat$compY2 <- NA_real_
  }

  list(imputed = dataComp2L,
       data_merged = outDat,
       n_censored = length(cens_ids),
       n_imputed  = if (is.null(dataComp2L)) 0L else nrow(dataComp2L),
       errors = errorData)
}
