#' Two-stage DTR learning/evaluation pipeline with cross-fitting
#'
#' @description
#' \code{simrun()} runs one replicate (indexed by \code{jj}) of a two-stage dynamic
#' treatment regime (DTR) learning and evaluation pipeline using outer cross-fitting
#' (default 5 folds). Within each fold, the function imputes censoring-related
#' composite outcomes on the training partition, constructs double scores
#' (propensity/prognostic scores) when requested, performs matching-based
#' pseudo-outcome construction for stage 2 and stage 1, fits random-forest
#' policies via \code{\link{rfdtr}}, applies the learned rules to the held-out
#' fold, and returns fold-level value summaries and benchmark comparisons.
#'
#' @details
#' The outer loop is a 5-fold split created by \code{caret::createFolds()} on the
#' stage-1 observed treatment \code{A1.var}. For each fold:
#' \enumerate{
#'   \item Split into training and test partitions.
#'   \item Standardize selected time variables in training (\code{OY.sd}, \code{Y1.sd}, \code{Y2.sd}).
#'   \item Impute censoring-related composite outcomes on the training fold via
#'         \code{\link{impute_censored_outcomes}} (internally calling stage-2 then stage-1 imputation).
#'   \item Optionally compute double scores on the imputed training fold via
#'         \code{\link{get_doublescores}}.
#'   \item Stage 2 (among \code{eta2==1}): build matched opposite-arm pseudo-outcomes via
#'         \code{\link{matchpotential_DTR}}, tune/fit a random-forest policy via
#'         \code{\link{rfdtr}}, and compute an estimated stage-2 optimal value by
#'         substituting matched outcomes for units assigned to a non-optimal observed treatment.
#'   \item Stage 1: construct a pseudo-outcome \code{Rtilde} (including the estimated
#'         stage-2 optimal value for stage-2 entrants), perform matching via
#'         \code{\link{matchpotential_DTR}}, and tune/fit the stage-1 policy via \code{\link{rfdtr}}.
#'   \item Apply learned policies to the held-out fold and compute fold-level summaries:
#'         value estimates (means and medians, capped at \code{cap_months}), agreement
#'         measures (\code{CCR1}, \code{CCR2}), and observed/uniform benchmark regimes.
#' }
#'
#' @param data A \code{data.frame} containing all required columns used by
#'  `Get_data()`, matching, scoring, and policy learning.
#'
#' @param id.var Name of unique subject identifier column.
#' @param eta2.var Name of stage-2 entry indicator (1 = enters stage 2).
#' @param Y1.var Name of stage-1 observed time/outcome component.
#' @param Y2.var Name of stage-2 observed time/outcome component.
#' @param delta.var Name of event indicator (1 = event, 0 = censored).
#' @param OY.var Name of overall observed time (used as a boundary in donor filtering).
#' @param A1.var Name of stage-1 observed treatment/action.
#' @param A2.var Name of stage-2 observed treatment/action.
#'
#' @param exact1.vars Character vector of variables used for exact matching at stage 1.
#'   Use \code{NULL} or \code{character(0)} to disable exact matching.
#' @param exact2.vars Character vector of variables used for exact matching at stage 2.
#'   Use \code{NULL} or \code{character(0)} to disable exact matching.
#'
#' @param names.var1 Character vector of covariate names used at stage 1 (scores/matching/model).
#' @param names.var2 Character vector of covariate names used at stage 2 (scores/matching/model).
#' @param usecov Logical; if \code{TRUE}, matching formulas are built from \code{names.var1}/\code{names.var2}.
#'   If \code{FALSE}, matching may rely on score variables (e.g., pg/ps) produced by the pipeline.
#'
#' @param cores Integer number of workers used for internal parallel steps.
#' @param tau Optional horizon(s) passed to score construction; in the shown implementation,
#'   \code{rtau} is set to \code{cap_months} for both stages.
#' @param sl.seed Seed passed to SuperLearner / score construction components.
#'
#' @param A.SL.library Candidate learners for treatment (propensity) score models.
#' @param Y.SL.library Candidate learners for outcome/censoring models used in \code{ComputeScores()}.
#' @param A.method Optional optimization method for propensity fitting in score construction.
#' @param Y.method Optional optimization method for outcome fitting in score construction.
#'
#' @param plotps Logical; if \code{TRUE}, generate propensity score overlap plots where implemented.
#' @param ngrid Integer; grid size for internal optimization routines used by score construction.
#' @param param.tune Optional tuning structure passed into score construction.
#' @param param.weights.fix,param.weights.init Optional fixed/initial weights passed to score construction.
#' @param optim.method Optimization method used in score construction (e.g., \code{"Nelder-Mead"}).
#' @param maxit Maximum iterations for score construction optimizers.
#'
#' @param usepenalty Logical; if \code{TRUE}, define \code{SL.glmnet.tune1}/\code{SL.glmnet.tune2} wrappers
#'   with \code{penalty.factor} for stage-specific propensity estimation.
#' @param runseed Integer vector of seeds; replicate \code{jj} uses \code{runseed[jj]}.
#'
#' @param useds Logical; if \code{TRUE}, compute/use double scores via \code{\link{get_doublescores}}.
#' @param Taus Optional matrix/data.frame of candidate horizons (legacy/compatibility).
#' @param adjustdelta1 Logical; if \code{TRUE}, modifies stage-1 censoring indicator as in
#'   \code{impute_censored_outcomes()} and \code{get_doublescores()}.
#' @param modeltype Character; model type for policy learning in \code{\link{rfdtr}} (e.g., \code{"ranger"}).
#' @param usecv Logical; if \code{TRUE}, tune policy learner with internal CV in \code{\link{rfdtr}}.
#' @param doublepg Logical; if \code{TRUE}, use arm-specific prognostic scores (pg0/pg1) where applicable.
#'
#' @param model.pg Character; prognostic model family (e.g., \code{"cox"} or \code{"aft"}) for score construction.
#' @param standardize Logical; if \code{TRUE}, standardize covariates for \code{glmnet} components.
#' @param superLearn Logical; if \code{TRUE}, use SuperLearner; if \code{FALSE}, use glmnet-only path
#'   where implemented by \code{ComputeScores()}.
#' @param pslink Character; link for propensity estimation in score construction (e.g., \code{"logit"}).
#' @param pglink Character; distribution/link for prognostic score modeling (project-specific; e.g., \code{"lognormal"}).
#'
#' @param distance Matching distance used by \code{\link{matchpotential_DTR}} (e.g., \code{"mahalanobis"}).
#' @param method Matching method used in \code{\link{matchpotential_DTR}} (e.g., \code{"nearest"} or \code{"full"}).
#' @param K Integer donor ratio for matching (number of donors per anchor).
#' @param replacement Logical; whether to allow donor reuse in nearest-neighbor matching.
#'
#' @param cap_months Numeric truncation cap (months). Value summaries are computed after capping at this horizon.
#' @param plotbalance Logical; if \code{TRUE} and \pkg{cobalt} is installed, prints Love plots of
#'   standardized mean differences for covariate balance after each matching step. Default \code{FALSE}.
#'
#' @return
#' A \code{data.frame} with one row per outer fold (default 5). Columns include:
#' \itemize{
#'   \item Fold-level value summaries (means and medians) under the learned DTR and under several
#'         alternative constructions produced inside the pipeline (e.g., \code{Ttot}, \code{T2}, \code{T2.rf},
#'         \code{TtotOSml}, \code{Ttot1LOSmed}, etc.).
#'   \item Agreement metrics: \code{CCR1} (stage-1 agreement) and \code{CCR2} (stage-2 agreement among \code{eta2==1}).
#'   \item Observed and uniform-regime benchmarks (means/medians) at stage 2 and overall.
#'   \item Diagnostics stored as list-columns, including RF best hyperparameters and predicted optimal actions.
#' }
#' The returned object also includes \code{iter} (replicate index) and \code{fold}.
#' If a fatal error occurs, \code{NULL} is returned.
#'
#' @seealso
##' `Get_data` `Get_data`,
#' \code{\link{impute_censored_outcomes}},
#' \code{\link{get_doublescores}},
#' \code{\link{matchpotential_DTR}},
#' \code{\link{rfdtr}}.
#'
#' @examples
#' \dontrun{
#' # Requires project-specific helpers: Get_data(), impute_censored_outcomes(),
#' # get_doublescores(), matchpotential_DTR(), rfdtr(), and my_score.Surv().
#' #
#' # res <- simrun(
#' #   data = your_data,
#' #   runseed = 2025,
#' #   cores = 2
#' # )
#' # head(res)
#' }
#'
#' @export

Drmatch <- function(
    data                       = data,
    id.var                     = 'PatientID',
    eta2.var                   = 'eta2',
    Y1.var                     = 'OS_time.1L',
    Y2.var                     = 'OS_time.2L',
    delta.var                  = 'deathInd.raw',
    OY.var                     = 'OS_time',
    A1.var                     = 'txgroup1L.sd',
    A2.var                     = 'txgroup2L0.sd',
    exact1.vars                = c('txgroup1L.sd','gender.sd', 'ECOG1st0.sd', 'ECOG1st1.sd'),
    exact2.vars                = c('txgroup1L.sd','gender.sd','ECOG2nd0.sd', 'ECOG2nd1.sd', 'txgroup2L0.sd'),
    names.var1                 = c('ageAt1L', 'gender.sd', 'ECOG1st0.sd', 'ECOG1st1.sd',
                                   'Albumin1st', 'Lymphocyte1st'),
    names.var2                 = c('ageAt1L', 'gender.sd', 'ECOG2nd0.sd', 'Lymphocyte2nd',
                                   'OS_time.1L', 'ECOG2nd1.sd','Albumin2nd', 'Lymphocyte2nd'),
    usecov                     = FALSE,
    cores                      = 5,
    tau                        = NULL,
    sl.seed                    = 1234,
    A.SL.library               = c("SL.glm", "SL.glmnet", "SL.ranger"),
    Y.SL.library               = c("LIB_COXlasso","LIB_COXen","LIB_AFTggamma", "LIB_RSF"),
    A.method                   = "method.NNloglik",
    Y.method                   = "ibll",
    plotps                     = FALSE,
    ngrid                      = 1000,
    param.tune                 = NULL,
    param.weights.fix          = NULL,
    param.weights.init         = NULL,
    optim.method               = 'Nelder-Mead',
    maxit                      = 10000,
    usepenalty                 = FALSE,
    runseed                    = 2025,
    useds                      = TRUE,
    Taus                       = NULL,
    adjustdelta1               = FALSE,
    modeltype                  = "ranger",
    usecv                      = TRUE,
    doublepg                   = TRUE,
    model.pg                   = "cox",
    standardize                = FALSE,
    superLearn                 = TRUE,
    pslink                     = "logit",
    pglink                     = "lognormal",
    distance                   = 'mahalanobis',
    method                     = 'nearest',
    K                          = 3,
    replacement                = TRUE,
    cap_months                 = 36,
    plotbalance                = FALSE
) {

  fit <- tryCatch({

    set.seed(runseed)

    # -------------------------
    # Full-data prep
    # -------------------------
    mldata <- data # Get_data(data)
    rtau <- c(cap_months, cap_months)

    # SL.glmnet.tune1 <- function(...) SL.glmnet(..., penalty.factor = penalty1)
    # SL.glmnet.tune2 <- function(...) SL.glmnet(..., penalty.factor = penalty2)
    # assign("SL.glmnet.tune1", SL.glmnet.tune1, envir = .GlobalEnv)
    # assign("SL.glmnet.tune2", SL.glmnet.tune2, envir = .GlobalEnv)

    SL.glmnet.tune1 <- function(...) SuperLearner::SL.glmnet(..., penalty.factor = penalty1)
    SL.glmnet.tune2 <- function(...) SuperLearner::SL.glmnet(..., penalty.factor = penalty2)

    if (usepenalty) {
      penalty2 <- rep(1, length(names.var2))
      penalty1 <- rep(1, length(names.var1))
      names(penalty2) <- names.var2
      names(penalty1) <- names.var1
      penalty2[c(1:3)] <- 0
      penalty1[c(1)] <- 0
      A.SL.library1 <- c("SL.glmnet.tune1", A.SL.library)
      A.SL.library2 <- c("SL.glmnet.tune2", A.SL.library)
    } else {
      penalty1 <- NULL
      penalty2 <- NULL
      A.SL.library1 <- A.SL.library
      A.SL.library2 <- A.SL.library
    }

    # standardize internal time variables on full training data
    mldata$OY.sd <- as.numeric(scale(mldata[[OY.var]]))
    mldata$Y1.sd <- as.numeric(scale(mldata[[Y1.var]]))
    mldata$Y2.sd <- as.numeric(scale(mldata[[Y2.var]]))

    # -------------------------
    # Imputation + double scores
    # -------------------------
    imptrain <- impute_censored_outcomes(
      data         = mldata,
      id.var       = id.var,
      eta2.var     = eta2.var,
      Y1.var       = Y1.var,
      Y2.var       = Y2.var,
      delta.var    = delta.var,
      OY.var       = OY.var,
      A1.var       = A1.var,
      A2.var       = A2.var,
      names.var1   = if (adjustdelta1) c(names.var1, A1.var, 'Y1.sd') else c(names.var1, A1.var, 'OY.sd'),
      names.var2   = c(names.var2, A1.var, A2.var, 'Y1.sd', 'Y2.sd'),
      exact1.vars  = exact1.vars,
      exact2.vars  = exact2.vars,
      usecov       = FALSE,
      useds        = FALSE,
      adjustdelta1 = adjustdelta1,
      cores        = cores,
      tau          = NULL,
      sl.seed      = 123,
      A.SL.library1 = A.SL.library1,
      A.SL.library2 = A.SL.library2,
      Y.SL.library  = Y.SL.library,
      A.method     = A.method,
      Y.method     = Y.method,
      param.weights.fix  = param.weights.fix,
      param.weights.init = param.weights.init,
      optim.method = optim.method,
      maxit        = 10000,
      penalty1     = penalty1,
      penalty2     = penalty2,
      ngrid        = 2000,
      pscens       = TRUE,
      pgcens       = FALSE,
      param.tune   = param.tune,
      plotps       = plotps,
      model.pg     = model.pg,
      standardize  = standardize,
      superLearn   = superLearn,
      pslink       = pslink,
      distance     = distance,
      method       = method,
      K            = K,
      replacement  = replacement
    )

    MLdata <- get_doublescores(
      data         = imptrain,
      id.var       = id.var,
      eta2.var     = eta2.var,
      Y1.var       = Y1.var,
      Y2.var       = Y2.var,
      delta.var    = delta.var,
      OY.var       = OY.var,
      A1.var       = A1.var,
      A2.var       = A2.var,
      names.var1   = names.var1,
      names.var2   = names.var2,
      Xtrt1        = NULL,
      Xtrt2        = NULL,
      useds        = TRUE,
      cores        = cores,
      tau          = NULL,
      sl.seed      = 123,
      A.SL.library1 = A.SL.library1,
      A.SL.library2 = A.SL.library2,
      Y.SL.library  = Y.SL.library,
      A.method     = A.method,
      Y.method     = Y.method,
      param.weights.fix  = param.weights.fix,
      param.weights.init = param.weights.init,
      optim.method = optim.method,
      maxit        = maxit,
      penalty1     = penalty1,
      penalty2     = penalty2,
      ngrid        = ngrid,
      censmod      = FALSE,
      doublepg     = TRUE,
      param.tune   = param.tune,
      adjustdelta1 = adjustdelta1,
      plotps       = plotps,
      model.pg     = model.pg,
      standardize  = standardize,
      superLearn   = superLearn,
      pslink       = pslink,
      pglink       = pglink
    )

    cl <- parallel::makeCluster(cores)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      try(foreach::registerDoSEQ(), silent = TRUE)
    }, add = TRUE)
    doSNOW::registerDoSNOW(cl)

    # -------------------------
    # Stage 2
    # -------------------------
    idx2 <- which(MLdata[[eta2.var]] == 1)

    xx <- exact2.vars[!exact2.vars %in% c(A1.var, A2.var)]
    exact2.vars2 <- if (length(xx) == 0) NULL else xx

    dat2 <- MLdata[idx2, , drop = FALSE]
    dat2$compY2[dat2$compY2 > rtau[2]] <- rtau[2]

    subHs2.match <- matchpotential_DTR(
      dat         = dat2,
      txgroup     = A2.var,
      exact_vars  = exact2.vars2,
      compY       = 'compY2',
      vec         = list(c("pg02","ps2"),
                         c("pg12","ps2")),
      Id          = id.var,
      method      = method,
      na_handling = 'drop',
      k           = K,
      replace     = TRUE,
      caliper     = NULL,
      distance    = "mahalanobis",
      compW       = "ipcw.R",
      plotbalance = plotbalance
    )

    gridpar2 <- expand.grid(
      mtry  = unique(c(floor(length(names.var2)),
                       round(seq(1, length(names.var2), length.out = 12)),
                       round(c(.75, 1) * length(names.var2)))),
      ntree = 1000,
      nodesize = c(2, 10)
    ) |>
      dplyr::filter(mtry > 0) |>
      unique()

    obs2 <- data.frame(
      A = subHs2.match$newTxt,
      subHs2.match[, names.var2, drop = FALSE],
      stringsAsFactors = FALSE
    )

    rfresult2 <- rfdtr(
      modeltype = modeltype,
      usecv     = usecv,
      sl.seed   = 123,
      metric    = 'policyval',
      A.obs     = subHs2.match[[A2.var]],
      Q.obs     = subHs2.match$compY2,
      Q.match   = subHs2.match$pairedCompY,
      obs       = obs2,
      W         = subHs2.match$match.weight,
      gridpar   = gridpar2
    )

    subHs2.match$estAopt.s2.ml <- rfresult2$estA.obs
    subHs2.match$estCompOSopt.s2.ml <- subHs2.match$compY2
    flip2 <- subHs2.match[[A2.var]] != subHs2.match$estAopt.s2.ml
    flip2[is.na(flip2)] <- FALSE
    subHs2.match$estCompOSopt.s2.ml[flip2] <- subHs2.match$pairedCompY[flip2]

    subb <- subHs2.match |>
      dplyr::select(
        dplyr::all_of(id.var),
        estAopt.s2.ml,
        estCompOSopt.s2.ml
      )

    # -------------------------
    # Stage 1
    # -------------------------
    subHs1 <- merge(MLdata, subb, by = id.var, all = TRUE)

    subHs1$Rtilde <- subHs1$compOY
    idx_eta2 <- subHs1[[eta2.var]] == 1
    idx_eta2[is.na(idx_eta2)] <- FALSE
    subHs1$Rtilde[idx_eta2] <- subHs1[[Y1.var]][idx_eta2] + subHs1$estCompOSopt.s2.ml[idx_eta2]
    subHs1 <- subHs1[!is.na(subHs1$Rtilde), , drop = FALSE]
    subHs1$Rtilde[subHs1$Rtilde > rtau[1]] <- rtau[1]

    xx <- exact1.vars[!exact1.vars %in% c(A1.var, A2.var)]
    exact1.vars1 <- if (length(xx) == 0) NULL else xx

    subHs1.match <- matchpotential_DTR(
      dat         = subHs1,
      txgroup     = A1.var,
      exact_vars  = exact1.vars1,
      compY       = 'Rtilde',
      vec         = list(c("pg01","ps1"),
                         c("pg11","ps1")),
      Id          = id.var,
      method      = method,
      na_handling = 'drop',
      k           = K,
      replace     = TRUE,
      caliper     = NULL,
      distance    = "mahalanobis",
      compW       = "ipcw.R",
      plotbalance = plotbalance
    )

    gridpar1 <- expand.grid(
      mtry  = unique(c(floor(length(names.var1)),
                       round(seq(1, length(names.var1), length.out = 12)),
                       round(c(.75, .90, 1) * length(names.var1)))),
      ntree = 1000,
      nodesize = c(2, 10)
    ) |>
      dplyr::filter(mtry > 0) |>
      unique()

    obs1 <- data.frame(
      A = subHs1.match$newTxt,
      subHs1.match[, names.var1, drop = FALSE],
      stringsAsFactors = FALSE
    )

    rfresult1 <- rfdtr(
      modeltype = modeltype,
      usecv     = usecv,
      sl.seed   = 123,
      metric    = 'policyval',
      A.obs     = subHs1.match[[A1.var]],
      Q.obs     = subHs1.match$Rtilde,
      Q.match   = subHs1.match$pairedCompY,
      obs       = obs1,
      W         = subHs1.match$match.weight,
      gridpar   = gridpar1
    )

    # -------------------------
    # Apparent fitted actions on training data
    # -------------------------
    A1.opt.train <- as.numeric(as.character(
      predict(rfresult1$model, data = MLdata[, names.var1, drop = FALSE])$predictions
    ))

    A2.opt.train <- rep(NA_real_, nrow(MLdata))
    ok2 <- MLdata[[eta2.var]] == 1 & complete.cases(MLdata[, names.var2, drop = FALSE])
    ok2[is.na(ok2)] <- FALSE
    if (any(ok2)) {
      A2.opt.train[ok2] <- as.numeric(as.character(
        predict(rfresult2$model, data = MLdata[ok2, names.var2, drop = FALSE])$predictions
      ))
    }

    mod = structure(
      list(
        stage1_model = rfresult1$model,
        stage2_model = rfresult2$model,
        stage1_fit   = rfresult1,
        stage2_fit   = rfresult2,
        names.var1   = names.var1,
        names.var2   = names.var2,
        id.var       = id.var,
        eta2.var     = eta2.var,
        A1.var       = A1.var,
        A2.var       = A2.var,
        cap_months   = cap_months,
        modeltype    = modeltype,
        train_actions = data.frame(
          id     = MLdata[[id.var]],
          eta2   = MLdata[[eta2.var]],
          A1.opt = A1.opt.train,
          A2.opt = A2.opt.train
        ),
        tuning = list(
          stage1 = rfresult1$best,
          stage2 = rfresult2$best
        ),
        data = MLdata,
        call = match.call()
      ),
      class = c("Drmatch", "list")
    )

  }, error = function(e) {
    message("Error in deployment fit: ", conditionMessage(e))
    NULL
  })

  fit
}

