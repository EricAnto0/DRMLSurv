
#' Predict optimal treatment decisions from a fitted Drmatch object
#'
#' Generates predicted optimal treatment assignments from a fitted
#' `Drmatch` object for new observations using stage-specific random
#' forest models fitted with the \pkg{ranger} package. Predictions may be
#' obtained for stage 1 only, stage 2 only, or both stages of a two-stage
#' treatment regime.
#'
#' For stage 1, predictions are produced for rows in `newdata` with complete
#' values on the stage-1 predictor set stored in `object$names.var1`. For
#' stage 2, predictions are produced only for rows that are eligible for
#' second-stage treatment, as indicated by `object$eta2.var == 1`, and that
#' also have complete values on the stage-2 predictor set stored in
#' `object$names.var2`.
#'
#' Rows that are ineligible for prediction because of missing required
#' covariates, or because they are not eligible for stage 2, receive `NA`
#' for the corresponding predicted treatment.
#'
#' @param object A fitted `Drmatch` object containing the trained
#'   stage-specific \pkg{ranger} random forest models and metadata needed for
#'   deployment, including `stage1_model`, `stage2_model`, `names.var1`,
#'   `names.var2`, and `eta2.var`.
#' @param newdata A data frame containing the predictor variables required for
#'   stage-1 and/or stage-2 prediction.
#' @param stage Character string indicating which stage predictions to return.
#'   Must be one of `"both"`, `"stage1"`, or `"stage2"`. The default is
#'   `"both"`.
#' @param ... Additional arguments passed through for S3 compatibility.
#'
#' @details
#' This function is an S3 `predict()` method for objects of class
#' `Drmatch`.
#'
#' When `stage = "stage1"`, the function checks that all variables listed in
#' `object$names.var1` are present in `newdata`. Predictions are returned in
#' the column `A1.opt`.
#'
#' When `stage = "stage2"`, the function checks that all variables listed in
#' `object$names.var2` are present in `newdata`, and also verifies that the
#' stage-2 eligibility indicator named by `object$eta2.var` is available in
#' `newdata`. Predictions are returned in the column `A2.opt`.
#'
#' When `stage = "both"`, both `A1.opt` and `A2.opt` are returned.
#'
#' Internally, stage-specific predictions are obtained from fitted
#' \pkg{ranger} models using the `predictions` component returned by
#' [`predict.ranger()`]:
#' \preformatted{
#' predict(model, data = newdata_subset)$predictions
#' }
#' For classification forests, these predicted class labels are coerced to
#' numeric values, so the fitted models are expected to predict treatment
#' classes coded as `-1` and `1`.
#'
#' @return
#' A data frame with one row per row of `newdata`. The returned data frame
#' always contains:
#'
#' \describe{
#'   \item{`row_id`}{Row index corresponding to the original row position in `newdata`.}
#' }
#'
#' Depending on `stage`, it also contains:
#'
#' \describe{
#'   \item{`A1.opt`}{Predicted optimal stage-1 treatment, or `NA` when stage-1 prediction is not available for that row.}
#'   \item{`A2.opt`}{Predicted optimal stage-2 treatment, or `NA` when stage-2 prediction is not available for that row.}
#' }
#'
#' @section Required components of `object`:
#' The fitted `Drmatch` object is expected to contain at least the
#' following elements:
#' \describe{
#'   \item{`stage1_model`}{A fitted \pkg{ranger} classification model for stage 1.}
#'   \item{`stage2_model`}{A fitted \pkg{ranger} classification model for stage 2.}
#'   \item{`names.var1`}{Character vector of predictor names required for stage 1.}
#'   \item{`names.var2`}{Character vector of predictor names required for stage 2.}
#'   \item{`eta2.var`}{Name of the stage-2 eligibility indicator in `newdata`.}
#' }
#'
#' @examples
#' ## Not run:
#' ## Suppose `fit` is a fitted Drmatch object
#' ## and `new_patients` is a data frame of candidate patients.
#' ##
#' ## Predict both stages
#' ## pred <- predict(fit, newdata = new_patients, stage = "both")
#' ##
#' ## Predict stage 1 only
#' ## pred1 <- predict(fit, newdata = new_patients, stage = "stage1")
#' ##
#' ## Predict stage 2 only
#' ## pred2 <- predict(fit, newdata = new_patients, stage = "stage2")
#'
#' @method predict Drmatch
#' @export



predict.Drmatch <- function(object, newdata, stage = c("both", "stage1", "stage2"), ...) {
  stage <- match.arg(stage)

  out <- data.frame(row_id = seq_len(nrow(newdata)))

  if (stage %in% c("both", "stage1")) {
    miss1 <- setdiff(object$names.var1, names(newdata))
    if (length(miss1) > 0) {
      stop("Missing stage-1 variables in `newdata`: ", paste(miss1, collapse = ", "))
    }

    ok1 <- complete.cases(newdata[, object$names.var1, drop = FALSE])
    A1.opt <- rep(NA_real_, nrow(newdata))

    if (any(ok1)) {
      A1.opt[ok1] <- as.numeric(as.character(
      #ranger:::predict.ranger(object$stage1_model, data = newdata[ok1, object$names.var1, drop = FALSE])$predictions
      predict(object$stage1_model, data = newdata[ok1, object$names.var1, drop = FALSE])$predictions
      ))
    }

    out$A1.opt <- A1.opt
  }

  if (stage %in% c("both", "stage2")) {
    miss2 <- setdiff(object$names.var2, names(newdata))
    if (length(miss2) > 0) {
      stop("Missing stage-2 variables in `newdata`: ", paste(miss2, collapse = ", "))
    }

    if (!object$eta2.var %in% names(newdata)) {
      stop("`newdata` must contain ", object$eta2.var, " to determine who is eligible for stage 2 prediction.")
    }

    ok2 <- newdata[[object$eta2.var]] == 1 & complete.cases(newdata[, object$names.var2, drop = FALSE])
    ok2[is.na(ok2)] <- FALSE

    A2.opt <- rep(NA_real_, nrow(newdata))
    if (any(ok2)) {
      A2.opt[ok2] <- as.numeric(as.character(
        #ranger:::predict.ranger(object$stage2_model, data = newdata[ok2, object$names.var2, drop = FALSE])$predictions
        predict(object$stage2_model, data = newdata[ok2, object$names.var2, drop = FALSE])$predictions
      ))
    }

    out$A2.opt <- A2.opt
  }

  out
}



#' Summarize policy performance metrics for a two-stage treatment regime
#'
#' Computes classification and optional value-based performance metrics for
#' estimated treatment decisions in a two-stage dynamic treatment regime.
#' The function compares estimated stage-1 and stage-2 treatment assignments
#' against observed assignments and returns stage-specific as well as joint
#' treatment-path metrics.
#'
#' Stage-1 and stage-2 summaries include accuracy, Matthews correlation
#' coefficient (MCC), sensitivity, specificity, positive predictive value
#' (PPV), negative predictive value (NPV), and F1 score. For the joint
#' two-stage treatment path, the function additionally returns multiclass MCC
#' together with macro-averaged and prevalence-weighted classification
#' summaries across the four possible treatment paths.
#'
#' Optional RMST-style value summaries can also be computed when corresponding
#' outcome columns are supplied.
#'
#' @param tmpData A data frame containing observed treatment assignments,
#'   stage-2 eligibility indicators, and optionally outcome/value columns.
#' @param estA1 A vector of estimated stage-1 treatment assignments. Must have
#'   length equal to `nrow(tmpData)`. Expected coding is `-1` and `1`.
#' @param estA2 A vector of estimated stage-2 treatment assignments. Must have
#'   length equal to `nrow(tmpData)`. Expected coding is `-1` and `1`.
#' @param obs1_var Character string giving the column name in `tmpData`
#'   containing observed stage-1 treatment assignments. Default is `"A1.var"`.
#' @param obs2_var Character string giving the column name in `tmpData`
#'   containing observed stage-2 treatment assignments. Default is `"A2.var"`.
#' @param eta2_var Character string giving the column name in `tmpData`
#'   indicating whether an individual proceeds to stage 2. Subjects with
#'   `tmpData[[eta2_var]] == 1` are included in stage-2 and joint-path metrics.
#'   Default is `"eta2"`.
#' @param Tc1_var Optional character string giving the column name for the
#'   stage-1 value under control or treatment option `-1`. Default is `NULL`.
#' @param Tt1_var Optional character string giving the column name for the
#'   stage-1 value under treatment option `1`. Default is `NULL`.
#' @param Tc2_var Optional character string giving the column name for the
#'   stage-2 value under control or treatment option `-1`. Default is `NULL`.
#' @param Tt2_var Optional character string giving the column name for the
#'   stage-2 value under treatment option `1`. Default is `NULL`.
#' @param Tc_total_var Optional character string giving the column name for the
#'   total outcome/value under control or treatment option `-1`. Default is
#'   `NULL`.
#' @param Tt_total_var Optional character string giving the column name for the
#'   total outcome/value under treatment option `1`. Default is `NULL`.
#' @param stage Character string indicating which stage(s) to summarize. Must be one of `"both"`, `"stage1"`, or `"stage2"`. Default is `"both"`.
#' @details
#' The function assumes binary treatment coding with values `-1` and `1`.
#' Internally, estimated and observed treatments are converted to factors with
#' levels `c(-1, 1)`.
#'
#' Stage-1 metrics are computed using all complete cases for the estimated and
#' observed stage-1 treatment assignments. Stage-2 and joint metrics are
#' computed only among subjects satisfying `eta2 == 1`.
#'
#' Joint treatment paths are defined as the four-level factor
#' `c("-1_-1", "-1_1", "1_-1", "1_1")`, corresponding to the sequence of
#' stage-1 and stage-2 treatment assignments.
#'
#' Macro-averaged metrics are calculated as the unweighted mean of class-wise
#' statistics from the joint confusion matrix. Weighted summaries use empirical
#' class prevalences from the observed joint treatment paths.
#'
#' Optional RMST-style summaries are computed when the corresponding column
#' names are supplied and found in `tmpData`. These use an internal scoring
#' rule that assigns subject-level values based on the estimated treatment.
#'
#' @return A named numeric vector containing some or all of the following
#'   components:
#'
#' \describe{
#'   \item{`Acc1L`}{Stage-1 accuracy.}
#'   \item{`MCC1L`}{Stage-1 Matthews correlation coefficient.}
#'   \item{`Sensitivity1L`}{Stage-1 sensitivity.}
#'   \item{`Specificity1L`}{Stage-1 specificity.}
#'   \item{`PPV1L`}{Stage-1 positive predictive value.}
#'   \item{`NPV1L`}{Stage-1 negative predictive value.}
#'   \item{`F11L`}{Stage-1 F1 score.}
#'   \item{`RMST1L`}{Optional stage-1 value summary.}
#'
#'   \item{`Acc2L`}{Stage-2 accuracy among subjects with `eta2 == 1`.}
#'   \item{`MCC2L`}{Stage-2 Matthews correlation coefficient.}
#'   \item{`Sensitivity2L`}{Stage-2 sensitivity.}
#'   \item{`Specificity2L`}{Stage-2 specificity.}
#'   \item{`PPV2L`}{Stage-2 positive predictive value.}
#'   \item{`NPV2L`}{Stage-2 negative predictive value.}
#'   \item{`F12L`}{Stage-2 F1 score.}
#'   \item{`RMST2L`}{Optional stage-2 value summary.}
#'
#'   \item{`AccTotal`}{Joint-path accuracy.}
#'   \item{`MCCTotal`}{Multiclass Matthews correlation coefficient for the joint path.}
#'   \item{`SensitivityTotalMacro`}{Macro-averaged sensitivity across joint classes.}
#'   \item{`SpecificityTotalMacro`}{Macro-averaged specificity across joint classes.}
#'   \item{`PPVTotalMacro`}{Macro-averaged positive predictive value across joint classes.}
#'   \item{`NPVTotalMacro`}{Macro-averaged negative predictive value across joint classes.}
#'   \item{`F1TotalMacro`}{Macro-averaged F1 score across joint classes.}
#'   \item{`SensitivityTotalWtd`}{Prevalence-weighted sensitivity across joint classes.}
#'   \item{`SpecificityTotalWtd`}{Prevalence-weighted specificity across joint classes.}
#'   \item{`PPVTotalWtd`}{Prevalence-weighted positive predictive value across joint classes.}
#'   \item{`NPVTotalWtd`}{Prevalence-weighted negative predictive value across joint classes.}
#'   \item{`F1TotalWtd`}{Prevalence-weighted F1 score across joint classes.}
#'   \item{`RMSTTotal`}{Optional total value summary.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 100
#' dat <- data.frame(
#'   A1.var = sample(c(-1, 1), n, replace = TRUE),
#'   A2.var = sample(c(-1, 1), n, replace = TRUE),
#'   eta2   = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7))
#' )
#'
#' estA1 <- sample(c(-1, 1), n, replace = TRUE)
#' estA2 <- sample(c(-1, 1), n, replace = TRUE)
#'
#' policy_summary_metrics(
#'   tmpData = dat,
#'   estA1 = estA1,
#'   estA2 = estA2
#' )
#'}
#' @seealso [summary.Drmatch()], [print.summary.Drmatch()],
#'   [caret::confusionMatrix()], [mltools::mcc()]
#'
#' @export

policy_summary_metrics <- function(
    tmpData,
    estA1,
    estA2,
    obs1_var      = "A1.var",
    obs2_var      = "A2.var",
    eta2_var      = "eta2",
    Tc1_var       = NULL,
    Tt1_var       = NULL,
    Tc2_var       = NULL,
    Tt2_var       = NULL,
    Tc_total_var  = NULL,
    Tt_total_var  = NULL,
    stage         = c("both", "stage1", "stage2")
) {
  stage <- match.arg(stage)

  has_cols <- function(dat, vars) {
    if (length(vars) == 0L) return(FALSE)
    vals <- unlist(vars, use.names = FALSE)
    if (length(vals) != length(vars)) return(FALSE)
    if (any(is.na(vals)) || any(!nzchar(vals))) return(FALSE)
    all(vals %in% names(dat))
  }

  safe_div <- function(num, den) {
    if (!is.finite(den) || den <= 0) return(NA_real_)
    num / den
  }

  my_score.Surv <- function(pred, Q, pQ) {
    pred <- as.numeric(as.character(pred))
    mean(ifelse(pred == 1, pQ, Q), na.rm = TRUE)
  }

  binary_metrics <- function(actual, pred, negative = -1, positive = 1) {
    actual <- suppressWarnings(as.numeric(as.character(actual)))
    pred   <- suppressWarnings(as.numeric(as.character(pred)))

    keep <- complete.cases(actual, pred) &
      actual %in% c(negative, positive) &
      pred   %in% c(negative, positive)

    if (!any(keep)) {
      return(list(
        acc  = NA_real_,
        mcc  = NA_real_,
        sens = NA_real_,
        spec = NA_real_,
        ppv  = NA_real_,
        npv  = NA_real_,
        f1   = NA_real_
      ))
    }

    actual <- actual[keep]
    pred   <- pred[keep]

    tp <- sum(actual == positive & pred == positive)
    tn <- sum(actual == negative & pred == negative)
    fp <- sum(actual == negative & pred == positive)
    fn <- sum(actual == positive & pred == negative)

    acc  <- safe_div(tp + tn, tp + tn + fp + fn)
    sens <- safe_div(tp, tp + fn)
    spec <- safe_div(tn, tn + fp)
    ppv  <- safe_div(tp, tp + fp)
    npv  <- safe_div(tn, tn + fn)
    f1   <- safe_div(2 * tp, 2 * tp + fp + fn)

    den_mcc <- sqrt(
      as.double(tp + fp) *
        as.double(tp + fn) *
        as.double(tn + fp) *
        as.double(tn + fn)
    )

    mcc <- if (!is.finite(den_mcc) || den_mcc == 0) {
      NA_real_
    } else {
      (tp * tn - fp * fn) / den_mcc
    }

    list(
      acc  = acc,
      mcc  = mcc,
      sens = sens,
      spec = spec,
      ppv  = ppv,
      npv  = npv,
      f1   = f1
    )
  }

  multiclass_mcc <- function(actual, pred, levs = NULL) {
    actual <- as.character(actual)
    pred   <- as.character(pred)

    keep <- complete.cases(actual, pred)
    if (!any(keep)) return(NA_real_)

    actual <- actual[keep]
    pred   <- pred[keep]

    if (is.null(levs)) {
      levs <- sort(unique(c(actual, pred)))
    }

    actual <- factor(actual, levels = levs)
    pred   <- factor(pred,   levels = levs)

    C <- as.matrix(table(actual, pred))
    C <- apply(C, c(1, 2), as.double)

    s <- sum(C)
    if (s == 0) return(NA_real_)

    c0  <- sum(diag(C))
    tk  <- rowSums(C)
    pk  <- colSums(C)
    den <- sqrt((s^2 - sum(pk^2)) * (s^2 - sum(tk^2)))

    if (!is.finite(den) || den == 0) return(NA_real_)
    (c0 * s - sum(tk * pk)) / den
  }

  multiclass_one_vs_rest <- function(actual, pred, levs) {
    actual <- as.character(actual)
    pred   <- as.character(pred)

    keep <- complete.cases(actual, pred)
    if (!any(keep)) {
      return(list(
        macro = list(sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_, f1 = NA_real_),
        weighted = list(sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_, f1 = NA_real_)
      ))
    }

    actual <- factor(actual[keep], levels = levs)
    pred   <- factor(pred[keep],   levels = levs)

    C <- as.matrix(table(actual, pred))
    total_n <- sum(C)

    if (total_n == 0) {
      return(list(
        macro = list(sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_, f1 = NA_real_),
        weighted = list(sens = NA_real_, spec = NA_real_, ppv = NA_real_, npv = NA_real_, f1 = NA_real_)
      ))
    }

    by_class <- lapply(seq_along(levs), function(j) {
      tp <- C[j, j]
      fn <- sum(C[j, ]) - tp
      fp <- sum(C[, j]) - tp
      tn <- total_n - tp - fn - fp

      list(
        sens = safe_div(tp, tp + fn),
        spec = safe_div(tn, tn + fp),
        ppv  = safe_div(tp, tp + fp),
        npv  = safe_div(tn, tn + fn),
        f1   = safe_div(2 * tp, 2 * tp + fp + fn)
      )
    })

    wts <- rowSums(C) / total_n

    macro <- list(
      sens = mean(vapply(by_class, `[[`, numeric(1), "sens"), na.rm = TRUE),
      spec = mean(vapply(by_class, `[[`, numeric(1), "spec"), na.rm = TRUE),
      ppv  = mean(vapply(by_class, `[[`, numeric(1), "ppv"),  na.rm = TRUE),
      npv  = mean(vapply(by_class, `[[`, numeric(1), "npv"),  na.rm = TRUE),
      f1   = mean(vapply(by_class, `[[`, numeric(1), "f1"),   na.rm = TRUE)
    )

    weighted <- list(
      sens = sum(vapply(by_class, `[[`, numeric(1), "sens") * wts, na.rm = TRUE),
      spec = sum(vapply(by_class, `[[`, numeric(1), "spec") * wts, na.rm = TRUE),
      ppv  = sum(vapply(by_class, `[[`, numeric(1), "ppv")  * wts, na.rm = TRUE),
      npv  = sum(vapply(by_class, `[[`, numeric(1), "npv")  * wts, na.rm = TRUE),
      f1   = sum(vapply(by_class, `[[`, numeric(1), "f1")   * wts, na.rm = TRUE)
    )

    list(macro = macro, weighted = weighted)
  }

  if (!is.data.frame(tmpData)) {
    stop("`tmpData` must be a data.frame.", call. = FALSE)
  }
  if (!obs1_var %in% names(tmpData)) {
    stop("`obs1_var` not found in `tmpData`.", call. = FALSE)
  }
  if (!obs2_var %in% names(tmpData)) {
    stop("`obs2_var` not found in `tmpData`.", call. = FALSE)
  }
  if (!eta2_var %in% names(tmpData)) {
    stop("`eta2_var` not found in `tmpData`.", call. = FALSE)
  }

  obs1 <- tmpData[[obs1_var]]
  obs2 <- tmpData[[obs2_var]]
  eta2 <- tmpData[[eta2_var]]

  estA1 <- suppressWarnings(as.numeric(as.character(estA1)))
  estA2 <- suppressWarnings(as.numeric(as.character(estA2)))

  if (length(estA1) != nrow(tmpData)) {
    stop("`estA1` must have length nrow(tmpData).", call. = FALSE)
  }
  if (length(estA2) != nrow(tmpData)) {
    stop("`estA2` must have length nrow(tmpData).", call. = FALSE)
  }

  out <- list(
    Acc1L = NA_real_,
    MCC1L = NA_real_,
    Sensitivity1L = NA_real_,
    Specificity1L = NA_real_,
    PPV1L = NA_real_,
    NPV1L = NA_real_,
    F11L = NA_real_,
    RMST1L = NA_real_,

    Acc2L = NA_real_,
    MCC2L = NA_real_,
    Sensitivity2L = NA_real_,
    Specificity2L = NA_real_,
    PPV2L = NA_real_,
    NPV2L = NA_real_,
    F12L = NA_real_,
    RMST2L = NA_real_,

    AccTotal = NA_real_,
    MCCTotal = NA_real_,
    SensitivityTotalMacro = NA_real_,
    SpecificityTotalMacro = NA_real_,
    PPVTotalMacro = NA_real_,
    NPVTotalMacro = NA_real_,
    F1TotalMacro = NA_real_,
    SensitivityTotalWtd = NA_real_,
    SpecificityTotalWtd = NA_real_,
    PPVTotalWtd = NA_real_,
    NPVTotalWtd = NA_real_,
    F1TotalWtd = NA_real_,
    RMSTTotal = NA_real_
  )

  if (stage %in% c("both", "stage1")) {
    m1 <- binary_metrics(actual = obs1, pred = estA1)

    out$Acc1L         <- m1$acc
    out$MCC1L         <- m1$mcc
    out$Sensitivity1L <- m1$sens
    out$Specificity1L <- m1$spec
    out$PPV1L         <- m1$ppv
    out$NPV1L         <- m1$npv
    out$F11L          <- m1$f1

    if (has_cols(tmpData, list(Tc1_var, Tt1_var))) {
      out$RMST1L <- my_score.Surv(
        pred = estA1,
        Q    = tmpData[[Tc1_var]],
        pQ   = tmpData[[Tt1_var]]
      )
    }
  }

  if (stage %in% c("both", "stage2")) {
    idx2 <- which(eta2 == 1)

    if (length(idx2) > 0) {
      m2 <- binary_metrics(actual = obs2[idx2], pred = estA2[idx2])

      out$Acc2L         <- m2$acc
      out$MCC2L         <- m2$mcc
      out$Sensitivity2L <- m2$sens
      out$Specificity2L <- m2$spec
      out$PPV2L         <- m2$ppv
      out$NPV2L         <- m2$npv
      out$F12L          <- m2$f1

      if (has_cols(tmpData, list(Tc2_var, Tt2_var))) {
        out$RMST2L <- my_score.Surv(
          pred = estA2[idx2],
          Q    = tmpData[[Tc2_var]][idx2],
          pQ   = tmpData[[Tt2_var]][idx2]
        )
      }
    }
  }

  if (stage == "both") {
    idx_joint <- which(eta2 == 1)
    joint_levels <- c("-1_-1", "-1_1", "1_-1", "1_1")

    if (length(idx_joint) > 0) {
      pred_joint <- paste(estA1[idx_joint], estA2[idx_joint], sep = "_")
      obs_joint  <- paste(obs1[idx_joint],  obs2[idx_joint],  sep = "_")

      keep_joint <- complete.cases(pred_joint, obs_joint) &
        pred_joint %in% joint_levels &
        obs_joint  %in% joint_levels

      if (any(keep_joint)) {
        pred_joint <- pred_joint[keep_joint]
        obs_joint  <- obs_joint[keep_joint]

        out$AccTotal <- mean(pred_joint == obs_joint, na.rm = TRUE)
        out$MCCTotal <- multiclass_mcc(obs_joint, pred_joint, levs = joint_levels)

        joint_stats <- multiclass_one_vs_rest(
          actual = obs_joint,
          pred   = pred_joint,
          levs   = joint_levels
        )

        out$SensitivityTotalMacro <- joint_stats$macro$sens
        out$SpecificityTotalMacro <- joint_stats$macro$spec
        out$PPVTotalMacro         <- joint_stats$macro$ppv
        out$NPVTotalMacro         <- joint_stats$macro$npv
        out$F1TotalMacro          <- joint_stats$macro$f1

        out$SensitivityTotalWtd   <- joint_stats$weighted$sens
        out$SpecificityTotalWtd   <- joint_stats$weighted$spec
        out$PPVTotalWtd           <- joint_stats$weighted$ppv
        out$NPVTotalWtd           <- joint_stats$weighted$npv
        out$F1TotalWtd            <- joint_stats$weighted$f1
      }
    }

    if (has_cols(tmpData, list(Tc_total_var, Tt_total_var))) {
      out$RMSTTotal <- my_score.Surv(
        pred = estA1,
        Q    = tmpData[[Tc_total_var]],
        pQ   = tmpData[[Tt_total_var]]
      )
    }
  }

  out
}

#' Summarize a fitted `Drmatch` object on new data
#'
#' Generates predicted treatment decisions from a fitted `Drmatch` object and
#' summarizes their performance against observed treatment assignments in a new
#' dataset. Depending on the requested stage, the function evaluates stage-1
#' decisions, stage-2 decisions, or both jointly.
#'
#' Internally, predicted treatment assignments are passed to
#' [`policy_summary_metrics()`] to compute stage-specific and joint-path
#' classification metrics, with optional value-based summaries when outcome
#' columns are supplied.
#'
#' @param object A fitted object of class `"Drmatch"`.
#' @param newdata A data frame on which to compute predicted treatment
#'   assignments and performance summaries.
#' @param stage Character string indicating which decision stage to summarize.
#'   One of `"both"`, `"stage1"`, or `"stage2"`. Default is `"both"`.
#' @param obs1_var Optional character string giving the observed stage-1
#'   treatment column in `newdata`. If `NULL`, the value is taken from
#'   `object$A1.var`.
#' @param obs2_var Optional character string giving the observed stage-2
#'   treatment column in `newdata`. If `NULL`, the value is taken from
#'   `object$A2.var`.
#' @param eta2_var Optional character string giving the stage-2 eligibility
#'   column in `newdata`. If `NULL`, the value is taken from `object$eta2.var`.
#' @param pred1_var Character string giving the column name in the prediction
#'   output corresponding to the estimated stage-1 treatment. Default is
#'   `"A1.opt"`.
#' @param pred2_var Character string giving the column name in the prediction
#'   output corresponding to the estimated stage-2 treatment. Default is
#'   `"A2.opt"`.
#' @param Tc1_var Optional character string giving the stage-1 value column for
#'   treatment option `-1`.
#' @param Tt1_var Optional character string giving the stage-1 value column for
#'   treatment option `1`.
#' @param Tc2_var Optional character string giving the stage-2 value column for
#'   treatment option `-1`.
#' @param Tt2_var Optional character string giving the stage-2 value column for
#'   treatment option `1`.
#' @param Tc_total_var Optional character string giving the total value column
#'   for treatment option `-1`.
#' @param Tt_total_var Optional character string giving the total value column
#'   for treatment option `1`.
#' @param as.data.frame Logical; if `TRUE`, the output is returned as a
#'   one-row data frame. Otherwise, a named numeric vector is returned.
#'   Default is `TRUE`.
#' @param overall_type Character string indicating the type of overall summary metrics to compute for the joint path. One of `"macro"` or `"weighted"`. Default is `"macro"`.
#' @param ... Additional arguments passed to [predict.Drmatch()].
#'
#' @details
#' When `stage = "stage1"`, stage-2 predictions are set to `NA` before summary
#' metrics are computed. When `stage = "stage2"`, the observed stage-1
#' treatment from `newdata` is used as the stage-1 component of the joint path,
#' while stage-2 predictions are taken from `predict(object, ...)`.
#'
#' This method requires `newdata` because performance is evaluated by comparing
#' predicted treatment decisions with observed treatments and, optionally, with
#' value columns provided in the new dataset.
#'
#' @return Either a one-row data frame or a named numeric vector containing the
#'   policy performance metrics returned by [policy_summary_metrics()].
#'
#' @examples
#' \dontrun{
#' # fit <- Drmatch(...)
#' # summary(fit, newdata = test_dat)
#' # summary(fit, newdata = test_dat, stage = "stage1")
#' # summary(fit, newdata = test_dat, stage = "stage2", as.data.frame = FALSE)
#' }
#'
#' @seealso [policy_summary_metrics()], [print.summary.Drmatch()],
#'   [predict.Drmatch()]
#'
#' @method summary Drmatch
#' @export
summary.Drmatch <- function(
    object,
    newdata,
    stage = c("both", "stage1", "stage2"),
    overall_type = c("macro", "weighted"),
    obs1_var = NULL,
    obs2_var = NULL,
    eta2_var = NULL,
    pred1_var = "A1.opt",
    pred2_var = "A2.opt",
    Tc1_var = NULL,
    Tt1_var = NULL,
    Tc2_var = NULL,
    Tt2_var = NULL,
    Tc_total_var = NULL,
    Tt_total_var = NULL,
    as.data.frame = TRUE,
    ...
) {
  stage <- match.arg(stage)
  overall_type <- match.arg(overall_type)

  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` must be supplied to `summary.Drmatch()`.", call. = FALSE)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data.frame.", call. = FALSE)
  }

  if (is.null(obs1_var)) obs1_var <- object$A1.var
  if (is.null(obs2_var)) obs2_var <- object$A2.var
  if (is.null(eta2_var)) eta2_var <- object$eta2.var

  if (is.null(obs1_var) || length(obs1_var) != 1L || !nzchar(obs1_var) || !obs1_var %in% names(newdata)) {
    stop("Could not resolve `obs1_var` in `newdata`.", call. = FALSE)
  }
  if (is.null(obs2_var) || length(obs2_var) != 1L || !nzchar(obs2_var) || !obs2_var %in% names(newdata)) {
    stop("Could not resolve `obs2_var` in `newdata`.", call. = FALSE)
  }
  if (is.null(eta2_var) || length(eta2_var) != 1L || !nzchar(eta2_var) || !eta2_var %in% names(newdata)) {
    stop("Could not resolve `eta2_var` in `newdata`.", call. = FALSE)
  }

  pred <- predict(object, newdata = newdata, stage = stage, ...)

  if (stage == "both") {
    if (!pred1_var %in% names(pred)) {
      stop("`pred1_var` not found in predictions for stage = 'both'.", call. = FALSE)
    }
    if (!pred2_var %in% names(pred)) {
      stop("`pred2_var` not found in predictions for stage = 'both'.", call. = FALSE)
    }
  }

  if (stage == "stage1") {
    if (!pred1_var %in% names(pred)) {
      stop("`pred1_var` not found in predictions for stage = 'stage1'.", call. = FALSE)
    }
    pred[[pred2_var]] <- rep(NA_real_, nrow(newdata))
  }

  if (stage == "stage2") {
    if (!pred2_var %in% names(pred)) {
      stop("`pred2_var` not found in predictions for stage = 'stage2'.", call. = FALSE)
    }
    pred[[pred1_var]] <- rep(NA_real_, nrow(newdata))
  }

  out <- policy_summary_metrics(
    tmpData       = newdata,
    estA1         = pred[[pred1_var]],
    estA2         = pred[[pred2_var]],
    obs1_var      = obs1_var,
    obs2_var      = obs2_var,
    eta2_var      = eta2_var,
    Tc1_var       = Tc1_var,
    Tt1_var       = Tt1_var,
    Tc2_var       = Tc2_var,
    Tt2_var       = Tt2_var,
    Tc_total_var  = Tc_total_var,
    Tt_total_var  = Tt_total_var,
    stage         = stage
  )

  metric_rows <- c(
    "Accuracy", "MCC", "Sensitivity", "Specificity",
    "PPV", "NPV", "F1", "RMST"
  )

  if (stage == "stage1") {
    res <- data.frame(
      `Stage 1` = c(
        out$Acc1L,
        out$MCC1L,
        out$Sensitivity1L,
        out$Specificity1L,
        out$PPV1L,
        out$NPV1L,
        out$F11L,
        out$RMST1L
      ),
      row.names = metric_rows,
      check.names = FALSE
    )
  } else if (stage == "stage2") {
    res <- data.frame(
      `Stage 2` = c(
        out$Acc2L,
        out$MCC2L,
        out$Sensitivity2L,
        out$Specificity2L,
        out$PPV2L,
        out$NPV2L,
        out$F12L,
        out$RMST2L
      ),
      row.names = metric_rows,
      check.names = FALSE
    )
  } else {
    overall_vals <- if (overall_type == "weighted") {
      c(
        out$AccTotal,
        out$MCCTotal,
        out$SensitivityTotalWtd,
        out$SpecificityTotalWtd,
        out$PPVTotalWtd,
        out$NPVTotalWtd,
        out$F1TotalWtd,
        out$RMSTTotal
      )
    } else {
      c(
        out$AccTotal,
        out$MCCTotal,
        out$SensitivityTotalMacro,
        out$SpecificityTotalMacro,
        out$PPVTotalMacro,
        out$NPVTotalMacro,
        out$F1TotalMacro,
        out$RMSTTotal
      )
    }

    res <- data.frame(
      `Stage 1` = c(
        out$Acc1L,
        out$MCC1L,
        out$Sensitivity1L,
        out$Specificity1L,
        out$PPV1L,
        out$NPV1L,
        out$F11L,
        out$RMST1L
      ),
      `Stage 2` = c(
        out$Acc2L,
        out$MCC2L,
        out$Sensitivity2L,
        out$Specificity2L,
        out$PPV2L,
        out$NPV2L,
        out$F12L,
        out$RMST2L
      ),
      Overall = overall_vals,
      row.names = metric_rows,
      check.names = FALSE
    )
  }
  if (!is.null(res) && !is.null(dim(res)) && length(dim(res)) == 2L && ncol(res) > 0L) {
    res <- res[!is.na(res[, 1]), , drop = FALSE]
  }

  return(res)
}
#' Print a summary.Drmatch object
#'
#' @param x An object returned by \code{summary.Drmatch()}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.summary.Drmatch <- function(x, ...) {
  print(x, ...)
  invisible(x)
}
