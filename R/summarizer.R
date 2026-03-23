
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
#'
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
    Tt_total_var  = NULL
) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required.")
  }
  if (!requireNamespace("mltools", quietly = TRUE)) {
    stop("Package 'mltools' is required.")
  }

  my_score.Surv <- function(pred, Q, pQ) {
    pred <- as.numeric(as.character(pred))
    mean(ifelse(pred == 1, pQ, Q), na.rm = TRUE)
  }

  multiclass_mcc <- function(actual, pred) {
    actual <- factor(actual)
    pred   <- factor(pred, levels = levels(actual))

    C <- as.matrix(table(actual, pred))
    s <- sum(C)
    c <- sum(diag(C))
    t_k <- rowSums(C)
    p_k <- colSums(C)

    numerator <- c * s - sum(t_k * p_k)
    denominator <- sqrt((s^2 - sum(p_k^2)) * (s^2 - sum(t_k^2)))

    if (denominator == 0) return(NA_real_)
    numerator / denominator
  }

  safe_byclass_extract <- function(cm, stat) {
    out <- cm$byClass[stat]
    if (length(out) == 0 || all(is.na(out))) return(NA_real_)
    unname(out[[1]])
  }

  has_cols <- function(dat, vars) {
    !any(vapply(vars, is.null, logical(1))) && all(vars %in% names(dat))
  }

  if (!obs1_var %in% names(tmpData)) stop("`obs1_var` not found in `tmpData`.")
  if (!obs2_var %in% names(tmpData)) stop("`obs2_var` not found in `tmpData`.")
  if (!eta2_var %in% names(tmpData)) stop("`eta2_var` not found in `tmpData`.")

  obs1 <- tmpData[[obs1_var]]
  obs2 <- tmpData[[obs2_var]]
  eta2 <- tmpData[[eta2_var]]

  estA1 <- as.numeric(as.character(estA1))
  estA2 <- as.numeric(as.character(estA2))

  if (length(estA1) != nrow(tmpData)) {
    stop("`estA1` must have length nrow(tmpData).")
  }
  if (length(estA2) != nrow(tmpData)) {
    stop("`estA2` must have length nrow(tmpData).")
  }

  idx2_use      <- which(eta2 == 1)
  idx_joint_use <- which(eta2 == 1)

  ## Stage 1
  pred1 <- factor(estA1, levels = c(-1, 1))
  obs1f <- factor(obs1,  levels = c(-1, 1))
  keep1 <- complete.cases(pred1, obs1f)

  cm1 <- caret::confusionMatrix(
    data      = pred1[keep1],
    reference = obs1f[keep1],
    positive  = "1"
  )

  rmst1 <- NA_real_
  if (has_cols(tmpData, c(Tc1_var, Tt1_var))) {
    rmst1 <- my_score.Surv(
      pred = pred1,
      Q    = tmpData[[Tc1_var]],
      pQ   = tmpData[[Tt1_var]]
    )
  }

  ## Stage 2
  acc2  <- NA_real_
  mcc2  <- NA_real_
  sens2 <- NA_real_
  spec2 <- NA_real_
  ppv2  <- NA_real_
  npv2  <- NA_real_
  f12   <- NA_real_
  rmst2 <- NA_real_

  if (length(idx2_use) > 0) {
    pred2 <- factor(estA2[idx2_use], levels = c(-1, 1))
    obs2f <- factor(obs2[idx2_use],  levels = c(-1, 1))
    keep2 <- complete.cases(pred2, obs2f)

    if (sum(keep2) > 0) {
      cm2 <- caret::confusionMatrix(
        data      = pred2[keep2],
        reference = obs2f[keep2],
        positive  = "1"
      )

      acc2  <- mean(pred2[keep2] == obs2f[keep2], na.rm = TRUE)
      mcc2  <- mltools::mcc(preds = pred2[keep2], actuals = obs2f[keep2])
      sens2 <- safe_byclass_extract(cm2, "Sensitivity")
      spec2 <- safe_byclass_extract(cm2, "Specificity")
      ppv2  <- safe_byclass_extract(cm2, "Pos Pred Value")
      npv2  <- safe_byclass_extract(cm2, "Neg Pred Value")
      f12   <- safe_byclass_extract(cm2, "F1")

      if (has_cols(tmpData, c(Tc2_var, Tt2_var))) {
        rmst2 <- my_score.Surv(
          pred = pred2,
          Q    = tmpData[[Tc2_var]][idx2_use],
          pQ   = tmpData[[Tt2_var]][idx2_use]
        )
      }
    }
  }

  ## Joint
  acc_total <- NA_real_
  mcc_total <- NA_real_
  sens_total_macro <- NA_real_
  spec_total_macro <- NA_real_
  ppv_total_macro  <- NA_real_
  npv_total_macro  <- NA_real_
  f1_total_macro   <- NA_real_
  sens_total_weighted <- NA_real_
  spec_total_weighted <- NA_real_
  ppv_total_weighted  <- NA_real_
  npv_total_weighted  <- NA_real_
  f1_total_weighted   <- NA_real_

  if (length(idx_joint_use) > 0) {
    pred_joint <- factor(
      paste(estA1[idx_joint_use], estA2[idx_joint_use], sep = "_"),
      levels = c("-1_-1", "-1_1", "1_-1", "1_1")
    )
    obs_joint <- factor(
      paste(obs1[idx_joint_use], obs2[idx_joint_use], sep = "_"),
      levels = c("-1_-1", "-1_1", "1_-1", "1_1")
    )

    keep_joint <- complete.cases(pred_joint, obs_joint)

    if (sum(keep_joint) > 0) {
      pred_joint <- pred_joint[keep_joint]
      obs_joint  <- obs_joint[keep_joint]

      cmtot <- caret::confusionMatrix(pred_joint, obs_joint)
      acc_total <- mean(pred_joint == obs_joint, na.rm = TRUE)
      mcc_total <- multiclass_mcc(actual = obs_joint, pred = pred_joint)

      byclass_total <- cmtot$byClass
      if (is.null(dim(byclass_total))) {
        byclass_total <- matrix(byclass_total, nrow = 1)
        colnames(byclass_total) <- names(cmtot$byClass)
        rownames(byclass_total) <- paste0("Class: ", levels(obs_joint)[1])
      }

      sens_total_macro <- mean(byclass_total[, "Sensitivity"],    na.rm = TRUE)
      spec_total_macro <- mean(byclass_total[, "Specificity"],    na.rm = TRUE)
      ppv_total_macro  <- mean(byclass_total[, "Pos Pred Value"], na.rm = TRUE)
      npv_total_macro  <- mean(byclass_total[, "Neg Pred Value"], na.rm = TRUE)
      f1_total_macro   <- mean(byclass_total[, "F1"],             na.rm = TRUE)

      joint_wts <- prop.table(table(obs_joint))
      joint_names <- sub("^Class: ", "", rownames(byclass_total))
      joint_wts <- as.numeric(joint_wts[joint_names])

      sens_total_weighted <- sum(byclass_total[, "Sensitivity"]    * joint_wts, na.rm = TRUE)
      spec_total_weighted <- sum(byclass_total[, "Specificity"]    * joint_wts, na.rm = TRUE)
      ppv_total_weighted  <- sum(byclass_total[, "Pos Pred Value"] * joint_wts, na.rm = TRUE)
      npv_total_weighted  <- sum(byclass_total[, "Neg Pred Value"] * joint_wts, na.rm = TRUE)
      f1_total_weighted   <- sum(byclass_total[, "F1"]             * joint_wts, na.rm = TRUE)
    }
  }

  ## Optional total RMST
  rmst_total <- NA_real_
  if (has_cols(tmpData, c(Tc_total_var, Tt_total_var))) {
    rmst_total <- my_score.Surv(
      pred = estA1,
      Q    = tmpData[[Tc_total_var]],
      pQ   = tmpData[[Tt_total_var]]
    )
  }

  c(
    Acc1L         = mean(pred1[keep1] == obs1f[keep1], na.rm = TRUE),
    MCC1L         = mltools::mcc(preds = pred1[keep1], actuals = obs1f[keep1]),
    Sensitivity1L = safe_byclass_extract(cm1, "Sensitivity"),
    Specificity1L = safe_byclass_extract(cm1, "Specificity"),
    PPV1L         = safe_byclass_extract(cm1, "Pos Pred Value"),
    NPV1L         = safe_byclass_extract(cm1, "Neg Pred Value"),
    F11L          = safe_byclass_extract(cm1, "F1"),
    RMST1L        = rmst1,

    Acc2L         = acc2,
    MCC2L         = mcc2,
    Sensitivity2L = sens2,
    Specificity2L = spec2,
    PPV2L         = ppv2,
    NPV2L         = npv2,
    F12L          = f12,
    RMST2L        = rmst2,

    AccTotal              = acc_total,
    MCCTotal              = mcc_total,
    SensitivityTotalMacro = sens_total_macro,
    SpecificityTotalMacro = spec_total_macro,
    PPVTotalMacro         = ppv_total_macro,
    NPVTotalMacro         = npv_total_macro,
    F1TotalMacro          = f1_total_macro,
    SensitivityTotalWtd   = sens_total_weighted,
    SpecificityTotalWtd   = spec_total_weighted,
    PPVTotalWtd           = ppv_total_weighted,
    NPVTotalWtd           = npv_total_weighted,
    F1TotalWtd            = f1_total_weighted,
    RMSTTotal             = rmst_total
  )
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

  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` must be supplied to `summary.Drmatch()`.")
  }

  pred <- predict(object, newdata = newdata, stage = stage, ...)

  if (!pred1_var %in% names(pred)) {
    stop("`pred1_var` not found in predictions.")
  }

  if (is.null(obs1_var)) {
    obs1_var <- object$A1.var
  }
  if (is.null(obs2_var)) {
    obs2_var <- object$A2.var
  }
  if (is.null(eta2_var)) {
    eta2_var <- object$eta2.var
  }

  if (stage == "stage1") {
    pred[[pred2_var]] <- NA_real_
  }

  if (stage == "stage2") {
    if (!pred2_var %in% names(pred)) {
      stop("`pred2_var` not found in predictions.")
    }
    pred[[pred1_var]] <- newdata[[obs1_var]]
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
    Tt_total_var  = Tt_total_var
  )

  if (isTRUE(as.data.frame)) {
    out <- as.data.frame(as.list(out), check.names = FALSE)
  }

  out
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
