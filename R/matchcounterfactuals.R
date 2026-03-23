#' Construct matched counterfactual outcomes for a (single-stage) treatment contrast
#'
#' @description
#' Builds matched “potential outcomes” for each unit under both treatment levels by matching each unit
#' to comparable units in the opposite treatment group, and then imputing the opposite-arm outcome using
#' matched donor outcomes. The function supports either:
#' \itemize{
#'   \item \strong{full matching} (one global match structure), or
#'   \item \strong{nearest-neighbor matching} via a two-pass strategy (ATT then ATC) so that every unit
#'   receives an opposite-arm imputation while respecting anchor definitions.
#' }
#'
#' The returned dataset includes (i) the observed outcome, (ii) an imputed opposite-arm “paired” outcome,
#' and (iii) constructed counterfactual outcomes for each arm for every unit.
#'
#' @details
#' \strong{Treatment coding.} Internally, the treatment indicator \code{txgroup} is recoded to \eqn{\{-1,+1\}},
#' where \code{+1} denotes treated and \code{-1} denotes control. Matching formulas are constructed using
#' \code{I(tx == 1)} or \code{I(tx == -1)} to define anchors in the nearest-neighbor passes.
#'
#' \strong{Outcome and weight variables.}
#' \itemize{
#'   \item \code{compY} is the (possibly composite) observed outcome used to form counterfactuals.
#'   \item \code{compW} is an auxiliary “weight-like” variable (default \code{"ipcw.R"}) used to compute
#'   a match quality/weight proxy. If \code{compW} is not present in the data, it is initialized as
#'   \code{compW <- compY} (so weights default to outcome differences).
#' }
#'
#' \strong{Matching covariates.} The matching covariates are supplied via \code{vec} and used to build:
#' \itemize{
#'   \item \code{fmla_bin}: full matching formula using all covariates in \code{vec}
#'   \item \code{fmla_bin_trt}: nearest-neighbor pass 1 formula (treated anchors)
#'   \item \code{fmla_bin_ctrl}: nearest-neighbor pass 2 formula (control anchors)
#' }
#' In your current implementation, \code{fmla_bin_trt} uses only \code{vec[[1]]} and \code{fmla_bin_ctrl}
#' uses only \code{vec[[2]]}. This implies you are allowing different covariate sets for the two passes
#' (which can be intentional), but you should ensure \code{vec} is structured as a list with elements
#' \code{vec[[1]]} and \code{vec[[2]]} if you rely on this behavior.
#'
#' \strong{Full matching branch.}
#' Runs \code{MatchIt::matchit(..., method="full", estimand="ATE")} and then, within each matched
#' \code{subclass}, computes the mean outcome in treated and control units. Each unit’s opposite-arm
#' imputation is the mean outcome in the opposite arm within its subclass. The same is done for the
#' weight-like variable \code{compW}.
#'
#' \strong{Nearest-neighbor branch (two-pass).}
#' To obtain an opposite-arm imputation for every unit (treated and control) while keeping anchors
#' consistent:
#' \enumerate{
#'   \item \emph{ATT pass (treated anchors).} Match treated units to controls (\code{estimand="ATT"}).
#'   For each treated anchor, compute the mean control outcome (and mean \code{compW}) among its matched
#'   controls; assign these as \code{pairedCompY_ATT_for_treated} and \code{paired.ipcw.R_ATT_for_treated}.
#'   \item \emph{ATC pass (control anchors).} Re-run matching with control anchors by switching the
#'   “treated” definition in the formula to \code{I(tx==-1)}. For each control anchor, compute the mean
#'   treated outcome (and mean \code{compW}) among its matched treated units; assign these as
#'   \code{pairedCompY_ATC_for_controls} and \code{paired.ipcw.R_ATC_for_controls}.
#'   \item Merge the two summaries back to a one-row-per-ID dataset, and choose the appropriate opposite-arm
#'   imputation depending on the unit’s observed treatment.
#' }
#'
#' \strong{Constructed counterfactuals.}
#' After matching, the function constructs two counterfactual outcome columns for each unit:
#' \itemize{
#'   \item \code{trt_cf1L}: the unit’s outcome under treatment. Equals observed \code{compY} if treated,
#'   otherwise equals the matched opposite-arm imputation.
#'   \item \code{ctrl_cf1L}: the unit’s outcome under control. Equals observed \code{compY} if control,
#'   otherwise equals the matched opposite-arm imputation.
#' }
#'
#' \strong{Match-weight proxy.}
#' A simple distance-like proxy is computed as \code{diff.wt = compW - paired.ipcw.R} and
#' \code{match.weight = abs(diff.wt)}. This is not a \pkg{MatchIt} weight; it is a user-defined measure
#' intended to quantify discrepancy between the unit and its opposite-arm donor set in \code{compW}.
#'
#' \strong{Handling missing opposite-arm matches.}
#' If a unit fails to obtain a valid opposite-arm imputation (e.g., due to empty donor sets),
#' behavior is controlled by \code{na_handling}:
#' \itemize{
#'   \item \code{"drop"}: remove affected units from the output.
#'   \item \code{"zero_weight"}: keep them but set \code{match.weight = 0}.
#' }
#'
#' @param dat A data.frame containing at minimum the treatment indicator \code{txgroup}, outcome \code{compY},
#' identifier \code{Id}, and matching covariates in \code{vec} (or \code{vec[[1]]}, \code{vec[[2]]} for nearest).
#'
#' @param txgroup Character scalar. Name of the treatment indicator column. Internally recoded to \eqn{\{-1,+1\}}.
#'
#' @param exact_vars Optional. Exact matching specification passed to \pkg{MatchIt}. Can be:
#' \code{NULL}, a one-sided formula (e.g., \code{~ sex + site}), or a character vector of column names.
#'
#' @param compY Character scalar. Outcome column used to construct counterfactuals and paired outcomes.
#'
#' @param vec Matching covariates. For \code{method="full"}, \code{vec} should be a character vector of covariate names.
#' For \code{method="nearest"}, your current implementation expects \code{vec} to be a list with
#' \code{vec[[1]]} (covariates used when treated are anchors) and \code{vec[[2]]} (covariates used when controls are anchors).
#'
#' @param Id Character scalar. Subject identifier column name.
#'
#' @param method Matching method: \code{"nearest"} (two-pass ATT/ATC) or \code{"full"}.
#'
#' @param k Integer. Match ratio (number of donors per anchor) for nearest-neighbor matching. Ignored for full matching.
#'
#' @param replace Logical. Whether matching is performed with replacement (nearest-neighbor only).
#'
#' @param caliper Optional numeric. Caliper passed to \code{MatchIt::matchit} (applies to both matching types when supplied).
#'
#' @param distance Distance specification passed to \code{MatchIt::matchit} (e.g., \code{"mahalanobis"} or a numeric vector).
#'
#' @param compW Character scalar. Name of the “weight-like” column used to compute \code{diff.wt}. If missing,
#' it is created as a copy of \code{compY}. Default is \code{"ipcw.R"}.
#'
#' @param na_handling How to handle units with missing paired outcomes/weights: \code{"drop"} or \code{"zero_weight"}.
#'
#' @param plotbalance Logical. If \code{TRUE} and the \pkg{cobalt} package is installed, prints a Love
#'   plot of standardized mean differences before and after matching, providing a visual balance
#'   diagnostic. For \code{method="nearest"}, separate Love plots are produced for the ATT and ATC
#'   passes. Default is \code{FALSE}.
#'
#' @return A data.frame containing (at least) the original columns and additional derived columns:
#' \itemize{
#'   \item \code{pairedCompY}: imputed opposite-arm outcome for each unit (mean within matched set)
#'   \item \code{paired.ipcw.R}: imputed opposite-arm \code{compW} value (mean within matched set)
#'   \item \code{trt_cf1L}, \code{ctrl_cf1L}: constructed counterfactual outcomes under treatment/control
#'   \item \code{diff.wt}, \code{match.weight}: discrepancy in \code{compW} and its absolute value
#'   \item \code{newTxt}: sign-coded label based on treatment and \code{diff.wt}
#' }
#' When \code{plotbalance = TRUE} and \pkg{cobalt} is available, a \code{"balance"} attribute is
#' attached to the returned data.frame containing the \code{cobalt::bal.tab()} result(s) for
#' programmatic inspection.
#'
#' @seealso \code{\link[MatchIt]{matchit}}, \code{\link[MatchIt]{match.data}},
#'   \code{\link[MatchIt]{get_matches}}, \code{\link[cobalt]{love.plot}},
#'   \code{\link[cobalt]{bal.tab}}
#' @export
matchpotential_DTR <- function(dat, txgroup, exact_vars, compY, vec, Id,
                               method = c("nearest","full"),
                               k = 3, replace = TRUE, caliper = NULL,
                               distance = "mahalanobis",
                               compW = "ipcw.R",
                               na_handling = c("drop","zero_weight"),
                               plotbalance = FALSE){

  method <- match.arg(method)
  na_handling <- match.arg(na_handling)

  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' is required for matchpotential_DTR().", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for matchpotential_DTR().", call. = FALSE)
  }

  subHs <- dat

  # basic checks
  need <- c(txgroup, compY, Id)
  miss <- setdiff(need, names(subHs))
  if (length(miss) > 0L) stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

  # internal {-1,+1} coding
  subHs[[txgroup]] <- ifelse(subHs[[txgroup]] == 1, 1, -1)

  # compW fallback
  if (!compW %in% names(subHs)) subHs[[compW]] <- subHs[[compY]]

  # exact matching arg
  match_vars <- if (is.null(exact_vars)) NULL else {
    if (inherits(exact_vars, "formula")) exact_vars
    else if (is.character(exact_vars)) stats::reformulate(exact_vars)
    else stop("exact_vars must be NULL, a one-sided formula, or a character vector", call. = FALSE)
  }

  tx <- txgroup
  y  <- compY
  wname <- compW

  # Build formulas:
  # - full matching uses all covariates in vec (character vector)
  # - nearest matching uses vec[[1]] and vec[[2]] (list of two character vectors)
  if (method == "full") {
    if (!is.character(vec)) stop("For method='full', vec must be a character vector of covariate names.", call. = FALSE)
    if (!all(vec %in% names(subHs))) stop("Some 'vec' covariates not found in dat.", call. = FALSE)
    fmla_bin <- stats::as.formula(paste0("I(", tx, "==1) ~ ", paste(vec, collapse = " + ")))

  } else {
    if (!is.list(vec) || length(vec) < 2) stop("For method='nearest', vec must be a list with vec[[1]] and vec[[2]].", call. = FALSE)
    if (!all(vec[[1]] %in% names(subHs))) stop("Some vec[[1]] covariates not found in dat.", call. = FALSE)
    if (!all(vec[[2]] %in% names(subHs))) stop("Some vec[[2]] covariates not found in dat.", call. = FALSE)

    fmla_bin_trt  <- stats::as.formula(paste0("I(", tx, "==1) ~ ",  paste(vec[[1]], collapse = " + ")))
    fmla_bin_ctrl <- stats::as.formula(paste0("I(", tx, "==-1) ~ ", paste(vec[[2]], collapse = " + ")))
  }

  if (method == "full") {
    # --------------------
    # FULL MATCHING (ATE)
    # --------------------
    m.out <- MatchIt::matchit(
      formula  = fmla_bin,
      data     = subHs,
      method   = "full",
      distance = distance,
      caliper  = caliper,
      exact    = match_vars,
      estimand = "ATE"
    )

    if (plotbalance) {
      if (requireNamespace("cobalt", quietly = TRUE)) {
        bal <- cobalt::bal.tab(m.out, stats = c("mean.diffs", "variance.ratios"),
                               thresholds = c(m = 0.1, v = 2))
        print(cobalt::love.plot(m.out, stats = "mean.diffs", thresholds = c(m = 0.1),
                                abs = TRUE, stars = "std",
                                title = "Covariate Balance After Full Matching"))
      } else {
        message("Install the 'cobalt' package to enable balance diagnostics (plotbalance=TRUE).")
      }
    }

    dm <- MatchIt::match.data(m.out)

    dm <- dm %>%
      dplyr::group_by(.data$subclass) %>%
      dplyr::mutate(
        y_t_bar = mean(ifelse(.data[[tx]] ==  1, .data[[y]],     NA_real_), na.rm = TRUE),
        y_c_bar = mean(ifelse(.data[[tx]] == -1, .data[[y]],     NA_real_), na.rm = TRUE),
        w_t_bar = mean(ifelse(.data[[tx]] ==  1, .data[[wname]], NA_real_), na.rm = TRUE),
        w_c_bar = mean(ifelse(.data[[tx]] == -1, .data[[wname]], NA_real_), na.rm = TRUE),
        pairedCompY   = ifelse(.data[[tx]] == 1, y_c_bar, y_t_bar),
        paired.ipcw.R = ifelse(.data[[tx]] == 1, w_c_bar, w_t_bar)
      ) %>%
      dplyr::ungroup()

    out <- dm %>%
      dplyr::group_by(.data[[Id]]) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()

  } else {
    # -------------------------------------------
    # NEAREST MATCHING (two-pass: ATT then ATC)
    # -------------------------------------------

    # PASS 1: ATT (anchors = treated)
    m.att <- MatchIt::matchit(
      formula   = fmla_bin_trt,
      data      = subHs,
      method    = "nearest",
      distance  = distance,
      ratio     = k,
      replace   = replace,
      caliper   = caliper,
      exact     = match_vars,
      estimand  = "ATT"
    )

    md_att <- MatchIt::get_matches(m.att, data = subHs) %>%
      dplyr::group_by(.data$subclass) %>%
      dplyr::mutate(
        y_c_bar = mean(ifelse(.data[[tx]] == -1, .data[[y]],     NA_real_), na.rm = TRUE),
        w_c_bar = mean(ifelse(.data[[tx]] == -1, .data[[wname]], NA_real_), na.rm = TRUE),
        pairedCompY_ATT_for_treated   = ifelse(.data[[tx]] ==  1, y_c_bar, NA_real_),
        paired.ipcw.R_ATT_for_treated = ifelse(.data[[tx]] ==  1, w_c_bar, NA_real_)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(.data[[Id]]) %>%
      dplyr::summarise(
        pairedCompY_ATT_for_treated   = suppressWarnings(mean(pairedCompY_ATT_for_treated,   na.rm = TRUE)),
        paired.ipcw.R_ATT_for_treated = suppressWarnings(mean(paired.ipcw.R_ATT_for_treated, na.rm = TRUE)),
        .groups = "drop"
      )

    # PASS 2: ATC (anchors = controls)
    m.atc <- MatchIt::matchit(
      formula   = fmla_bin_ctrl,
      data      = subHs,
      method    = "nearest",
      distance  = distance,
      ratio     = k,
      replace   = replace,
      caliper   = caliper,
      exact     = match_vars,
      estimand  = "ATT"
    )

    if (plotbalance) {
      if (requireNamespace("cobalt", quietly = TRUE)) {
        bal_att <- cobalt::bal.tab(m.att, stats = c("mean.diffs", "variance.ratios"),
                                   thresholds = c(m = 0.1, v = 2))
        bal_atc <- cobalt::bal.tab(m.atc, stats = c("mean.diffs", "variance.ratios"),
                                   thresholds = c(m = 0.1, v = 2))
        print(cobalt::love.plot(m.att, stats = "mean.diffs", thresholds = c(m = 0.1),
                                abs = TRUE, stars = "std",
                                title = "Covariate Balance After Nearest Matching (ATT)"))
        print(cobalt::love.plot(m.atc, stats = "mean.diffs", thresholds = c(m = 0.1),
                                abs = TRUE, stars = "std",
                                title = "Covariate Balance After Nearest Matching (ATC)"))
      } else {
        message("Install the 'cobalt' package to enable balance diagnostics (plotbalance=TRUE).")
      }
    }

    md_atc <- MatchIt::get_matches(m.atc, data = subHs) %>%
      dplyr::group_by(.data$subclass) %>%
      dplyr::mutate(
        y_t_bar = mean(ifelse(.data[[tx]] ==  1, .data[[y]],     NA_real_), na.rm = TRUE),
        w_t_bar = mean(ifelse(.data[[tx]] ==  1, .data[[wname]], NA_real_), na.rm = TRUE),
        pairedCompY_ATC_for_controls   = ifelse(.data[[tx]] == -1, y_t_bar, NA_real_),
        paired.ipcw.R_ATC_for_controls = ifelse(.data[[tx]] == -1, w_t_bar, NA_real_)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(.data[[Id]]) %>%
      dplyr::summarise(
        pairedCompY_ATC_for_controls   = suppressWarnings(mean(pairedCompY_ATC_for_controls,   na.rm = TRUE)),
        paired.ipcw.R_ATC_for_controls = suppressWarnings(mean(paired.ipcw.R_ATC_for_controls, na.rm = TRUE)),
        .groups = "drop"
      )

    # base table: one row per unit (try to keep original columns)
    out <- subHs %>%
      dplyr::distinct() %>%
      dplyr::left_join(md_att, by = Id) %>%
      dplyr::left_join(md_atc, by = Id) %>%
      dplyr::mutate(
        pairedCompY = ifelse(.data[[tx]] ==  1,
                             pairedCompY_ATT_for_treated,
                             pairedCompY_ATC_for_controls),
        paired.ipcw.R = ifelse(.data[[tx]] ==  1,
                               paired.ipcw.R_ATT_for_treated,
                               paired.ipcw.R_ATC_for_controls)
      )
  }

  # Construct counterfactuals
  out$trt_cf1L  <- (as.integer(out[[tx]] ==  1) * out[[y]] +
                      as.integer(out[[tx]] != 1) * out$pairedCompY)
  out$ctrl_cf1L <- (as.integer(out[[tx]] == -1) * out[[y]] +
                      as.integer(out[[tx]] != -1) * out$pairedCompY)

  # Match-weight proxy
  out$diff.wt      <- out[[wname]] - out$paired.ipcw.R
  out$match.weight <- abs(out$diff.wt)

  # Handle residual NAs
  bad <- is.na(out$paired.ipcw.R) | is.na(out$pairedCompY) | is.na(out[[wname]])
  if (any(bad)) {
    nbad <- sum(bad)
    if (na_handling == "drop") {
      out <- out[!bad, , drop = FALSE]
      message(sprintf("Dropped %d ID(s) with no valid opposite-arm match.", nbad))
    } else {
      out$match.weight[bad] <- 0
      message(sprintf("Kept %d ID(s) with no valid opposite-arm match; set match.weight = 0.", nbad))
    }
  }

  # Labels: robust to NA/zero
  sgn <- sign(out$diff.wt); sgn[is.na(sgn)] <- 0
  out$newTxt <- out[[tx]] * sgn
  idx0 <- (!is.na(out$diff.wt)) & (out$diff.wt == 0)
  if (any(idx0)) out$newTxt[idx0] <- out[[tx]][idx0]

  # Ensure numeric
  out$pairedCompY   <- as.numeric(out$pairedCompY)
  out$paired.ipcw.R <- as.numeric(out$paired.ipcw.R)

  if (plotbalance && requireNamespace("cobalt", quietly = TRUE)) {
    if (method == "full") {
      attr(out, "balance") <- bal
    } else {
      attr(out, "balance") <- list(ATT = bal_att, ATC = bal_atc)
    }
  }

  out
}
