# --------------------------------------------------------------------------
# Measurement Invariance Testing
# Configural -> Metric -> Scalar -> Strict
# --------------------------------------------------------------------------

#' Test Measurement Invariance Across Groups
#'
#' Fits a sequence of increasingly constrained multi-group CFA models and
#' compares them via likelihood ratio tests (or scaled difference tests for
#' robust estimators). This tests whether the factor structure is stable
#' across biological groups (e.g., patients vs. controls).
#'
#' @section Invariance Levels:
#' \describe{
#'   \item{Configural}{Same factor structure in both groups (loadings free)}
#'   \item{Metric (weak)}{Equal factor loadings across groups — required for
#'     meaningful group comparisons on the latent factor}
#'   \item{Scalar (strong)}{Equal loadings AND intercepts — required for
#'     comparing latent means across groups}
#'   \item{Strict}{Equal loadings, intercepts, AND residual variances}
#' }
#'
#' @param fit A \code{biomimic_fit} object, or a \code{biomimic_selection}
#'   object (in which case \code{group_var} is required).
#' @param group_var Character, name of the grouping variable in the data.
#'   Required when \code{fit} is a selection object. For \code{biomimic_fit}
#'   objects of type \code{"multigroup"}, extracted automatically.
#' @param estimator Character, lavaan estimator (default \code{"MLR"}).
#' @param K Integer, number of factors. If \code{NULL}, uses the fitted K.
#'
#' @return A list of class \code{biomimic_invariance} with components:
#' \describe{
#'   \item{fits}{Named list of lavaan fit objects (configural, metric, scalar, strict)}
#'   \item{comparisons}{Data frame with LRT comparisons between successive models}
#'   \item{fit_indices}{Data frame with GOF indices for each level}
#'   \item{verdict}{Character vector with plain-language assessment per level}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fit_sem(sel, type = "multigroup", group_var = "group")
#' inv <- test_invariance(fit)   # group_var extracted automatically
#' print(inv)
#' }
#'
#' @export
test_invariance <- function(fit, group_var = NULL,
                             estimator = "MLR", K = NULL) {

    # Accept both biomimic_fit and biomimic_selection
    if (inherits(fit, "biomimic_fit")) {
        selection <- fit$selection
        # Extract group_var from lavaan fit if multigroup
        if (is.null(group_var)) {
            lav_group <- fit$lavaan_fit@Data@group
            if (length(lav_group) > 0 && nchar(lav_group) > 0) {
                group_var <- lav_group
            }
        }
        if (is.null(estimator) || estimator == "MLR") {
            estimator <- fit$lavaan_fit@Options$estimator
        }
    } else if (inherits(fit, "biomimic_selection")) {
        selection <- fit
    } else {
        stop("fit must be a biomimic_fit or biomimic_selection object.")
    }

    if (is.null(group_var)) {
        stop("group_var is required for invariance testing.")
    }
    stopifnot(is.character(group_var), length(group_var) == 1)
    stopifnot(group_var %in% names(selection$data))

    if (is.null(K)) K <- max(selection$K, 1L)
    vars <- selection$variables

    .log_info("test_invariance: K=%d, group_var='%s', %d indicators",
              K, group_var, length(vars))

    # Build CFA syntax
    if (K == 1L) {
        model_syntax <- paste0("F1 =~ ", paste(vars, collapse = " + "))
    } else {
        groups <- selection$indicator_groups
        lines <- character()
        for (k in seq_len(K)) {
            factor_vars <- vars[groups[[k]]]
            lines <- c(lines,
                       paste0("F", k, " =~ ", paste(factor_vars, collapse = " + ")))
        }
        model_syntax <- paste(lines, collapse = "\n")
    }

    # Fit configural (no constraints across groups)
    fit_configural <- lavaan::cfa(model_syntax, data = selection$data,
                                   group = group_var, estimator = estimator)

    # Fit metric (equal loadings)
    fit_metric <- lavaan::cfa(model_syntax, data = selection$data,
                               group = group_var, estimator = estimator,
                               group.equal = "loadings")

    # Fit scalar (equal loadings + intercepts)
    fit_scalar <- lavaan::cfa(model_syntax, data = selection$data,
                               group = group_var, estimator = estimator,
                               group.equal = c("loadings", "intercepts"))

    # Fit strict (equal loadings + intercepts + residual VARIANCES).
    # The lavaan keyword for equal residual variances is "residuals";
    # "residual.covariances" would constrain only the (zero, fixed) off-diagonal
    # residual covariances and leave strict identical to scalar.
    fit_strict <- lavaan::cfa(model_syntax, data = selection$data,
                               group = group_var, estimator = estimator,
                               group.equal = c("loadings", "intercepts",
                                               "residuals"))

    fits <- list(
        configural = fit_configural,
        metric     = fit_metric,
        scalar     = fit_scalar,
        strict     = fit_strict
    )

    # --- LRT comparisons ---
    comparisons <- data.frame(
        comparison = character(0),
        chisq_diff = numeric(0),
        df_diff    = numeric(0),
        pvalue     = numeric(0),
        stringsAsFactors = FALSE
    )

    pairs <- list(
        c("configural", "metric"),
        c("metric", "scalar"),
        c("scalar", "strict")
    )

    for (pair in pairs) {
        lrt <- tryCatch(
            lavaan::lavTestLRT(fits[[pair[1]]], fits[[pair[2]]]),
            error = function(e) NULL
        )
        if (!is.null(lrt) && nrow(lrt) >= 2) {
            comparisons <- rbind(comparisons, data.frame(
                comparison = paste(pair[1], "->", pair[2]),
                chisq_diff = lrt[["Chisq diff"]][2],
                df_diff    = lrt[["Df diff"]][2],
                pvalue     = lrt[["Pr(>Chisq)"]][2],
                stringsAsFactors = FALSE
            ))
        }
    }

    # --- Fit indices table ---
    indices_names <- c("chisq", "df", "rmsea", "cfi", "tli", "srmr")
    fit_indices <- data.frame(level = names(fits), stringsAsFactors = FALSE)
    for (idx in indices_names) {
        fit_indices[[idx]] <- vapply(fits, function(f) {
            fm <- lavaan::fitMeasures(f, idx)
            unname(fm)
        }, numeric(1))
    }

    # --- Delta CFI criterion (Cheung & Rensvold, 2002) ---
    delta_cfi <- diff(fit_indices$cfi)
    verdicts <- character(3)
    comp_names <- c("Metric", "Scalar", "Strict")
    for (i in seq_len(3)) {
        p_val <- if (i <= nrow(comparisons)) comparisons$pvalue[i] else NA
        d_cfi <- delta_cfi[i]
        if (is.na(p_val)) {
            verdicts[i] <- "Could not compute"
        } else if (p_val > 0.05 || abs(d_cfi) < 0.01) {
            verdicts[i] <- "SUPPORTED"
        } else {
            verdicts[i] <- "NOT SUPPORTED"
        }
    }
    names(verdicts) <- comp_names

    for (nm in comp_names) {
        .log_info("test_invariance: %s -> %s", nm, verdicts[nm])
    }

    structure(
        list(
            fits         = fits,
            comparisons  = comparisons,
            fit_indices  = fit_indices,
            verdict      = verdicts
        ),
        class = "biomimic_invariance"
    )
}


#' Print Method for biomimic_invariance
#'
#' @param x A \code{biomimic_invariance} object.
#' @param ... Ignored.
#' @export
print.biomimic_invariance <- function(x, ...) {
    cat("=== Measurement Invariance Test ===\n")
    cat("Tests whether the protein-factor structure is stable across groups.\n\n")

    # Fit indices table
    cat("--- Fit Indices by Invariance Level ---\n")
    fi <- x$fit_indices
    cat(sprintf("  %-12s  %7s  %3s  %6s  %6s  %6s  %6s\n",
                "Level", "Chi-sq", "df", "RMSEA", "CFI", "TLI", "SRMR"))
    for (i in seq_len(nrow(fi))) {
        cat(sprintf("  %-12s  %7.1f  %3d  %6.3f  %6.3f  %6.3f  %6.3f\n",
                    fi$level[i], fi$chisq[i], fi$df[i], fi$rmsea[i],
                    fi$cfi[i], fi$tli[i], fi$srmr[i]))
    }
    cat("\n")

    # LRT comparisons
    cat("--- Likelihood Ratio Tests ---\n")
    for (i in seq_len(nrow(x$comparisons))) {
        row <- x$comparisons[i, ]
        cat(sprintf("  %s: delta-Chi2 = %.2f, delta-df = %d, p = %s\n",
                    row$comparison, row$chisq_diff, row$df_diff,
                    .fmt_pval(row$pvalue)))
    }
    cat("\n")

    # Verdicts
    cat("--- Assessment ---\n")
    meanings <- c(
        Metric = "Equal loadings? (needed for group comparisons)",
        Scalar = "Equal intercepts? (needed for latent mean comparisons)",
        Strict = "Equal residuals? (strongest form)"
    )
    for (nm in names(x$verdict)) {
        icon <- if (x$verdict[nm] == "SUPPORTED") "[PASS]" else "[FAIL]"
        cat(sprintf("  %s %s: %s\n    %s\n",
                    icon, nm, x$verdict[nm], meanings[nm]))
    }
    cat("\nCriteria: LRT p > 0.05 or |delta-CFI| < 0.01 => supported.\n")

    invisible(x)
}
