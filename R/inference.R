# --------------------------------------------------------------------------
# Inference Tools for biomimic_fit
# Biology-friendly summaries, confidence intervals, hypothesis tests
# --------------------------------------------------------------------------

#' Summary Method for biomimic_fit
#'
#' Produces a biology-friendly summary of the fitted SEM model, translating
#' SEM parameters into domain-relevant language for biologists and
#' bioinformaticians.
#'
#' @param object A \code{biomimic_fit} object.
#' @param standardized Logical, whether to include standardized estimates
#'   (default \code{TRUE}).
#' @param ci Logical, whether to include confidence intervals
#'   (default \code{TRUE}).
#' @param ci.level Numeric, confidence level (default 0.95).
#' @param ... Ignored.
#'
#' @return Invisibly returns a list of class \code{biomimic_summary} with
#'   components: \code{loadings}, \code{structural}, \code{direct_effects},
#'   \code{gof}, \code{factor_labels}, \code{indicator_labels}.
#'
#' @export
summary.biomimic_fit <- function(object, standardized = TRUE, ci = TRUE,
                                  ci.level = 0.95, ...) {

    fit <- object$lavaan_fit
    pe <- lavaan::parameterEstimates(fit, standardized = standardized,
                                      ci = ci, level = ci.level)

    # --- Loadings (measurement model) ---
    loadings <- pe[pe$op == "=~", , drop = FALSE]
    loadings$indicator <- .bio_label(loadings$rhs, object$indicator_labels)
    loadings$factor <- .bio_label(loadings$lhs, .factor_label_map(object))

    # --- Structural coefficients (regressions on latent variables) ---
    structural <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]
    structural$factor <- .bio_label(structural$lhs, .factor_label_map(object))
    structural$predictor <- structural$rhs

    # --- Direct effects (regressions on observed variables) ---
    direct_effects <- pe[pe$op == "~" & !grepl("^F[0-9]+$", pe$lhs), ,
                         drop = FALSE]
    if (nrow(direct_effects) > 0) {
        direct_effects$indicator <- .bio_label(direct_effects$lhs,
                                                object$indicator_labels)
    }

    # --- Residual variances ---
    residuals <- pe[pe$op == "~~" & pe$lhs == pe$rhs &
                    !grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]
    residuals$indicator <- .bio_label(residuals$lhs, object$indicator_labels)

    # --- Factor (co)variances ---
    factor_var <- pe[pe$op == "~~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]

    # --- Goodness of fit ---
    gof <- goodness_of_fit(object)

    result <- structure(
        list(
            loadings         = loadings,
            structural       = structural,
            direct_effects   = direct_effects,
            residuals        = residuals,
            factor_var       = factor_var,
            gof              = gof,
            model_type       = object$model_type,
            factor_labels    = object$factor_labels,
            indicator_labels = object$indicator_labels,
            K                = object$K,
            n                = lavaan::nobs(fit),
            estimator        = fit@Options$estimator
        ),
        class = "biomimic_summary"
    )

    print(result)
    invisible(result)
}


#' Print Method for biomimic_summary
#'
#' @param x A \code{biomimic_summary} object.
#' @param digits Number of decimal places (default 3).
#' @param ... Ignored.
#' @export
print.biomimic_summary <- function(x, digits = 3, ...) {

    # --- Header ---
    cat("=== biomimic SEM Results ===\n")
    cat(sprintf("Model: %s | Estimator: %s | n = %d\n\n",
                toupper(x$model_type), x$estimator, x$n))

    # --- Goodness of Fit ---
    cat("--- Model Fit ---\n")
    gof <- x$gof
    cat(sprintf("  Chi-sq = %.1f (df = %d, p = %s)\n",
                gof$chisq, gof$df, .fmt_pval(gof$pvalue)))
    cat(sprintf("  RMSEA  = %.3f [%.3f, %.3f]  %s\n",
                gof$rmsea, gof$rmsea.ci.lower, gof$rmsea.ci.upper,
                .interpret_rmsea(gof$rmsea)))
    cat(sprintf("  CFI    = %.3f  %s\n", gof$cfi, .interpret_cfi(gof$cfi)))
    cat(sprintf("  SRMR   = %.3f  %s\n", gof$srmr, .interpret_srmr(gof$srmr)))
    cat("\n")

    # --- Structural (group effects on latent pathways) ---
    if (nrow(x$structural) > 0) {
        cat("--- Group Effects on Latent Pathways ---\n")
        for (i in seq_len(nrow(x$structural))) {
            row <- x$structural[i, ]
            sig <- .signif_stars(row$pvalue)
            ci_str <- if ("ci.lower" %in% names(row)) {
                sprintf("[%.3f, %.3f]", row$ci.lower, row$ci.upper)
            } else ""
            cat(sprintf("  %s -> %s: est = %.3f, SE = %.3f, z = %.2f, p = %s %s  %s\n",
                        row$predictor, row$factor,
                        row$est, row$se, row$z, .fmt_pval(row$pvalue),
                        sig, ci_str))
        }
        cat("\n")
    }

    # --- Loadings (protein-pathway associations) ---
    cat("--- Protein-Pathway Associations (Factor Loadings) ---\n")
    for (i in seq_len(nrow(x$loadings))) {
        row <- x$loadings[i, ]
        sig <- .signif_stars(row$pvalue)
        std_str <- if ("std.all" %in% names(row)) {
            sprintf(" (std = %.3f)", row$std.all)
        } else ""
        cat(sprintf("  %s <- %s: %.3f%s %s\n",
                    row$indicator, row$factor, row$est, std_str, sig))
    }
    cat("\n")

    # --- Direct effects ---
    if (nrow(x$direct_effects) > 0) {
        cat("--- Direct Effects (not mediated by latent pathway) ---\n")
        for (i in seq_len(nrow(x$direct_effects))) {
            row <- x$direct_effects[i, ]
            sig <- .signif_stars(row$pvalue)
            cat(sprintf("  %s -> %s: est = %.3f, p = %s %s\n",
                        row$rhs, row$indicator, row$est,
                        .fmt_pval(row$pvalue), sig))
        }
        cat("\n")
    }

    cat("Signif. codes: '***' < 0.001 '**' < 0.01 '*' < 0.05\n")
    cat("Access the full lavaan object via $lavaan_fit for advanced output.\n")

    invisible(x)
}


#' Extract Coefficients from biomimic_fit
#'
#' Extracts parameter estimates in a tidy data frame with biology-friendly
#' labels.
#'
#' @param object A \code{biomimic_fit} object.
#' @param type Character, which parameters to extract:
#'   \code{"loadings"}, \code{"structural"}, \code{"direct"}, \code{"all"}
#'   (default \code{"all"}).
#' @param standardized Logical, include standardized estimates (default TRUE).
#' @param ... Ignored.
#'
#' @return A data frame with columns: \code{parameter}, \code{from}, \code{to},
#'   \code{estimate}, \code{se}, \code{z}, \code{pvalue}, \code{ci_lower},
#'   \code{ci_upper}, and optionally \code{std_estimate}.
#'
#' @export
coef.biomimic_fit <- function(object, type = c("all", "loadings",
                                                "structural", "direct"),
                               standardized = TRUE, ...) {
    type <- match.arg(type)
    fit <- object$lavaan_fit

    pe <- lavaan::parameterEstimates(fit, standardized = standardized,
                                      ci = TRUE, level = 0.95)

    # Filter by type
    if (type == "loadings") {
        pe <- pe[pe$op == "=~", , drop = FALSE]
    } else if (type == "structural") {
        pe <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]
    } else if (type == "direct") {
        pe <- pe[pe$op == "~" & !grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]
    } else {
        pe <- pe[pe$op %in% c("=~", "~"), , drop = FALSE]
    }

    if (nrow(pe) == 0) return(pe)

    # Build tidy output
    out <- data.frame(
        parameter = paste(pe$lhs, pe$op, pe$rhs),
        type      = ifelse(pe$op == "=~", "loading",
                    ifelse(grepl("^F[0-9]+$", pe$lhs), "structural", "direct")),
        from      = pe$rhs,
        to        = pe$lhs,
        estimate  = pe$est,
        se        = pe$se,
        z         = pe$z,
        pvalue    = pe$pvalue,
        ci_lower  = pe$ci.lower,
        ci_upper  = pe$ci.upper,
        stringsAsFactors = FALSE
    )

    if (standardized && "std.all" %in% names(pe)) {
        out$std_estimate <- pe$std.all
    }

    # Add biology labels
    out$from_label <- .bio_label(out$from, c(object$indicator_labels,
                                              .factor_label_map(object)))
    out$to_label   <- .bio_label(out$to, c(object$indicator_labels,
                                             .factor_label_map(object)))

    rownames(out) <- NULL
    out
}


#' Confidence Intervals for biomimic_fit
#'
#' @param object A \code{biomimic_fit} object.
#' @param parm Not used (for S3 compatibility).
#' @param level Confidence level (default 0.95).
#' @param ... Ignored.
#'
#' @return A data frame with confidence intervals for all model parameters.
#' @export
confint.biomimic_fit <- function(object, parm, level = 0.95, ...) {
    pe <- lavaan::parameterEstimates(object$lavaan_fit, ci = TRUE,
                                      level = level)
    pe <- pe[pe$op %in% c("=~", "~"), , drop = FALSE]

    out <- data.frame(
        parameter = paste(pe$lhs, pe$op, pe$rhs),
        estimate  = pe$est,
        ci_lower  = pe$ci.lower,
        ci_upper  = pe$ci.upper,
        level     = level,
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}


#' Likelihood Ratio Test Between Nested biomimic Models
#'
#' Compares two nested \code{biomimic_fit} objects using the likelihood ratio
#' test (or its robust equivalent when MLR is used).
#'
#' @param fit_full A \code{biomimic_fit} object (the less restricted model).
#' @param fit_restricted A \code{biomimic_fit} object (the more restricted model).
#'
#' @return A data frame with the test statistic, df, and p-value.
#'
#' @export
lr_test <- function(fit_full, fit_restricted) {
    stopifnot(inherits(fit_full, "biomimic_fit"),
              inherits(fit_restricted, "biomimic_fit"))

    anova_result <- lavaan::lavTestLRT(fit_restricted$lavaan_fit,
                                        fit_full$lavaan_fit)

    data.frame(
        chisq_diff = anova_result[["Chisq diff"]][2],
        df_diff    = anova_result[["Df diff"]][2],
        pvalue     = anova_result[["Pr(>Chisq)"]][2],
        stringsAsFactors = FALSE
    )
}


#' Print Method for biomimic_fit
#'
#' @param x A \code{biomimic_fit} object.
#' @param ... Ignored.
#' @export
print.biomimic_fit <- function(x, ...) {
    cat("biomimic SEM fit\n")
    cat(sprintf("  Model: %s (K = %d)\n", x$model_type, x$K))
    cat(sprintf("  Indicators: %d\n", length(x$selection$variables)))
    cat(sprintf("  Estimator: %s\n", x$lavaan_fit@Options$estimator))
    converged <- lavaan::lavInspect(x$lavaan_fit, "converged")
    cat(sprintf("  Converged: %s\n", converged))
    if (converged) {
        gof <- goodness_of_fit(x)
        cat(sprintf("  RMSEA = %.3f, CFI = %.3f, SRMR = %.3f\n",
                    gof$rmsea, gof$cfi, gof$srmr))
    }
    cat("\nUse summary() for biology-friendly results.\n")
    cat("Access lavaan object via $lavaan_fit.\n")
    invisible(x)
}


# ===========================================================================
# Internal helpers
# ===========================================================================

#' Map variable names to biology-friendly labels
#' @noRd
.bio_label <- function(names, label_map) {
    if (is.null(label_map)) return(names)
    out <- label_map[names]
    # Fall back to original name if no mapping
    missing <- is.na(out)
    out[missing] <- names[missing]
    unname(out)
}

#' Build factor label map from biomimic_fit
#' @noRd
.factor_label_map <- function(fit) {
    K <- fit$K
    fnames <- paste0("F", seq_len(K))
    setNames(fit$factor_labels, fnames)
}

#' Format p-value for display
#' @noRd
.fmt_pval <- function(p) {
    ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "< 0.001",
    ifelse(p < 0.01, sprintf("%.3f", p),
    sprintf("%.3f", p))))
}

#' Significance stars
#' @noRd
.signif_stars <- function(p) {
    ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
    ifelse(p < 0.05, "*", ""))))
}
