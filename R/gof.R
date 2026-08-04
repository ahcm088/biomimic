# --------------------------------------------------------------------------
# Goodness-of-Fit Diagnostics
# RMSEA, CFI, TLI, SRMR, residual analysis, Mahalanobis distance
# --------------------------------------------------------------------------

#' Goodness-of-Fit Indices
#'
#' Extracts fit indices from a \code{biomimic_fit} object with plain-language
#' interpretation guidelines.
#'
#' @param fit A \code{biomimic_fit} object.
#'
#' @return A list of class \code{biomimic_gof} with fit indices and
#'   interpretation.
#'
#' @export
goodness_of_fit <- function(fit) {
    stopifnot(inherits(fit, "biomimic_fit"))

    fm <- lavaan::fitMeasures(fit$lavaan_fit,
                               c("chisq", "df", "pvalue",
                                 "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                                 "rmsea.pvalue",
                                 "cfi", "tli",
                                 "srmr",
                                 "aic", "bic"))

    structure(
        list(
            chisq          = unname(fm["chisq"]),
            df             = unname(fm["df"]),
            pvalue         = unname(fm["pvalue"]),
            rmsea          = unname(fm["rmsea"]),
            rmsea.ci.lower = unname(fm["rmsea.ci.lower"]),
            rmsea.ci.upper = unname(fm["rmsea.ci.upper"]),
            rmsea.pvalue   = unname(fm["rmsea.pvalue"]),
            cfi            = unname(fm["cfi"]),
            tli            = unname(fm["tli"]),
            srmr           = unname(fm["srmr"]),
            aic            = unname(fm["aic"]),
            bic            = unname(fm["bic"])
        ),
        class = "biomimic_gof"
    )
}


#' Print Method for biomimic_gof
#' @param x A \code{biomimic_gof} object.
#' @param ... Ignored.
#' @export
print.biomimic_gof <- function(x, ...) {
    cat("=== Model Fit Assessment ===\n\n")

    cat(sprintf("Chi-squared = %.2f (df = %d, p = %s)\n",
                x$chisq, x$df, .fmt_pval(x$pvalue)))
    cat(sprintf("  Note: sensitive to sample size; prefer RMSEA/CFI/SRMR.\n\n"))

    cat(sprintf("RMSEA = %.3f  90%% CI [%.3f, %.3f]\n",
                x$rmsea, x$rmsea.ci.lower, x$rmsea.ci.upper))
    cat(sprintf("  %s\n", .interpret_rmsea(x$rmsea)))
    cat(sprintf("  Close-fit test (H0: RMSEA <= 0.05): p = %s\n\n",
                .fmt_pval(x$rmsea.pvalue)))

    cat(sprintf("CFI  = %.3f  %s\n", x$cfi, .interpret_cfi(x$cfi)))
    cat(sprintf("TLI  = %.3f  %s\n", x$tli, .interpret_tli(x$tli)))
    cat(sprintf("SRMR = %.3f  %s\n\n", x$srmr, .interpret_srmr(x$srmr)))

    cat(sprintf("AIC = %.1f | BIC = %.1f\n", x$aic, x$bic))
    cat("  (lower is better; use for model comparison)\n")

    invisible(x)
}


#' Residual Analysis for biomimic_fit
#'
#' Computes standardized residuals from the fitted covariance matrix.
#' Large residuals indicate local misfit (pairs of variables whose
#' covariance is poorly explained by the model).
#'
#' @param fit A \code{biomimic_fit} object.
#' @param type Character, \code{"cor.bollen"} (default) for standardized
#'   residual correlations or \code{"raw"} for raw residual covariances.
#'
#' @return A list with:
#' \describe{
#'   \item{residuals}{Matrix of residual (co)variances or correlations}
#'   \item{summary}{Data frame of the largest absolute residuals}
#' }
#'
#' @export
residual_analysis <- function(fit, type = c("cor.bollen", "raw")) {
    type <- match.arg(type)
    stopifnot(inherits(fit, "biomimic_fit"))

    resid_mat <- lavaan::lavResiduals(fit$lavaan_fit, type = type)
    cov_resid <- resid_mat$cov

    # Extract upper triangle as data frame of pairs
    p <- nrow(cov_resid)
    vnames <- rownames(cov_resid)
    pairs <- data.frame(
        var1     = character(0),
        var2     = character(0),
        residual = numeric(0),
        stringsAsFactors = FALSE
    )

    for (i in seq_len(p - 1)) {
        for (j in (i + 1):p) {
            pairs <- rbind(pairs, data.frame(
                var1     = vnames[i],
                var2     = vnames[j],
                residual = cov_resid[i, j],
                stringsAsFactors = FALSE
            ))
        }
    }

    pairs$abs_residual <- abs(pairs$residual)
    pairs <- pairs[order(-pairs$abs_residual), ]

    # Add biology labels
    pairs$var1_label <- .bio_label(pairs$var1, fit$indicator_labels)
    pairs$var2_label <- .bio_label(pairs$var2, fit$indicator_labels)

    rownames(pairs) <- NULL

    list(
        residuals = cov_resid,
        summary   = pairs
    )
}


#' Mahalanobis Distance for Outlier Detection
#'
#' Computes the Mahalanobis distance for each observation based on the
#' model-implied covariance matrix. Observations with large distances
#' are potential outliers that may unduly influence the factor model.
#'
#' @param fit A \code{biomimic_fit} object.
#' @param p_threshold Numeric, p-value threshold for flagging outliers
#'   (default 0.001 based on chi-squared distribution).
#'
#' @return A data frame with columns: \code{observation}, \code{mahal_dist},
#'   \code{pvalue}, \code{outlier} (logical flag).
#'
#' @export
mahalanobis_distance <- function(fit, p_threshold = 0.001) {
    stopifnot(inherits(fit, "biomimic_fit"))

    # Get observed data (indicator variables only)
    vars <- fit$selection$variables
    dat <- fit$selection$data[, vars, drop = FALSE]
    dat <- as.matrix(dat)

    # Model-implied covariance (ensure matching variable order)
    sigma_hat <- lavaan::lavInspect(fit$lavaan_fit, "cov.ov")
    sigma_hat <- sigma_hat[vars, vars]
    mu_hat <- lavaan::lavInspect(fit$lavaan_fit, "mean.ov")

    # If means not estimated, use column means
    if (length(mu_hat) == 0) {
        mu_hat <- colMeans(dat)
    } else {
        mu_hat <- mu_hat[vars]
    }

    # Mahalanobis distance
    n <- nrow(dat)
    p <- length(vars)
    sigma_inv <- solve(sigma_hat)

    # Centre the data
    dat_c <- sweep(dat, 2, mu_hat)
    d2 <- rowSums((dat_c %*% sigma_inv) * dat_c)

    # P-value from chi-squared distribution
    pval <- stats::pchisq(d2, df = p, lower.tail = FALSE)

    data.frame(
        observation = seq_len(n),
        mahal_dist  = d2,
        pvalue      = pval,
        outlier     = pval < p_threshold,
        stringsAsFactors = FALSE
    )
}


# ===========================================================================
# Interpretation helpers
# ===========================================================================

#' @noRd
.interpret_rmsea <- function(rmsea) {
    if (is.na(rmsea)) return("")
    if (rmsea <= 0.05) return("[Excellent fit]")
    if (rmsea <= 0.08) return("[Acceptable fit]")
    if (rmsea <= 0.10) return("[Mediocre fit]")
    "[Poor fit]"
}

#' @noRd
.interpret_cfi <- function(cfi) {
    if (is.na(cfi)) return("")
    if (cfi >= 0.95) return("[Excellent]")
    if (cfi >= 0.90) return("[Acceptable]")
    "[Below threshold]"
}

#' @noRd
.interpret_tli <- function(tli) {
    if (is.na(tli)) return("")
    if (tli >= 0.95) return("[Excellent]")
    if (tli >= 0.90) return("[Acceptable]")
    "[Below threshold]"
}

#' @noRd
.interpret_srmr <- function(srmr) {
    if (is.na(srmr)) return("")
    if (srmr <= 0.05) return("[Excellent fit]")
    if (srmr <= 0.08) return("[Acceptable fit]")
    "[Poor fit]"
}
