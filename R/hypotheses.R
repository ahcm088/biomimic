# --------------------------------------------------------------------------
# Unified Hypothesis Testing and Diagnostics
# Dispatches per model type so the user always calls the same function
# --------------------------------------------------------------------------

#' Test Hypotheses for a Fitted SEM Model
#'
#' Runs the pre-defined set of hypothesis tests appropriate for the model
#' type. The user does not need to know which tests apply — the function
#' selects them automatically based on the model specification.
#'
#' @section Tests by model type:
#' \describe{
#'   \item{MIMIC (K=1)}{
#'     Wald tests for each covariate effect on the latent factor (H0: B_k = 0).
#'     Omnibus test for any group effect.}
#'   \item{Multi-factor MIMIC (K>=2)}{
#'     Per-factor Wald tests (H0: B_group,k = 0 for each k).
#'     Omnibus test across all factors.
#'     Factor correlation test (H0: factors are uncorrelated).}
#'   \item{Multi-group CFA}{
#'     Measurement invariance hierarchy: configural -> metric -> scalar -> strict.
#'     Each level tested via scaled LRT (delta-chi2) and delta-CFI.}
#'   \item{MIMIC + direct effects}{
#'     Structural effect tests (H0: B_group = 0 via latent pathway).
#'     Direct effect tests (H0: gamma_j = 0 for each flagged protein).
#'     Comparison of indirect vs. direct effect magnitudes.}
#' }
#'
#' @param fit A \code{biomimic_fit} object.
#' @param level Significance level for hypothesis tests (default 0.05).
#' @param ... Additional arguments (currently unused).
#'
#' @return A list of class \code{biomimic_hypotheses} with components:
#' \describe{
#'   \item{model_type}{Character, which model was tested}
#'   \item{tests}{Data frame of individual test results}
#'   \item{omnibus}{Omnibus test result (if applicable)}
#'   \item{invariance}{Invariance results (multigroup only)}
#'   \item{level}{Significance level used}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fit_sem(sel, type = "mimic")
#' test_hypotheses(fit)
#'
#' fit_mg <- fit_sem(sel, type = "multigroup", group_var = "group")
#' test_hypotheses(fit_mg)  # automatically runs invariance hierarchy
#' }
#'
#' @export
test_hypotheses <- function(fit, level = 0.05, ...) {
    stopifnot(inherits(fit, "biomimic_fit"))

    .log_info("test_hypotheses: model_type=%s, level=%.2f", fit$model_type, level)

    result <- switch(fit$model_type,
        mimic       = .test_mimic(fit, level),
        multi_mimic = .test_multi_mimic(fit, level),
        multigroup  = .test_multigroup(fit, level),
        direct      = .test_direct(fit, level),
        stop("Unknown model type: ", fit$model_type)
    )

    result$model_type <- fit$model_type
    result$level <- level
    class(result) <- "biomimic_hypotheses"

    # Log summary of results
    if (!is.null(result$tests) && nrow(result$tests) > 0) {
        n_sig <- sum(result$tests$significant)
        .log_info("test_hypotheses: %d/%d tests significant",
                  n_sig, nrow(result$tests))
    }

    result
}


#' Run Model Diagnostics
#'
#' Runs a comprehensive diagnostic suite appropriate for the model type.
#' Includes goodness-of-fit, residual analysis, outlier detection, and
#' model-specific checks.
#'
#' @section Diagnostics by model type:
#' \describe{
#'   \item{All models}{Goodness-of-fit (RMSEA, CFI, SRMR), standardized
#'     residual correlations, Mahalanobis outlier detection.}
#'   \item{MIMIC}{Additionally: reliability (omega), explained variance per
#'     indicator (R-squared).}
#'   \item{Multi-group CFA}{Additionally: per-group fit indices, factor
#'     congruence (Tucker's phi).}
#'   \item{MIMIC + direct effects}{Additionally: comparison of models with
#'     and without direct effects (LRT).}
#' }
#'
#' @param fit A \code{biomimic_fit} object.
#' @param p_threshold P-value threshold for outlier detection (default 0.001).
#' @param ... Additional arguments (currently unused).
#'
#' @return A list of class \code{biomimic_diagnostics} with components that
#'   vary by model type.
#'
#' @export
diagnose <- function(fit, p_threshold = 0.001, ...) {
    stopifnot(inherits(fit, "biomimic_fit"))

    .log_info("diagnose: running diagnostics for %s model", fit$model_type)

    # Common diagnostics (all models)
    gof <- goodness_of_fit(fit)
    resid <- residual_analysis(fit)
    mahal <- mahalanobis_distance(fit, p_threshold = p_threshold)

    # R-squared per indicator
    r2 <- tryCatch(
        lavaan::lavInspect(fit$lavaan_fit, "rsquare"),
        error = function(e) NULL
    )

    # Model-specific diagnostics
    specific <- switch(fit$model_type,
        mimic       = .diagnose_mimic(fit),
        multi_mimic = .diagnose_multi_mimic(fit),
        multigroup  = .diagnose_multigroup(fit),
        direct      = .diagnose_direct(fit),
        list()
    )

    result <- c(
        list(
            gof              = gof,
            residuals        = resid,
            mahalanobis      = mahal,
            r_squared        = r2,
            n_outliers       = sum(mahal$outlier),
            largest_residual = resid$summary[1, , drop = FALSE],
            model_type       = fit$model_type
        ),
        specific
    )

    class(result) <- "biomimic_diagnostics"
    result
}


# ===========================================================================
# Print methods
# ===========================================================================

#' @export
print.biomimic_hypotheses <- function(x, ...) {
    cat(sprintf("=== Hypothesis Tests (%s model) ===\n\n",
                toupper(x$model_type)))

    # Invariance hierarchy (multigroup)
    if (!is.null(x$invariance)) {
        print(x$invariance)
        return(invisible(x))
    }

    # Individual tests
    if (!is.null(x$tests) && nrow(x$tests) > 0) {
        cat("--- Individual Tests ---\n")
        for (i in seq_len(nrow(x$tests))) {
            row <- x$tests[i, ]
            pval <- row$pvalue
            sig <- if (is.na(pval)) "" else .signif_stars(pval)
            verdict <- if (is.na(pval)) "p-value not available"
                       else if (pval < x$level) "REJECT H0"
                       else "Do not reject H0"
            cat(sprintf("  %s\n", row$hypothesis))
            cat(sprintf("    est = %.3f, z = %.2f, p = %s %s -> %s\n",
                        row$estimate,
                        if (is.na(row$z)) NA_real_ else row$z,
                        .fmt_pval(pval), sig, verdict))
        }
        cat("\n")
    }

    # Omnibus test
    if (!is.null(x$omnibus)) {
        cat("--- Omnibus Test ---\n")
        om <- x$omnibus
        cat(sprintf("  %s\n", om$hypothesis))
        cat(sprintf("  Chi2 = %.2f, df = %d, p = %s -> %s\n",
                    om$chisq, om$df, .fmt_pval(om$pvalue),
                    if (om$pvalue < x$level) "REJECT" else "Do not reject"))
        cat("\n")
    }

    # Mediation decomposition (direct effects model)
    if (!is.null(x$mediation)) {
        cat("--- Effect Decomposition ---\n")
        med <- x$mediation
        for (i in seq_len(nrow(med))) {
            row <- med[i, ]
            cov_lab <- if (!is.null(row$covariate)) {
                sprintf(" [%s]", row$covariate)
            } else ""
            cat(sprintf("  %s%s:\n", row$protein, cov_lab))
            cat(sprintf("    Indirect (via pathway): %.3f\n", row$indirect))
            cat(sprintf("    Direct:                 %.3f (p = %s)\n",
                        row$direct, .fmt_pval(row$direct_p)))
            pct <- if (abs(row$indirect) + abs(row$direct) > 0) {
                100 * abs(row$indirect) / (abs(row$indirect) + abs(row$direct))
            } else 0
            cat(sprintf("    Proportion mediated:    %.0f%%\n", pct))
        }
        cat("\n")
    }

    cat(sprintf("Significance level: %.2f\n", x$level))
    invisible(x)
}


#' @export
print.biomimic_diagnostics <- function(x, ...) {
    cat(sprintf("=== Model Diagnostics (%s) ===\n\n", toupper(x$model_type)))

    # GOF summary
    cat("--- Goodness of Fit ---\n")
    gof <- x$gof
    cat(sprintf("  RMSEA = %.3f %s | CFI = %.3f %s | SRMR = %.3f %s\n",
                gof$rmsea, .interpret_rmsea(gof$rmsea),
                gof$cfi, .interpret_cfi(gof$cfi),
                gof$srmr, .interpret_srmr(gof$srmr)))
    cat("\n")

    # R-squared
    if (!is.null(x$r_squared)) {
        cat("--- Variance Explained per Protein (R-squared) ---\n")
        r2 <- sort(x$r_squared, decreasing = TRUE)
        n_show <- min(10, length(r2))
        for (i in seq_len(n_show)) {
            bar <- paste(rep("|", round(r2[i] * 20)), collapse = "")
            cat(sprintf("  %-15s %.3f %s\n", names(r2)[i], r2[i], bar))
        }
        if (length(r2) > n_show) cat(sprintf("  ... and %d more\n", length(r2) - n_show))
        cat(sprintf("  Median R2 = %.3f\n\n", stats::median(r2)))
    }

    # Residuals
    cat("--- Largest Residual Correlations ---\n")
    n_show <- min(5, nrow(x$residuals$summary))
    for (i in seq_len(n_show)) {
        row <- x$residuals$summary[i, ]
        flag <- if (abs(row$residual) > 0.1) " [!]" else ""
        cat(sprintf("  %s <-> %s: %.3f%s\n",
                    row$var1_label, row$var2_label, row$residual, flag))
    }
    cat("\n")

    # Outliers
    cat(sprintf("--- Outlier Detection ---\n  %d outliers detected (Mahalanobis, p < 0.001)\n",
                x$n_outliers))
    if (x$n_outliers > 0 && x$n_outliers <= 10) {
        outlier_rows <- x$mahalanobis[x$mahalanobis$outlier, ]
        cat(sprintf("  Observations: %s\n",
                    paste(outlier_rows$observation, collapse = ", ")))
    }
    cat("\n")

    # Model-specific
    if (!is.null(x$reliability)) {
        cat(sprintf("--- Reliability ---\n  Omega = %.3f\n\n", x$reliability))
    }

    if (!is.null(x$congruence)) {
        cat(sprintf("--- Factor Congruence (Tucker's phi) ---\n  phi = %.3f %s\n\n",
                    x$congruence,
                    if (x$congruence >= 0.95) "[Essentially identical]"
                    else if (x$congruence >= 0.85) "[Fair similarity]"
                    else "[Poor similarity]"))
    }

    if (!is.null(x$direct_vs_nodirect)) {
        cat("--- Direct Effects Necessity ---\n")
        lrt <- x$direct_vs_nodirect
        cat(sprintf("  LRT (with vs. without direct effects): Chi2 = %.2f, df = %d, p = %s\n",
                    lrt$chisq_diff, lrt$df_diff, .fmt_pval(lrt$pvalue)))
        cat(sprintf("  -> %s\n\n",
                    if (lrt$pvalue < 0.05) "Direct effects improve the model significantly"
                    else "Direct effects are not necessary"))
    }

    invisible(x)
}


# ===========================================================================
# Model-specific test implementations
# ===========================================================================

#' @noRd
.test_mimic <- function(fit, level) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, ci = TRUE, level = 1 - level)
    struct <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]

    tests <- data.frame(
        hypothesis = paste0("H0: effect of '", struct$rhs, "' on ",
                            .bio_label(struct$lhs, .factor_label_map(fit)),
                            " = 0"),
        parameter  = paste(struct$lhs, "~", struct$rhs),
        estimate   = struct$est,
        se         = struct$se,
        z          = struct$z,
        pvalue     = struct$pvalue,
        ci_lower   = struct$ci.lower,
        ci_upper   = struct$ci.upper,
        significant = struct$pvalue < level,
        stringsAsFactors = FALSE
    )

    # Omnibus: joint test of all structural paths
    omnibus <- NULL
    if (nrow(struct) > 1) {
        constraints <- paste0(struct$lhs, " ~ ", struct$rhs, " == 0")
        omnibus <- tryCatch({
            wt <- lavaan::lavTestWald(fit$lavaan_fit, constraints = constraints)
            list(hypothesis = "H0: all covariate effects on latent pathway = 0",
                 chisq = wt$stat, df = wt$df, pvalue = wt$p.value)
        }, error = function(e) NULL)
    }

    list(tests = tests, omnibus = omnibus, mediation = NULL, invariance = NULL)
}


#' @noRd
.test_multi_mimic <- function(fit, level) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, ci = TRUE, level = 1 - level)
    struct <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]

    tests <- data.frame(
        hypothesis = paste0("H0: effect of '", struct$rhs, "' on ",
                            .bio_label(struct$lhs, .factor_label_map(fit)),
                            " = 0"),
        parameter  = paste(struct$lhs, "~", struct$rhs),
        estimate   = struct$est,
        se         = struct$se,
        z          = struct$z,
        pvalue     = struct$pvalue,
        ci_lower   = struct$ci.lower,
        ci_upper   = struct$ci.upper,
        significant = struct$pvalue < level,
        stringsAsFactors = FALSE
    )

    # Omnibus across all factors
    omnibus <- NULL
    if (nrow(struct) > 1) {
        constraints <- paste0(struct$lhs, " ~ ", struct$rhs, " == 0")
        omnibus <- tryCatch({
            wt <- lavaan::lavTestWald(fit$lavaan_fit, constraints = constraints)
            list(hypothesis = "H0: all covariate effects across all pathways = 0",
                 chisq = wt$stat, df = wt$df, pvalue = wt$p.value)
        }, error = function(e) NULL)
    }

    # Factor correlations
    factor_cors <- pe[pe$op == "~~" & pe$lhs != pe$rhs &
                      grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]

    list(tests = tests, omnibus = omnibus, factor_cors = factor_cors,
         mediation = NULL, invariance = NULL)
}


#' @noRd
.test_multigroup <- function(fit, level) {
    # For multigroup, the main test IS the invariance hierarchy
    inv <- test_invariance(fit)
    list(tests = NULL, omnibus = NULL, mediation = NULL, invariance = inv)
}


#' @noRd
.test_direct <- function(fit, level) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, ci = TRUE, level = 1 - level)

    # Structural paths (via latent)
    struct <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]

    # Direct effect paths
    direct <- pe[pe$op == "~" & !grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]

    # All tests
    all_tests <- rbind(
        data.frame(
            hypothesis = paste0("H0: latent effect of '", struct$rhs, "' on ",
                                .bio_label(struct$lhs, .factor_label_map(fit)),
                                " = 0"),
            parameter  = paste(struct$lhs, "~", struct$rhs),
            estimate   = struct$est,
            se         = struct$se,
            z          = struct$z,
            pvalue     = struct$pvalue,
            ci_lower   = struct$ci.lower,
            ci_upper   = struct$ci.upper,
            significant = struct$pvalue < level,
            stringsAsFactors = FALSE
        ),
        if (nrow(direct) > 0) data.frame(
            hypothesis = paste0("H0: direct effect on '",
                                .bio_label(direct$lhs, fit$indicator_labels),
                                "' from '", direct$rhs, "' = 0"),
            parameter  = paste(direct$lhs, "~", direct$rhs),
            estimate   = direct$est,
            se         = direct$se,
            z          = direct$z,
            pvalue     = direct$pvalue,
            ci_lower   = direct$ci.lower,
            ci_upper   = direct$ci.upper,
            significant = direct$pvalue < level,
            stringsAsFactors = FALSE
        ) else NULL
    )

    # Mediation decomposition. Pair the indirect path (B_{F~c} * lambda_v) and
    # the direct path (gamma_{v,c}) by the SAME covariate c (not by position),
    # so each row is the decomposition for one (indicator, covariate) pair --
    # in particular the group covariate, as in eq:mediation.
    mediation <- NULL
    if (nrow(direct) > 0 && nrow(struct) > 0) {
        loadings <- pe[pe$op == "=~", , drop = FALSE]
        med_rows <- list()
        for (di in seq_len(nrow(direct))) {
            v     <- direct$lhs[di]
            cov_c <- direct$rhs[di]
            lam   <- loadings[loadings$rhs == v, , drop = FALSE]
            st    <- struct[struct$rhs == cov_c, , drop = FALSE]
            if (nrow(lam) == 0 || nrow(st) == 0) next

            indirect <- st$est[1] * lam$est[1]   # B_{F~c} * lambda_v
            direct_e <- direct$est[di]           # gamma_{v,c}

            med_rows[[length(med_rows) + 1]] <- data.frame(
                protein    = .bio_label(v, fit$indicator_labels),
                variable   = v,
                covariate  = cov_c,
                indirect   = indirect,
                direct     = direct_e,
                direct_p   = direct$pvalue[di],
                total      = indirect + direct_e,
                stringsAsFactors = FALSE
            )
        }
        if (length(med_rows) > 0) {
            mediation <- do.call(rbind, med_rows)
        }
    }

    list(tests = all_tests, omnibus = NULL, mediation = mediation,
         invariance = NULL)
}


# ===========================================================================
# Model-specific diagnostic implementations
# ===========================================================================

#' @noRd
.diagnose_mimic <- function(fit) {
    # McDonald's omega reliability. Use the STANDARDIZED solution so the
    # (endogenous) factor variance is normalised to 1 and the marker-scaling
    # of the loadings is removed: omega = (sum l_std)^2 /
    # ((sum l_std)^2 + sum theta_std). Using raw marker-scaled loadings and
    # ignoring the factor disturbance variance (implicitly phi=1) is wrong for
    # a MIMIC, where the factor is endogenous (F ~ covariates).
    reliability <- tryCatch({
        ss <- lavaan::standardizedSolution(fit$lavaan_fit)
        loadings  <- ss[ss$op == "=~", "est.std"]
        residuals <- ss[ss$op == "~~" & ss$lhs == ss$rhs &
                        !grepl("^F", ss$lhs), "est.std"]
        sum_lam <- sum(loadings)
        sum_lam^2 / (sum_lam^2 + sum(residuals))
    }, error = function(e) NULL)

    list(reliability = reliability)
}


#' @noRd
.diagnose_multi_mimic <- function(fit) {
    .diagnose_mimic(fit)  # Same diagnostics apply
}


#' @noRd
.diagnose_multigroup <- function(fit) {
    # Tucker's congruence coefficient between group loading vectors
    congruence <- tryCatch({
        lam_list <- lavaan::lavInspect(fit$lavaan_fit, "est")
        if (is.list(lam_list) && length(lam_list) >= 2) {
            lam1 <- lam_list[[1]]$lambda
            lam2 <- lam_list[[2]]$lambda
            # Tucker's phi = sum(l1*l2) / sqrt(sum(l1^2) * sum(l2^2))
            sum(lam1 * lam2) / sqrt(sum(lam1^2) * sum(lam2^2))
        } else NULL
    }, error = function(e) NULL)

    list(congruence = congruence)
}


#' @noRd
.diagnose_direct <- function(fit) {
    base_diag <- .diagnose_mimic(fit)

    # LRT: model with direct effects vs. without
    direct_vs_nodirect <- tryCatch({
        # Fit the restricted model (no direct effects)
        sel <- fit$selection
        model_nodirect <- build_lavaan_model(sel, type = "mimic")
        fit_nodirect <- lavaan::sem(model_nodirect, data = sel$data,
                                     estimator = fit$lavaan_fit@Options$estimator)
        # Under a scaled estimator (MLR), lavTestLRT automatically returns the
        # Satorra-Bentler scaled chi-square difference, so "Chisq diff" below is
        # the correctly-scaled statistic (not the naive difference).
        lrt <- lavaan::lavTestLRT(fit_nodirect, fit$lavaan_fit)
        if (nrow(lrt) >= 2) {
            list(chisq_diff = lrt[["Chisq diff"]][2],
                 df_diff    = lrt[["Df diff"]][2],
                 pvalue     = lrt[["Pr(>Chisq)"]][2])
        } else NULL
    }, error = function(e) NULL)

    c(base_diag, list(direct_vs_nodirect = direct_vs_nodirect))
}
