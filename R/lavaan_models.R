# --------------------------------------------------------------------------
# lavaan Model Specification Generators
# Stage 3 of the biomimic pipeline
# --------------------------------------------------------------------------

#' Build lavaan Model Syntax from Selection
#'
#' Generates lavaan model syntax from the output of Stage 2. Supports
#' four model types: MIMIC, multi-factor MIMIC, multi-group CFA, and
#' MIMIC with direct effects.
#'
#' @section Identification:
#' Each model type requires specific conditions for parameter identification.
#' The function enforces these automatically:
#'
#' \describe{
#'   \item{MIMIC (K=1)}{Scale set via marker variable (first loading fixed
#'     to 1 by lavaan). Requires p >= 3 indicators. Identified when residual
#'     covariance is diagonal (lavaan default).
#'     Reference: Stapleton (1978).}
#'   \item{Multi-factor MIMIC (K>=2)}{Simple structure (block-diagonal Lambda)
#'     from Stage 2. Each factor needs >= 3 indicators for robust
#'     identification; with only 2 indicators, identification requires nonzero
#'     factor covariance (phi_kl != 0). Factor covariances are freely estimated
#'     (lavaan default).
#'     Reference: Reilly (1995), Marques (2018, Theorem 2).}
#'   \item{Multi-group CFA}{Configural model identified per group (marker
#'     variable). Equality constraints (metric, scalar, strict) cannot hinder
#'     identification of an already-identified model (Reilly, 1995,
#'     Observation 1). Reference: Marques (2018, Theorem 1).}
#'   \item{MIMIC + direct effects}{Same as MIMIC K=1 plus freely estimated
#'     direct paths. Identified when residual covariance is diagonal.
#'     Reference: Stapleton (1978).}
#' }
#'
#' @param selection A `biomimic_selection` object from \code{top_k()}.
#' @param type Character, model type:
#'   \describe{
#'     \item{`"mimic"`}{MIMIC model: latent factor(s) with covariate effects
#'       via structural regression. Default.}
#'     \item{`"multigroup"`}{Multi-group CFA for measurement invariance
#'       testing. Requires `group_var`.}
#'     \item{`"direct"`}{MIMIC with direct effects: adds X -> Y paths for
#'       variables flagged by SNR_Gamma in Stage 1.}
#'     \item{`"multi_mimic"`}{Multi-factor MIMIC (K >= 2) with simple
#'       structure from \code{discover_factors()}.}
#'   }
#' @param group_var Character, name of grouping variable for multi-group
#'   models (required when `type = "multigroup"`).
#' @param covariate_names Character vector of covariate names to include in
#'   structural regression (for MIMIC models). If NULL, uses all columns
#'   of X except the intercept.
#'
#' @return A character string containing lavaan model syntax.
#'
#' @references
#' Reilly, T. (1995). A necessary and sufficient condition for identification
#' of confirmatory factor analysis models of factor complexity one.
#' \emph{Sociological Methods & Research}, 23(4), 421-441.
#'
#' Stapleton, D. C. (1978). Analyzing political participation data with a
#' MIMIC model. \emph{Sociological Methodology}, 9, 52-74.
#'
#' Marques, A. H. C. (2018). Multiple factor analysis model with scale
#' mixture of normal distributions in the latent factors.
#' MSc dissertation, Universidade Federal de Pernambuco.
#'
#' @examples
#' \dontrun{
#' sel <- top_k(vb_fit, k = 15)
#'
#' # MIMIC model
#' model <- build_lavaan_model(sel, type = "mimic")
#' fit <- lavaan::sem(model, data = sel$data, estimator = "MLR")
#'
#' # Multi-group CFA
#' model_mg <- build_lavaan_model(sel, type = "multigroup", group_var = "group")
#' fit_mg <- lavaan::sem(model_mg, data = sel$data, group = "group")
#'
#' # MIMIC with direct effects
#' model_de <- build_lavaan_model(sel, type = "direct")
#' fit_de <- lavaan::sem(model_de, data = sel$data, estimator = "MLR")
#' }
#'
#' @export
build_lavaan_model <- function(selection,
                               type = c("mimic", "multigroup", "direct",
                                        "multi_mimic"),
                               group_var = NULL,
                               covariate_names = NULL) {
    type <- match.arg(type)
    stopifnot(inherits(selection, "biomimic_selection"))

    # --- Identification checks ---
    .check_identification_preconditions(selection, type)

    switch(type,
        mimic       = .build_mimic(selection, covariate_names),
        multigroup  = .build_multigroup(selection, group_var),
        direct      = .build_direct(selection, covariate_names),
        multi_mimic = .build_multi_mimic(selection, covariate_names)
    )
}


#' Run the Full biomimic Pipeline
#'
#' Convenience function that runs all three stages: VB-ARD screening,
#' top-k selection, and lavaan estimation. Returns a \code{biomimic_fit}
#' object with the VB-ARD fit attached.
#'
#' @param Y Numeric matrix (n x p) of observed variables.
#' @param X Numeric matrix (n x q) of covariates (should include intercept).
#' @param K Integer, number of latent factors (default 1).
#' @param k Integer, panel size for variable selection (default 15).
#' @param model_type Character, Stage 3 model type (default \code{"mimic"}).
#' @param estimator Character, lavaan estimator (default \code{"MLR"}).
#' @param group_var Character, grouping variable (for multigroup models).
#' @param indicator_labels Named character vector mapping variable names
#'   to biology-friendly labels. If \code{NULL}, original names are used.
#' @param factor_labels Character vector of factor names. If \code{NULL},
#'   defaults are generated.
#' @param gamma_anchor Identification anchors for the \code{B}/\code{Gamma}
#'   split. \code{"auto"} (the default) runs Stage 1 twice: a first pass to
#'   rank indicators by loading SNR, then \code{\link{select_gamma_anchor}} to
#'   pick \code{K} anchors, then a second pass with their direct effects held
#'   at zero. \code{NULL} skips anchoring and reproduces the unanchored,
#'   prior-regularised model, whose \code{B}/\code{Gamma} split is not
#'   identified by the likelihood (see \code{\link{vb_ard}}). A character or
#'   integer vector supplies the anchors directly, in one pass.
#' @param gamma_threshold SNR threshold for flagging direct effects, passed to
#'   \code{\link{top_k}}: a numeric value (default 3.0) or
#'   \code{"bootstrap"} for a parametric-bootstrap calibration that controls
#'   the family-wise error rate at 5 percent over the selected panel
#'   (see \code{\link{calibrate_gamma_threshold}}; adds \code{n_boot} refits
#'   of the screening model to the runtime).
#' @param ... Additional arguments passed to \code{vb_ard()}.
#'
#' @return A \code{biomimic_fit} object. The VB-ARD fit is accessible via
#'   \code{$selection$fit}.
#'
#' @export
biomimic <- function(Y, X, K = 1L, k = 15L,
                     model_type = "mimic",
                     estimator = "MLR",
                     group_var = NULL,
                     indicator_labels = NULL,
                     factor_labels = NULL,
                     gamma_anchor = "auto",
                     gamma_threshold = 3.0, ...) {
    # Stage 1
    if (identical(gamma_anchor, "auto")) {
        pass1   <- vb_ard(Y, X, K = K, ...)
        anchors <- select_gamma_anchor(pass1, n_anchor = K)
        vb_fit  <- vb_ard(Y, X, K = K, gamma_anchor = anchors, ...)
    } else {
        vb_fit <- vb_ard(Y, X, K = K, gamma_anchor = gamma_anchor, ...)
    }

    # Stage 2
    selection <- top_k(vb_fit, k = k, gamma_threshold = gamma_threshold)

    # Stage 3
    fit_sem(selection, type = model_type, estimator = estimator,
            group_var = group_var,
            indicator_labels = indicator_labels,
            factor_labels = factor_labels)
}


# --- Internal model builders ------------------------------------------------

#' @noRd
.build_mimic <- function(selection, covariate_names) {
    # Build MIMIC syntax: F =~ y1 + y2 + ... ; F ~ x1 + x2 + ...
    vars <- selection$variables
    covs <- covariate_names %||% .get_covariate_names(selection)

    measurement <- paste0("F1 =~ ", paste(vars, collapse = " + "))
    structural  <- paste0("F1 ~ ", paste(covs, collapse = " + "))

    paste(measurement, structural, sep = "\n")
}

#' @noRd
.build_multigroup <- function(selection, group_var) {
    if (is.null(group_var)) {
        stop("group_var is required for multigroup CFA models.")
    }
    # Build CFA syntax (no structural part -- groups handled by lavaan)
    vars <- selection$variables
    paste0("F1 =~ ", paste(vars, collapse = " + "))
}

#' @noRd
.build_direct <- function(selection, covariate_names) {
    # Build MIMIC + direct effects syntax
    vars <- selection$variables
    covs <- covariate_names %||% .get_covariate_names(selection)
    direct_vars <- intersect(selection$direct_effects, vars)

    # Identification condition (iii): the marker indicator (first in the
    # measurement line, whose loading lavaan fixes to 1) must NOT carry a
    # direct effect. Reorder so a pure (non-direct) indicator is the marker.
    pure <- setdiff(vars, direct_vars)
    if (length(pure) > 0 && vars[1] %in% direct_vars) {
        vars <- c(pure[1], setdiff(vars, pure[1]))
    }

    lines <- character()
    # Measurement model
    lines <- c(lines, paste0("F1 =~ ", paste(vars, collapse = " + ")))
    # Structural model
    lines <- c(lines, paste0("F1 ~ ", paste(covs, collapse = " + ")))
    # Direct effects (indicator regressed on the covariates)
    for (v in direct_vars) {
        lines <- c(lines, paste0(v, " ~ ", paste(covs, collapse = " + ")))
    }

    paste(lines, collapse = "\n")
}

#' @noRd
.build_multi_mimic <- function(selection, covariate_names) {
    # Build multi-factor MIMIC with simple structure
    groups <- selection$indicator_groups
    covs <- covariate_names %||% .get_covariate_names(selection)
    vars <- selection$variables
    K <- length(groups)

    lines <- character()
    # Measurement model: one line per factor
    for (k in seq_len(K)) {
        factor_vars <- vars[groups[[k]]]
        lines <- c(lines,
                   paste0("F", k, " =~ ", paste(factor_vars, collapse = " + ")))
    }
    # Structural model: each factor regressed on covariates
    for (k in seq_len(K)) {
        lines <- c(lines,
                   paste0("F", k, " ~ ", paste(covs, collapse = " + ")))
    }

    paste(lines, collapse = "\n")
}

#' @noRd
.get_covariate_names <- function(selection) {
    # Extract covariate names from X, excluding intercept
    xnames <- colnames(selection$fit$X)
    if (!is.null(xnames)) {
        xnames <- xnames[!tolower(xnames) %in% c("intercept", "(intercept)")]
    }
    xnames
}


# ===========================================================================
# Identification checks
# ===========================================================================

#' @noRd
.check_identification_preconditions <- function(selection, type) {
    vars <- selection$variables
    p <- length(vars)

    .log_debug("check_identification: type=%s, p=%d indicators", type, p)

    # All models: need at least 3 indicators for K=1
    if (type %in% c("mimic", "multigroup", "direct") && p < 3) {
        stop(sprintf(
            "Model requires at least 3 indicators for identification. Got %d.\n",
            p),
            "  Increase k in top_k() or add more variables.",
            call. = FALSE)
    }

    # MIMIC + direct effects: B is identified only if the pure indicators (no
    # direct path) have loadings of rank K, which needs at least K of them.
    # We keep the stricter |S| <= p - 2 as the operating rule: it equals the
    # rank requirement at K = 2 and is deliberately conservative at K = 1,
    # where a single pure indicator would suffice in principle but leaves no
    # slack if that indicator turns out to carry a direct effect after all.
    if (type == "direct") {
        n_direct <- length(intersect(selection$direct_effects, vars))
        n_pure   <- p - n_direct
        K_sel    <- selection$K %||% 1L
        if (n_pure < K_sel) {
            warning(sprintf(
                paste0("Direct-effect model: %d of %d indicators carry a ",
                       "direct effect, leaving %d pure indicator(s) for %d ",
                       "factor(s). Identification of B requires at least K ",
                       "pure indicators with linearly independent loadings; ",
                       "the model is under-identified."),
                n_direct, p, n_pure, K_sel), call. = FALSE)
        } else if (n_pure < 2) {
            warning(sprintf(
                paste0("Direct-effect model: only %d pure indicator(s) ",
                       "remain (|S| = %d of %d). This meets the bare rank ",
                       "condition but leaves no redundancy; >= 2 is ",
                       "recommended."),
                n_pure, n_direct, p), call. = FALSE)
        }
    }

    # Multi-factor MIMIC: check each factor has enough indicators
    if (type == "multi_mimic") {
        groups <- selection$indicator_groups
        if (is.null(groups)) {
            stop("indicator_groups is required for multi_mimic models.\n",
                 "  Fit with K >= 2 or use discover_factors().",
                 call. = FALSE)
        }
        K <- length(groups)
        for (k in seq_len(K)) {
            n_ind <- length(groups[[k]])
            if (n_ind < 2) {
                stop(sprintf(
                    "Factor %d has %d indicator(s). Minimum is 2.\n", k, n_ind),
                    "  Rerun top_k() with a larger k or adjust K.",
                    call. = FALSE)
            }
            if (n_ind == 2) {
                warning(sprintf(
                    paste0("Factor %d has only 2 indicators. ",
                           "Identification requires nonzero factor ",
                           "covariance (phi_kl != 0).\n",
                           "  Consider using >= 3 indicators per factor ",
                           "for robust identification.\n",
                           "  Reference: Reilly (1995), Example 2."),
                    k), call. = FALSE)
            }
        }
        # Total indicators must be at least 2*K + 1 for over-identification
        min_total <- 2 * K + 1
        if (p < min_total) {
            warning(sprintf(
                paste0("Total indicators (%d) is marginal for K=%d factors. ",
                       "The model may be just-identified or ",
                       "under-identified.\n",
                       "  Recommended: at least %d indicators for %d factors."),
                p, K, 3 * K, K), call. = FALSE)
        }
    }
}


#' Check Identification Status of a Fitted Model
#'
#' Evaluates whether the model parameters are identified by examining
#' the degrees of freedom, the information matrix rank, and model-specific
#' conditions. This function should be called after fitting to verify
#' that the estimation is meaningful.
#'
#' @section Identification Rules:
#' \describe{
#'   \item{Necessary condition (t-Rule)}{The number of free parameters must
#'     not exceed p(p+1)/2, the number of non-redundant elements in the
#'     observed covariance matrix. Equivalently, degrees of freedom >= 0.}
#'   \item{MIMIC K=1}{Identified with p >= 3 indicators and marker variable
#'     normalization (Stapleton, 1978).}
#'   \item{Multi-factor MIMIC K>=2}{Simple structure (block-diagonal Lambda)
#'     with >= 3 indicators per factor guarantees identification. With 2
#'     indicators per factor, requires nonzero factor covariance
#'     (Reilly, 1995; Marques, 2018).}
#'   \item{Multi-group CFA}{Each group must be separately identified.
#'     Equality constraints (metric, scalar) preserve identification
#'     (Reilly, 1995, Observation 1).}
#'   \item{Rank Rule (Reilly, 1995)}{For factor complexity 1 models,
#'     construct binary matrix R from normal equations and verify
#'     rank(R) equals the number of free parameters. This is a necessary
#'     and sufficient condition.}
#' }
#'
#' @param fit A \code{biomimic_fit} object.
#'
#' @return A list of class \code{biomimic_identification} with components:
#' \describe{
#'   \item{identified}{Logical, whether the model appears identified}
#'   \item{df}{Degrees of freedom (>= 0 required)}
#'   \item{n_free}{Number of free parameters}
#'   \item{n_equations}{Number of non-redundant covariance elements p(p+1)/2}
#'   \item{converged}{Whether lavaan converged}
#'   \item{information_rank}{Rank of the information matrix (if available)}
#'   \item{warnings}{Character vector of identification warnings}
#'   \item{model_type}{Character, model type}
#' }
#'
#' @export
check_identification <- function(fit) {
    stopifnot(inherits(fit, "biomimic_fit"))

    lav <- fit$lavaan_fit
    warnings <- character()
    identified <- TRUE

    # Basic counts
    pe <- lavaan::parameterEstimates(lav)
    n_free <- lavaan::fitMeasures(lav, "npar")
    df <- lavaan::fitMeasures(lav, "df")
    p <- length(fit$selection$variables)
    n_equations <- p * (p + 1) / 2

    converged <- lavaan::lavInspect(lav, "converged")

    # t-Rule: necessary condition
    if (df < 0) {
        identified <- FALSE
        warnings <- c(warnings, sprintf(
            "FAIL: Negative degrees of freedom (%d). Model has more free parameters (%d) than equations (%d).",
            df, n_free, n_equations))
    } else if (df == 0) {
        warnings <- c(warnings,
            "NOTE: Model is just-identified (df = 0). No overidentifying restrictions to test.")
    }

    # Convergence check
    if (!converged) {
        identified <- FALSE
        warnings <- c(warnings,
            "FAIL: lavaan did not converge. This may indicate non-identification or numerical issues.")
    }

    # Information matrix rank (if available)
    info_rank <- NULL
    info_check <- tryCatch({
        vcov_mat <- lavaan::vcov(lav)
        info_rank <- qr(vcov_mat)$rank
        n_free_int <- as.integer(n_free)
        if (info_rank < n_free_int) {
            identified <- FALSE
            warnings <- c(warnings, sprintf(
                "FAIL: Information matrix rank (%d) < number of free parameters (%d). Some parameters are not identified.",
                info_rank, n_free_int))
        }
        TRUE
    }, error = function(e) {
        warnings <<- c(warnings, paste(
            "WARN: Could not compute information matrix rank:", e$message))
        FALSE
    })

    # Model-specific checks
    mt <- fit$model_type

    if (mt == "mimic" || mt == "direct") {
        if (p < 3) {
            identified <- FALSE
            warnings <- c(warnings,
                "FAIL: MIMIC K=1 requires at least 3 indicators (Reilly, 1995).")
        }
    }

    if (mt == "multi_mimic") {
        groups <- fit$selection$indicator_groups
        if (!is.null(groups)) {
            for (k in seq_along(groups)) {
                n_ind <- length(groups[[k]])
                if (n_ind < 2) {
                    identified <- FALSE
                    warnings <- c(warnings, sprintf(
                        "FAIL: Factor %d has %d indicator(s). Minimum is 2.", k, n_ind))
                } else if (n_ind == 2) {
                    warnings <- c(warnings, sprintf(
                        "WARN: Factor %d has only 2 indicators. Requires nonzero factor covariance for identification (Reilly, 1995).",
                        k))
                }
            }
        }
    }

    # Negative variance estimates (Heywood cases)
    resid_var <- pe[pe$op == "~~" & pe$lhs == pe$rhs, ]
    neg_var <- resid_var[resid_var$est < 0, ]
    if (nrow(neg_var) > 0) {
        warnings <- c(warnings, sprintf(
            "WARN: %d negative variance estimate(s) detected (Heywood case). Variables: %s. This may indicate misspecification or non-identification.",
            nrow(neg_var), paste(neg_var$lhs, collapse = ", ")))
    }

    if (length(warnings) == 0) {
        warnings <- "All identification checks passed."
    }

    result <- list(
        identified       = identified,
        df               = as.integer(df),
        n_free           = as.integer(n_free),
        n_equations      = as.integer(n_equations),
        converged        = converged,
        information_rank = info_rank,
        warnings         = warnings,
        model_type       = mt
    )
    class(result) <- "biomimic_identification"
    result
}


#' @export
print.biomimic_identification <- function(x, ...) {
    cat("=== Identification Check ===\n")
    cat(sprintf("Model type: %s\n", x$model_type))
    cat(sprintf("Free parameters: %d\n", x$n_free))
    cat(sprintf("Covariance equations: %d\n", x$n_equations))
    cat(sprintf("Degrees of freedom: %d\n", x$df))
    cat(sprintf("Converged: %s\n", x$converged))
    if (!is.null(x$information_rank)) {
        cat(sprintf("Information matrix rank: %d / %d\n",
                    x$information_rank, x$n_free))
    }
    cat(sprintf("\nStatus: %s\n\n",
                if (x$identified) "IDENTIFIED" else "NOT IDENTIFIED"))
    for (w in x$warnings) {
        cat(sprintf("  %s\n", w))
    }
    invisible(x)
}
