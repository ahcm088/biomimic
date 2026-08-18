# --------------------------------------------------------------------------
# Variable Selection and Factor Structure Discovery
# Stage 2 of the biomimic pipeline
# --------------------------------------------------------------------------

#' Rank Variables by Signal-to-Noise Ratio
#'
#' Extracts and ranks variables by their posterior SNR from a fitted VB-ARD
#' model. For K=1, returns a single ranking. For K>1, returns the maximum
#' SNR across factors for each variable.
#'
#' @param fit A fitted \code{vb_ard} object.
#' @param type Character, which SNR to rank by: \code{"lambda"} (loading signal,
#'   default) or \code{"gamma"} (direct effect signal).
#'
#' @return A data.frame with columns: \code{variable}, \code{snr}, \code{rank},
#'   sorted by descending SNR. For K>1, additional columns \code{snr_k1},
#'   \code{snr_k2}, etc.
#'
#' @export
rank_variables <- function(fit, type = c("lambda", "gamma")) {
    type <- match.arg(type)
    stopifnot(inherits(fit, "vb_ard"))

    if (type == "lambda") {
        snr_mat <- fit$snr_lambda
    } else {
        if (is.null(fit$snr_gamma)) {
            stop("No Gamma SNR available. Fit with use_gamma=TRUE.")
        }
        snr_mat <- fit$snr_gamma
        # Exclude intercept / constant covariate columns: the intercept is not
        # a "direct effect" candidate (and with centered Y its coefficient is
        # ~0), so it must not enter the direct-effect ranking.
        const_cols <- which(apply(fit$X, 2,
                                  function(col) length(unique(col)) == 1L))
        if (length(const_cols) > 0 && length(const_cols) < ncol(snr_mat)) {
            snr_mat <- snr_mat[, -const_cols, drop = FALSE]
        }
    }

    keep_rows <- rep(TRUE, nrow(snr_mat))
    if (type == "gamma" && length(fit$anchor_idx)) {
        # Anchored indicators have gamma fixed at zero by construction, so they
        # are not direct-effect candidates and carry an undefined (NA) SNR.
        keep_rows[fit$anchor_idx] <- FALSE
        snr_mat <- snr_mat[keep_rows, , drop = FALSE]
    }

    p <- nrow(snr_mat)
    K_or_q <- ncol(snr_mat)

    # Aggregate: max SNR across columns
    snr_max <- apply(snr_mat, 1, max)

    df <- data.frame(
        variable = fit$var_names[keep_rows],
        snr = snr_max,
        stringsAsFactors = FALSE
    )

    # Add per-column SNRs for K > 1
    if (K_or_q > 1) {
        for (k in seq_len(K_or_q)) {
            col_name <- if (type == "lambda") {
                paste0("snr_k", k)
            } else {
                cn <- colnames(snr_mat)[k]
                if (is.null(cn)) paste0("snr_x", k) else paste0("snr_", cn)
            }
            df[[col_name]] <- snr_mat[, k]
        }
    }

    df$rank <- rank(-df$snr, ties.method = "first")
    df <- df[order(df$rank), ]
    rownames(df) <- NULL
    df
}


#' Select Top-k Variables
#'
#' Selects the top-k variables from a VB-ARD fit based on loading SNR,
#' and prepares data for Stage 3.
#'
#' @param fit A fitted \code{vb_ard} object.
#' @param k Integer, number of variables to select (default 15).
#' @param gamma_threshold SNR threshold for flagging direct effects. Either a
#'   numeric value (default 3.0) or \code{"bootstrap"}, which calibrates the
#'   threshold by parametric bootstrap under the fitted model with no direct
#'   effects (see \code{\link{calibrate_gamma_threshold}}); flagging then
#'   controls the family-wise error rate at 5 percent over the selected
#'   panel. The
#'   numeric default is adequate for cohesive panels (it sits at the top of
#'   the calibrated range observed in practice), but it is not
#'   scale-invariant across designs; the bootstrap is the principled choice.
#' @param n_boot Integer, bootstrap refits when
#'   \code{gamma_threshold = "bootstrap"} (default 199).
#' @param n_cores Integer, workers for the bootstrap (default all but one
#'   core; honours \code{BIOMIMIC_NCORES}). Ignored for numeric thresholds.
#'
#' @return A list of class \code{"biomimic_selection"} with components:
#' \describe{
#'   \item{variables}{Character vector of selected variable names}
#'   \item{ranking}{Full ranking data.frame from \code{rank_variables()}}
#'   \item{data}{Data frame with selected Y columns and X covariates}
#'   \item{K}{Number of latent factors}
#'   \item{indicator_groups}{List of indicator groups (for K >= 2), or NULL}
#'   \item{direct_effects}{Character vector of variables with high SNR_Gamma}
#'   \item{gamma_threshold}{The numeric threshold actually applied}
#'   \item{gamma_calibration}{The
#'     \code{\link{calibrate_gamma_threshold}} result, or NULL if a numeric
#'     threshold was supplied}
#'   \item{fit}{The original VB-ARD fit}
#' }
#'
#' @export
top_k <- function(fit, k = 15L, gamma_threshold = 3.0,
                  n_boot = 199L, n_cores = NULL) {
    stopifnot(inherits(fit, "vb_ard"))
    stopifnot(k >= 1L)

    .log_info("top_k: selecting %d variables from %d candidates", k, fit$p)

    ranking <- rank_variables(fit, type = "lambda")
    selected_vars <- ranking$variable[seq_len(min(k, nrow(ranking)))]

    # Anchored indicators must (a) be in the panel and (b) come first, so that
    # Stage 3 fixes its marker loading on an indicator whose direct effect is
    # already constrained to zero. That is condition (iii) of the direct-effect
    # identification result: the marker must be a pure indicator.
    anchors <- fit$gamma_anchor
    if (length(anchors)) {
        missing_anchor <- setdiff(anchors, selected_vars)
        if (length(missing_anchor)) {
            # Make room by dropping the weakest non-anchor indicators.
            drop_n <- length(missing_anchor)
            keepers <- setdiff(selected_vars, anchors)
            keepers <- keepers[seq_len(max(0L, length(keepers) - drop_n))]
            selected_vars <- c(intersect(selected_vars, anchors),
                               missing_anchor, keepers)
            .log_info("top_k: %d anchor(s) forced into the panel", drop_n)
        }
        selected_vars <- c(anchors, setdiff(selected_vars, anchors))
    }

    # Identify direct effects among selected variables
    direct_effects <- character(0)
    gamma_calibration <- NULL
    if (!is.null(fit$snr_gamma)) {
        if (identical(gamma_threshold, "bootstrap")) {
            gamma_calibration <- calibrate_gamma_threshold(
                fit, variables = selected_vars,
                n_boot = n_boot, n_cores = n_cores)
            gamma_threshold <- gamma_calibration$threshold
        }
        stopifnot(is.numeric(gamma_threshold),
                  length(gamma_threshold) == 1L)
        gamma_ranking <- rank_variables(fit, type = "gamma")
        # Only flag selected variables
        sel_gamma <- gamma_ranking[gamma_ranking$variable %in% selected_vars, ]
        direct_effects <- sel_gamma$variable[sel_gamma$snr > gamma_threshold]
        .log_info("top_k: %d direct effects flagged (gamma_threshold=%.2f%s)",
                  length(direct_effects), gamma_threshold,
                  if (is.null(gamma_calibration)) "" else ", bootstrap")
    }

    # Build data frame for lavaan
    sel_idx <- match(selected_vars, fit$var_names)
    Y_sel <- fit$Y[, sel_idx, drop = FALSE]
    colnames(Y_sel) <- selected_vars

    # Combine Y and X (without intercept) into a data frame
    X_df <- as.data.frame(fit$X)
    # Remove intercept column if present
    intercept_cols <- which(apply(fit$X, 2, function(col) all(col == 1)))
    if (length(intercept_cols) > 0) {
        X_df <- X_df[, -intercept_cols, drop = FALSE]
    }

    data <- cbind(as.data.frame(Y_sel), X_df)

    # Discover factor structure if K >= 2
    indicator_groups <- NULL
    if (fit$K >= 2) {
        # Use only selected variables for factor structure
        mu_lam_sel <- fit$mu_lambda[sel_idx, , drop = FALSE]
        indicator_groups <- .suggest_indicator_groups(mu_lam_sel, fit$K)
    }

    structure(
        list(
            variables       = selected_vars,
            ranking         = ranking,
            data            = data,
            K               = fit$K,
            indicator_groups = indicator_groups,
            direct_effects  = direct_effects,
            gamma_threshold = if (is.null(fit$snr_gamma)) NULL
                              else gamma_threshold,
            gamma_calibration = gamma_calibration,
            fit             = fit
        ),
        class = "biomimic_selection"
    )
}


#' Discover Factor Structure via SVD + Varimax
#'
#' For K >= 2, discovers simple structure (block-diagonal Lambda) by
#' applying varimax rotation to the VB-ARD posterior mean loadings.
#' Assigns each variable to the factor on which it loads most strongly.
#'
#' @param fit A fitted \code{vb_ard} object with K >= 2.
#' @param K Integer, number of factors. If NULL, uses \code{fit$K}.
#'
#' @return A list with components:
#' \describe{
#'   \item{indicator_groups}{List of K integer vectors, each containing
#'     column indices of variables assigned to that factor}
#'   \item{rotated_loadings}{Matrix of varimax-rotated loadings (p x K)}
#'   \item{assignments}{Integer vector of factor assignments (length p)}
#' }
#'
#' @export
discover_factors <- function(fit, K = NULL) {
    stopifnot(inherits(fit, "vb_ard"))

    if (is.null(K)) K <- fit$K
    stopifnot(K >= 2L)

    .log_info("discover_factors: K=%d, applying SVD + varimax", K)

    if (ncol(fit$mu_lambda) < K) {
        stop("fit$K < K. Refit with K >= ", K)
    }

    # vb_ard already applied the varimax rotation and stored the rotated
    # loadings (and per-variable factor assignments). Reuse them so the
    # grouping matches the rotated frame used everywhere downstream; only
    # rotate here as a fallback for objects that lack the stored rotation.
    if (!is.null(fit$rotated_lambda) && !is.null(fit$rotation_matrix)) {
        rotated <- fit$rotated_lambda[, seq_len(K), drop = FALSE]
    } else {
        mu_lam <- fit$mu_lambda[, seq_len(K), drop = FALSE]
        rot <- stats::varimax(mu_lam, normalize = TRUE)
        rotated <- mu_lam %*% rot$rotmat
    }

    indicator_groups <- .suggest_indicator_groups(rotated, K)
    assignments <- .get_assignments(rotated, K)

    list(
        indicator_groups = indicator_groups,
        rotated_loadings = rotated,
        assignments      = assignments
    )
}


# ===========================================================================
# Internal helpers
# ===========================================================================

#' @noRd
.suggest_indicator_groups <- function(loadings_matrix, K) {
    # Assign each variable to the factor with highest absolute loading
    assignments <- .get_assignments(loadings_matrix, K)

    groups <- vector("list", K)
    for (k in seq_len(K)) {
        groups[[k]] <- which(assignments == k)
    }

    # Ensure no empty groups
    for (k in seq_len(K)) {
        if (length(groups[[k]]) == 0) {
            warning(sprintf("Factor %d has no assigned variables.", k))
        }
    }

    groups
}

#' @noRd
.get_assignments <- function(loadings_matrix, K) {
    apply(abs(loadings_matrix[, seq_len(K), drop = FALSE]), 1, which.max)
}
