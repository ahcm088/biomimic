# --------------------------------------------------------------------------
# Parametric-bootstrap calibration of the direct-effect SNR threshold
# --------------------------------------------------------------------------

#' Calibrate the Direct-Effect SNR Threshold by Parametric Bootstrap
#'
#' Computes a data-driven threshold for flagging direct effects, replacing the
#' fixed default of \code{top_k()}. Data are simulated from the fitted model
#' under the global null of no direct effects (\eqn{\Gamma = 0}) with the
#' structural coefficients \code{B}, loadings \code{Lambda}, residual
#' variances \code{Psi} and the covariate matrix \code{X} all held at their
#' fitted values; each simulated dataset is refitted and the maximum
#' direct-effect SNR over the panel is recorded. The threshold is the
#' \code{level} quantile of these maxima, so flagging any panel variable above
#' it controls the family-wise error rate at \code{1 - level}.
#'
#' @section Why a parametric bootstrap and not a permutation:
#' Permuting the covariate destroys its structural effect \code{B} along with
#' any direct effects, but the null hypothesis of interest is
#' \eqn{\gamma_j = 0} \emph{with the latent pathway intact}. When the fitted
#' latent pathway under-absorbs a covariate effect (which happens when the
#' loading scale is mis-estimated), the residue leaks into the direct effects
#' of loaded variables in proportion to \code{B} and \eqn{\lambda_j}. A
#' permutation null cannot contain that leakage because it has no \code{B};
#' the parametric bootstrap null does, because \code{B} is kept at its fitted
#' value.
#'
#' @param fit A \code{vb_ard} object fitted with \code{use_gamma = TRUE}.
#' @param variables Character vector: the panel whose family-wise error is to
#'   be controlled (typically the \code{top_k()} selection). Anchored
#'   indicators are dropped automatically (their direct effects are fixed at
#'   zero and carry no SNR).
#' @param n_boot Integer, number of bootstrap refits (default 199).
#' @param level Numeric, quantile of the per-replicate maxima used as the
#'   threshold (default 0.95, controlling FWER at 5 percent).
#' @param n_cores Integer, workers for the PSOCK cluster. Default uses
#'   \code{detectCores() - 1}, honouring the \code{BIOMIMIC_NCORES}
#'   environment variable; set to 1 for serial execution.
#' @param max_iter Integer, CAVI iteration cap for each refit (default 300).
#' @param seed Integer base seed; replicate \code{b} uses \code{seed + b},
#'   so results are reproducible and independent of the worker schedule.
#'
#' @return An object of class \code{"biomimic_gamma_calibration"}: a list with
#'   \code{threshold} (the calibrated value), \code{level}, \code{boot_max}
#'   (the per-replicate maxima), \code{n_boot} (effective replicates),
#'   \code{variables} and \code{columns}.
#'
#' @details Each replicate refits the full VB-ARD model, so the cost is
#'   \code{n_boot} times one \code{vb_ard()} fit; the replicates run in
#'   parallel. The refit reproduces the original fit's settings
#'   (\code{K}, \code{use_ard}, anchors).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- simulate_mimic(n = 150, p_total = 40, p_signal = 15, K = 1,
#'                     has_direct = TRUE)
#' fit <- vb_ard(d$Y, d$X, K = 1)
#' sel <- top_k(fit, k = 15)
#' cal <- calibrate_gamma_threshold(fit, sel$variables,
#'                                  n_boot = 49, n_cores = 1)
#' cal
#' sel2 <- top_k(fit, k = 15, gamma_threshold = cal$threshold)
#' }
#'
#' @seealso \code{\link{top_k}}, \code{\link{vb_ard}}
#' @export
calibrate_gamma_threshold <- function(fit, variables,
                                      n_boot = 199L, level = 0.95,
                                      n_cores = NULL, max_iter = 300L,
                                      seed = 1L) {
    stopifnot(inherits(fit, "vb_ard"))
    if (is.null(fit$snr_gamma)) {
        stop("Fit has no direct effects. Refit with use_gamma = TRUE.",
             call. = FALSE)
    }
    stopifnot(is.character(variables), length(variables) >= 1L)
    miss <- setdiff(variables, fit$var_names)
    if (length(miss)) {
        stop("variables not found in the fit: ",
             paste(miss, collapse = ", "), call. = FALSE)
    }
    if (length(fit$anchor_idx)) {
        variables <- setdiff(variables, fit$var_names[fit$anchor_idx])
        if (!length(variables)) {
            stop("All requested variables are anchors; nothing to calibrate.",
                 call. = FALSE)
        }
    }

    # Non-constant covariate columns, matching the flagging statistic in
    # rank_variables(type = "gamma").
    const <- which(apply(fit$X, 2,
                         function(col) length(unique(col)) == 1L))
    cols <- setdiff(seq_len(fit$q), const)
    if (!length(cols)) {
        stop("No non-constant covariate columns to calibrate against.",
             call. = FALSE)
    }

    X   <- fit$X
    B   <- fit$mu_B
    Lam <- fit$mu_lambda
    psi <- 1 / fit$E_psi_inv
    n <- fit$n; p <- fit$p; K <- fit$K
    use_ard <- fit$use_ard
    anchors <- fit$gamma_anchor
    vnames  <- fit$var_names

    one_boot <- function(b) {
        set.seed(seed + b)
        eta <- matrix(stats::rnorm(n * K), n, K)
        h   <- X %*% B + eta
        Ys  <- h %*% t(Lam) +
               matrix(stats::rnorm(n * p, sd = rep(sqrt(psi), each = n)),
                      n, p)
        colnames(Ys) <- vnames
        f <- tryCatch(
            vb_ard(Ys, X, K = K, use_ard = use_ard, use_gamma = TRUE,
                   gamma_anchor = anchors, max_iter = max_iter),
            error = function(e) NULL)
        if (is.null(f) || is.null(f$snr_gamma)) return(NA_real_)
        s <- f$snr_gamma[variables, cols, drop = FALSE]
        s <- s[rowSums(is.finite(s)) > 0L, , drop = FALSE]
        if (!nrow(s)) return(NA_real_)
        max(s, na.rm = TRUE)
    }

    .log_info("calibrate_gamma_threshold: %d refits (panel=%d, cols=%d)",
              n_boot, length(variables), length(cols))
    res <- .biomimic_parallel(seq_len(n_boot), one_boot,
                              export = c("X", "B", "Lam", "psi", "n", "p",
                                         "K", "use_ard", "anchors", "vnames",
                                         "variables", "cols", "max_iter",
                                         "seed"),
                              envir = environment(), n_cores = n_cores)
    boot_max <- unlist(res)
    boot_max <- boot_max[is.finite(boot_max)]
    # Require at least half the requested refits, and never fewer than 20
    # unless fewer than 20 were requested in the first place.
    if (length(boot_max) < max(ceiling(n_boot / 2), min(n_boot, 20L))) {
        stop(sprintf(
            "Only %d of %d bootstrap refits succeeded; not enough to calibrate.",
            length(boot_max), n_boot), call. = FALSE)
    }

    threshold <- as.numeric(stats::quantile(boot_max, level))
    .log_info("calibrate_gamma_threshold: threshold=%.2f (level=%.2f, B=%d)",
              threshold, level, length(boot_max))

    structure(
        list(threshold = threshold,
             level     = level,
             boot_max  = boot_max,
             n_boot    = length(boot_max),
             variables = variables,
             columns   = colnames(fit$X)[cols]),
        class = "biomimic_gamma_calibration"
    )
}


#' Print method for biomimic_gamma_calibration objects
#' @param x A \code{biomimic_gamma_calibration} object.
#' @param ... Ignored.
#' @export
print.biomimic_gamma_calibration <- function(x, ...) {
    cat("Direct-effect SNR threshold (parametric bootstrap)\n")
    cat(sprintf("  Threshold: %.2f (FWER %.0f%% over %d panel variables)\n",
                x$threshold, 100 * (1 - x$level), length(x$variables)))
    cat(sprintf("  Bootstrap refits: %d | covariate columns: %s\n",
                x$n_boot, paste(x$columns, collapse = ", ")))
    cat(sprintf("  Null maxima: q50=%.2f q95=%.2f max=%.2f\n",
                stats::median(x$boot_max),
                stats::quantile(x$boot_max, 0.95), max(x$boot_max)))
    invisible(x)
}
