# --------------------------------------------------------------------------
# Utility functions for biomimic
# --------------------------------------------------------------------------

# Null-coalescing operator (internal, no Rd page)
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Compute SNR from posterior mean and standard deviation
#'
#' @param mu Posterior mean (scalar or vector).
#' @param sigma Posterior standard deviation (scalar or vector).
#' @return Absolute SNR = |mu| / sigma.
#' @noRd
.compute_snr <- function(mu, sigma) {
    abs(mu) / pmax(sigma, .Machine$double.eps)
}

#' KL divergence between two Gamma distributions (rate parameterisation)
#'
#' KL( Gamma(aq, bq) || Gamma(ap, bp) ) with the second parameter the RATE
#' (density proportional to x^{a-1} exp(-b x)).
#' @noRd
.kl_gamma_rate <- function(aq, bq, ap, bp) {
    (aq - ap) * digamma(aq) - lgamma(aq) + lgamma(ap) +
        ap * (log(bq) - log(bp)) + aq * (bp - bq) / bq
}

#' Woodbury matrix identity for K=1
#'
#' Computes the inverse of (I_K + A' D^(-1) A) efficiently when K=1
#' (reduces to scalar division).
#'
#' @param lambda_vec Numeric vector of loadings (length p).
#' @param psi_inv_vec Numeric vector of precision values (length p).
#' @return Scalar: the (1,1) element of the inverse.
#' @noRd
.woodbury_k1 <- function(lambda_vec, psi_inv_vec) {
    1.0 / (1.0 + sum(lambda_vec^2 * psi_inv_vec))
}

#' Woodbury matrix identity for general K
#'
#' Computes the inverse of (I_K + Lambda' diag(psi_inv) Lambda).
#'
#' @param Lambda Numeric matrix (p x K) of loadings.
#' @param psi_inv Numeric vector (length p) of precision values.
#' @return K x K matrix.
#' @noRd
.woodbury <- function(Lambda, psi_inv) {
    # M = I_K + Lambda' diag(psi_inv) Lambda
    K <- ncol(Lambda)
    M <- diag(K) + crossprod(Lambda * psi_inv, Lambda)
    solve(M)
}

#' Print method for vb_ard objects
#' @param x A `vb_ard` object.
#' @param ... Ignored.
#' @export
print.vb_ard <- function(x, ...) {
    cat("VB-ARD fit\n")
    cat(sprintf("  Dimensions: n=%d, p=%d, q=%d, K=%d\n",
                x$n, x$p, x$q, x$K))
    cat(sprintf("  Converged: %s (iter=%d)\n",
                x$converged, x$n_iter))
    cat(sprintf("  Final ELBO: %.2f\n", x$elbo[x$n_iter]))
    if (!is.null(x$snr_gamma)) {
        cat("  Direct effects (Gamma): estimated\n")
    }
    invisible(x)
}

#' Print method for biomimic_selection objects
#' @param x A `biomimic_selection` object.
#' @param ... Ignored.
#' @importFrom utils head
#' @export
print.biomimic_selection <- function(x, ...) {
    cat(sprintf("biomimic selection: %d variables selected\n",
                length(x$variables)))
    cat("  Variables:", paste(utils::head(x$variables, 10), collapse = ", "))
    if (length(x$variables) > 10) cat(", ...")
    cat("\n")
    if (length(x$direct_effects) > 0) {
        cat(sprintf("  Direct effects flagged: %d variables\n",
                    length(x$direct_effects)))
    }
    if (x$K > 1) {
        cat(sprintf("  Factor structure: K=%d, groups: %s\n",
                    x$K,
                    paste(vapply(x$indicator_groups, length, integer(1)),
                          collapse = ", ")))
    }
    invisible(x)
}

