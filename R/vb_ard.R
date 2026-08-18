# --------------------------------------------------------------------------
# VB-ARD: Variational Bayes with Automatic Relevance Determination
# Stage 1 of the biomimic pipeline
# --------------------------------------------------------------------------

#' Fit a VB-ARD Factor Model for Variable Screening
#'
#' Fits a Bayesian factor model with Automatic Relevance Determination (ARD)
#' priors using Coordinate Ascent Variational Inference (CAVI). The model
#' includes optional direct-effect paths (Gamma) from covariates to observed
#' variables.
#'
#' @section Model:
#' The generative model is:
#' \deqn{y_i = \Lambda h_i + \Gamma x_i + \varepsilon_i, \quad
#'       \varepsilon_i \sim N(0, \Psi)}
#' \deqn{h_i = B' x_i + \eta_i, \quad \eta_i \sim N(0, I_K)}
#'
#' where \eqn{\Lambda} has ARD priors \eqn{\lambda_{jk} \sim N(0, \tau_{jk}^{-1})},
#' \eqn{\tau_{jk} \sim \text{Gamma}(a_\tau, b_\tau)}, and \eqn{\Gamma} has ARD priors
#' \eqn{\gamma_{jl} \sim N(0, \xi_{jl}^{-1})},
#' \eqn{\xi_{jl} \sim \text{Gamma}(a_\xi, b_\xi)}.
#'
#' @param Y Numeric matrix (n x p) of observed variables.
#' @param X Numeric matrix (n x q) of covariates (should include intercept).
#' @param K Integer, number of latent factors (default 1).
#' @param use_ard Logical, whether to use ARD priors on loadings
#'   (default TRUE). Recommended TRUE for every K: with ARD off, the
#'   loading precisions stay fixed at a vague constant and the
#'   direct-effect screening statistic degrades sharply for K >= 2
#'   (in simulation, its null false-positive rate at the default
#'   threshold rises several-fold).
#' @param use_gamma Logical, whether to estimate direct effects Gamma
#'   (default TRUE).
#' @param max_iter Maximum number of CAVI iterations (default 500).
#' @param tol Convergence tolerance for relative ELBO change (default 1e-6).
#' @param a_tau,b_tau Shape and rate hyperparameters for ARD prior on Lambda
#'   (default 1e-3 each).
#' @param a_psi,b_psi Shape and rate hyperparameters for prior on inverse psi
#'   (default 1e-3 each).
#' @param a_xi,b_xi Shape and rate hyperparameters for ARD prior on Gamma
#'   (default 1e-3 each).
#' @param gamma_anchor Optional character vector of variable names (or integer
#'   column indices of \code{Y}) whose direct effects are held at exactly zero.
#'   These \emph{anchor} indicators identify the split between the structural
#'   coefficients \code{B} and the direct effects \code{Gamma}. Without them
#'   the split is not identified by the likelihood (see Details), and the
#'   separation rests on the ARD prior alone. Supply at least \code{K} anchors
#'   with linearly independent loading vectors; \code{\link{select_gamma_anchor}}
#'   picks them from a first-pass fit. Default \code{NULL} reproduces the
#'   unanchored (prior-regularised) behaviour.
#' @param verbose Logical, print progress (default FALSE).
#'
#' @details
#' \strong{Identification of B versus Gamma.} Marginalising the latent scores
#' gives \eqn{y_i \mid x_i \sim N((\Lambda B' + \Gamma) x_i,\ \Lambda\Lambda' +
#' \Psi)}, so the likelihood depends on \code{B} and \code{Gamma} only through
#' the sum \eqn{\Lambda B' + \Gamma}. For any \eqn{\Delta} (q x K) the pair
#' \eqn{(B + \Delta,\ \Gamma - \Lambda\Delta')} has an identical likelihood:
#' a flat direction of dimension \eqn{qK}. \code{Lambda} itself is identified
#' (up to rotation) from the conditional covariance and is unaffected.
#'
#' Fixing \eqn{\gamma_j = 0} for a set \eqn{A} of anchors makes
#' \eqn{\Pi_A = \Lambda_A B'} and pins \code{B} uniquely whenever
#' \eqn{\Lambda_A} has rank \code{K} -- hence the requirement of at least
#' \code{K} anchors with linearly independent loadings.
#'
#' @return A list of class \code{"vb_ard"} with components:
#' \describe{
#'   \item{mu_lambda}{Posterior mean of Lambda (p x K)}
#'   \item{sigma_lambda}{Posterior covariance of Lambda (list of p matrices, each K x K)}
#'   \item{mu_gamma}{Posterior mean of Gamma (p x q), or NULL if use_gamma=FALSE}
#'   \item{sigma_gamma}{Posterior covariance of Gamma (list of p matrices), or NULL}
#'   \item{mu_B}{Posterior mean of B (q x K)}
#'   \item{mu_eta}{Posterior mean of latent factors (n x K)}
#'   \item{sigma_eta}{Posterior covariance of latent factors (K x K, shared)}
#'   \item{snr_lambda}{Signal-to-noise ratio for each loading (p x K)}
#'   \item{snr_gamma}{Signal-to-noise ratio for direct effects (p x q), or NULL}
#'   \item{elbo}{Vector of ELBO values per iteration}
#'   \item{converged}{Logical, whether CAVI converged}
#'   \item{n_iter}{Number of iterations run}
#'   \item{K}{Number of latent factors}
#'   \item{p}{Number of observed variables}
#'   \item{n}{Number of observations}
#'   \item{q}{Number of covariates}
#'   \item{var_names}{Variable names from colnames(Y)}
#'   \item{Y}{Original data matrix}
#'   \item{X}{Original covariate matrix}
#' }
#'
#' @examples
#' set.seed(42)
#' n <- 200; p <- 50; q <- 2
#' X <- cbind(1, rbinom(n, 1, 0.5))
#' B <- c(0, 1)
#' Lambda <- c(rep(1, 15), rep(0, 35))
#' h <- X %*% B + rnorm(n)
#' Y <- h %*% t(Lambda) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
#'
#' fit <- vb_ard(Y, X, K = 1)
#' head(rank_variables(fit))
#'
#' @export
vb_ard <- function(Y, X, K = 1L,
                   use_ard = TRUE,
                   use_gamma = TRUE,
                   max_iter = 500L,
                   tol = 1e-6,
                   a_tau = 1e-3, b_tau = 1e-3,
                   a_psi = 1e-3, b_psi = 1e-3,
                   a_xi = 1e-3, b_xi = 1e-3,
                   gamma_anchor = NULL,
                   verbose = FALSE) {

    # --- Input validation ---
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    stopifnot(is.numeric(Y), is.numeric(X))
    stopifnot(nrow(Y) == nrow(X))
    stopifnot(K >= 1L)

    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(X)
    var_names <- colnames(Y) %||% paste0("V", seq_len(p))

    # --- Resolve the identification anchors ---
    # Without at least K indicators whose direct effects are held at exactly
    # zero, B and Gamma are not identified by the likelihood: for any
    # Delta (q x K), (B + Delta, Gamma - Lambda Delta') gives the same
    # reduced-form mean (Lambda B' + Gamma) X' and the same covariance.
    # Anchoring K rows of Gamma at zero removes that flat direction.
    anchor_idx <- integer(0)
    if (!is.null(gamma_anchor)) {
        if (!use_gamma) {
            stop("gamma_anchor requires use_gamma = TRUE.", call. = FALSE)
        }
        anchor_idx <- if (is.character(gamma_anchor)) {
            miss <- setdiff(gamma_anchor, var_names)
            if (length(miss)) {
                stop("gamma_anchor not found in colnames(Y): ",
                     paste(miss, collapse = ", "), call. = FALSE)
            }
            match(gamma_anchor, var_names)
        } else {
            idx <- as.integer(gamma_anchor)
            if (any(idx < 1L | idx > p)) {
                stop("gamma_anchor indices out of range [1, ", p, "].",
                     call. = FALSE)
            }
            idx
        }
        anchor_idx <- sort(unique(anchor_idx))
        if (length(anchor_idx) < K) {
            warning(sprintf(
                paste0("Only %d anchor(s) supplied for K = %d factors. ",
                       "Identification of B requires at least K anchors ",
                       "whose loading vectors are linearly independent."),
                length(anchor_idx), K), call. = FALSE)
        }
        if (length(anchor_idx) >= p) {
            stop("gamma_anchor cannot cover every indicator.", call. = FALSE)
        }
    }

    # Center the indicators. The observation model has no per-indicator
    # intercept (y_ij = lambda_j' h_i + gamma_j' x_i + eps), which is exactly
    # valid only for mean-zero indicators; centering also keeps the indicator
    # means out of the ARD-penalised direct effects. Idempotent when Y is
    # already centered upstream. Column means are stored for reconstruction.
    Y_mean <- colMeans(Y)
    Y <- sweep(Y, 2, Y_mean, "-")

    .log_info("vb_ard: n=%d, p=%d, q=%d, K=%d, use_ard=%s, use_gamma=%s",
              n, p, q, K, use_ard, use_gamma)

    # --- Precompute ---
    XtX <- crossprod(X)  # q x q

    # --- Initialisation ---
    .log_debug("vb_ard: initialising CAVI (SVD-based)")
    init <- .cavi_init(Y, X, n, p, q, K, use_ard, use_gamma,
                       a_tau, b_tau, a_psi, b_psi, a_xi, b_xi)

    mu_lambda    <- init$mu_lambda      # p x K
    sigma_lambda <- init$sigma_lambda   # list of p (K x K)
    mu_B         <- init$mu_B           # q x K
    sigma_B      <- NULL                # q x q row covariance of q(B)
    sigma_B_sq   <- 100                 # prior variance of B rows
    mu_eta       <- init$mu_eta         # n x K
    sigma_eta    <- init$sigma_eta      # K x K
    E_psi_inv    <- init$E_psi_inv      # p
    E_tau        <- init$E_tau          # p x K
    mu_gamma     <- init$mu_gamma       # p x q or NULL
    sigma_gamma  <- init$sigma_gamma    # list of p (q x q) or NULL
    E_xi         <- init$E_xi           # p x q or NULL

    elbo_vec <- numeric(max_iter)
    converged <- FALSE

    for (it in seq_len(max_iter)) {

        # --- Update q(h) ---
        res_eta <- .cavi_update_eta(Y, X, mu_lambda, sigma_lambda, mu_B,
                                    mu_gamma, E_psi_inv, n, p, K, q)
        mu_eta    <- res_eta$mu_eta
        sigma_eta <- res_eta$sigma_eta

        # --- Update q(Lambda_j) ---
        res_lam <- .cavi_update_lambda(Y, X, mu_eta, sigma_eta, mu_gamma,
                                       E_psi_inv, E_tau, n, p, K)
        mu_lambda    <- res_lam$mu_lambda
        sigma_lambda <- res_lam$sigma_lambda

        # --- Update q(Gamma_j) if requested ---
        if (use_gamma) {
            res_gam <- .cavi_update_gamma(Y, X, mu_eta, mu_lambda, E_psi_inv,
                                          E_xi, n, p, q, XtX, anchor_idx)
            mu_gamma    <- res_gam$mu_gamma
            sigma_gamma <- res_gam$sigma_gamma
        }

        # --- Update q(B) ---
        res_B   <- .cavi_update_B(X, mu_eta, XtX, q, K, sigma_B_sq = sigma_B_sq)
        mu_B    <- res_B$mu_B
        sigma_B <- res_B$sigma_B

        # --- Update q(tau) (ARD) ---
        if (use_ard) {
            E_tau <- .cavi_update_tau(mu_lambda, sigma_lambda, a_tau, b_tau,
                                     p, K)
        }

        # --- Update q(xi) (ARD on Gamma) ---
        if (use_gamma) {
            E_xi <- .cavi_update_xi(mu_gamma, sigma_gamma, a_xi, b_xi, p, q,
                                    anchor_idx)
        }

        # --- Update q(psi_j^{-1}) ---
        E_psi_inv <- .cavi_update_psi(Y, X, mu_eta, sigma_eta, mu_lambda,
                                      sigma_lambda, mu_gamma, sigma_gamma,
                                      mu_B, a_psi, b_psi, n, p, K, q)

        # --- ELBO ---
        elbo_vec[it] <- .compute_elbo(Y, X, mu_eta, sigma_eta, mu_lambda,
                                      sigma_lambda, mu_gamma, sigma_gamma,
                                      mu_B, sigma_B, E_psi_inv, E_tau, E_xi,
                                      a_tau, b_tau, a_psi, b_psi,
                                      a_xi, b_xi, sigma_B_sq,
                                      n, p, q, K, use_ard, use_gamma,
                                      anchor_idx)

        if (verbose && (it %% 50 == 0 || it == 1)) {
            message(sprintf("Iter %d: ELBO = %.2f", it, elbo_vec[it]))
        }

        # Convergence check (relative change)
        if (it > 1) {
            rel_change <- abs((elbo_vec[it] - elbo_vec[it - 1]) /
                              (abs(elbo_vec[it - 1]) + 1e-10))
            if (rel_change < tol) {
                converged <- TRUE
                if (verbose) message(sprintf("Converged at iteration %d", it))
                break
            }
        }
    }

    n_iter <- it
    elbo_vec <- elbo_vec[seq_len(n_iter)]

    if (converged) {
        .log_info("vb_ard: converged at iteration %d, ELBO=%.2f", n_iter,
                  elbo_vec[n_iter])
    } else {
        .log_warn("vb_ard: did NOT converge after %d iterations", n_iter)
    }

    # --- Identification: sign (K=1) and varimax rotation (K>=2) ---
    # For K>=2 the loadings are identified only up to an orthogonal rotation
    # (Lambda h = (Lambda Q)(Q' h)). We resolve it with varimax and rotate the
    # WHOLE posterior consistently -- means AND covariances -- so the loadings,
    # factor scores, B and the SNR all live in the same (rotated) frame; the
    # rotated loadings are then used for ranking and for Stage-2 factor grouping.
    rotated_lambda  <- mu_lambda
    rotation_matrix <- NULL
    factor_groups   <- NULL

    if (K >= 2L) {
        .log_info("vb_ard: applying varimax rotation (K=%d)", K)
        rot <- varimax(mu_lambda, normalize = TRUE)
        rotation_matrix <- rot$rotmat                 # orthogonal Q (K x K)
        rotated_lambda <- mu_lambda %*% rotation_matrix
        mu_eta <- mu_eta %*% rotation_matrix
        mu_B   <- mu_B %*% rotation_matrix
        sigma_lambda <- lapply(sigma_lambda, function(S)
            crossprod(rotation_matrix, S %*% rotation_matrix))   # Q' Sigma_j Q

        factor_groups <- apply(abs(rotated_lambda), 1, which.max)

        # Post-rotation sign identification (dominant loading per factor > 0);
        # sign flip S acts on loadings/scores/B and on the covariances as S Sigma S
        # (the diagonal, hence the SDs, is preserved).
        sign_vec <- rep(1, K)
        for (k in seq_len(K)) {
            dom <- which.max(abs(rotated_lambda[, k]))
            if (rotated_lambda[dom, k] < 0) sign_vec[k] <- -1
        }
        if (any(sign_vec < 0)) {
            S_sign <- diag(sign_vec, K)
            rotated_lambda <- rotated_lambda %*% S_sign
            mu_eta <- mu_eta %*% S_sign
            mu_B   <- mu_B %*% S_sign
            sign_outer <- tcrossprod(sign_vec)
            sigma_lambda <- lapply(sigma_lambda, function(S) S * sign_outer)
        }
        mu_lambda <- rotated_lambda   # adopt rotated frame as primary
    } else {
        dominant_j <- which.max(abs(mu_lambda[, 1]))
        if (mu_lambda[dominant_j, 1] < 0) {
            mu_lambda[, 1] <- -mu_lambda[, 1]
            mu_eta[, 1]    <- -mu_eta[, 1]
            mu_B[, 1]      <- -mu_B[, 1]
        }
        rotated_lambda <- mu_lambda
    }

    # --- Compute SNRs (means and SDs now in the same coordinate frame) ---
    sd_lambda <- .extract_posterior_sd(sigma_lambda, p, K)
    snr_lambda <- .compute_snr(mu_lambda, sd_lambda)
    dimnames(snr_lambda) <- list(var_names, paste0("K", seq_len(K)))

    snr_gamma <- NULL
    if (use_gamma) {
        sd_gamma <- .extract_posterior_sd_rect(sigma_gamma, p, q)
        snr_gamma <- .compute_snr(mu_gamma, sd_gamma)
        dimnames(snr_gamma) <- list(var_names, colnames(X))
        # Anchored rows have mu = 0 and sd = 0: the SNR is undefined, not zero.
        # NA keeps them out of the direct-effect ranking instead of ranking
        # them last on a fabricated value.
        if (length(anchor_idx)) snr_gamma[anchor_idx, ] <- NA_real_
    }

    structure(
        list(
            mu_lambda    = mu_lambda,
            sigma_lambda = sigma_lambda,
            rotated_lambda  = rotated_lambda,
            rotation_matrix = rotation_matrix,
            factor_groups   = factor_groups,
            mu_gamma     = mu_gamma,
            sigma_gamma  = sigma_gamma,
            mu_B         = mu_B,
            sigma_B      = sigma_B,
            mu_eta       = mu_eta,
            sigma_eta    = sigma_eta,
            snr_lambda   = snr_lambda,
            snr_gamma    = snr_gamma,
            elbo         = elbo_vec,
            converged    = converged,
            n_iter       = n_iter,
            K            = K,
            p            = p,
            n            = n,
            q            = q,
            var_names    = var_names,
            Y            = Y,
            Y_mean       = Y_mean,
            X            = X,
            use_ard      = use_ard,
            use_gamma    = use_gamma,
            gamma_anchor = if (length(anchor_idx)) var_names[anchor_idx]
                           else NULL,
            anchor_idx   = anchor_idx,
            E_psi_inv    = E_psi_inv,
            E_tau        = E_tau,
            E_xi         = E_xi,
            sigma_B_sq   = sigma_B_sq
        ),
        class = "vb_ard"
    )
}


# ===========================================================================
# CAVI Initialisation
# ===========================================================================

#' @noRd
.cavi_init <- function(Y, X, n, p, q, K, use_ard, use_gamma,
                       a_tau, b_tau, a_psi, b_psi, a_xi, b_xi) {

    # Lambda: SVD-based initialisation
    Y_resid <- Y - X %*% solve(crossprod(X), crossprod(X, Y))  # residuals
    svd_out <- svd(Y_resid, nu = K, nv = K)
    mu_lambda <- svd_out$v[, seq_len(K), drop = FALSE] *
                 rep(svd_out$d[seq_len(K)], each = p) / sqrt(n)
    dim(mu_lambda) <- c(p, K)

    sigma_lambda <- replicate(p, diag(0.01, K), simplify = FALSE)

    # B: OLS of h on X (initialise h from SVD scores)
    mu_eta <- svd_out$u[, seq_len(K), drop = FALSE] * sqrt(n)
    dim(mu_eta) <- c(n, K)
    mu_B <- solve(crossprod(X), crossprod(X, mu_eta))  # q x K

    sigma_eta <- diag(1, K)

    # Psi: residual variances
    Y_hat <- mu_eta %*% t(mu_lambda)
    resid <- Y - Y_hat
    psi_init <- pmax(colMeans(resid^2), 0.01)
    # E[psi_j^{-1}] initialised as 1/var
    E_psi_inv <- 1.0 / psi_init

    # ARD precisions
    E_tau <- matrix(a_tau / b_tau, nrow = p, ncol = K)
    if (!use_ard) {
        E_tau <- matrix(1e-3, nrow = p, ncol = K)
    }

    # Gamma (direct effects)
    mu_gamma <- NULL
    sigma_gamma <- NULL
    E_xi <- NULL
    if (use_gamma) {
        mu_gamma <- matrix(0, nrow = p, ncol = q)
        sigma_gamma <- replicate(p, diag(0.01, q), simplify = FALSE)
        E_xi <- matrix(a_xi / b_xi, nrow = p, ncol = q)
    }

    list(mu_lambda = mu_lambda, sigma_lambda = sigma_lambda,
         mu_B = mu_B, mu_eta = mu_eta, sigma_eta = sigma_eta,
         E_psi_inv = E_psi_inv, E_tau = E_tau,
         mu_gamma = mu_gamma, sigma_gamma = sigma_gamma, E_xi = E_xi)
}


# ===========================================================================
# CAVI Updates
# ===========================================================================

#' @noRd
.cavi_update_eta <- function(Y, X, mu_lambda, sigma_lambda, mu_B,
                             mu_gamma, E_psi_inv, n, p, K, q) {
    # E[Lambda' Psi^{-1} Lambda] = sum_j psi_j^{-1} (mu_j mu_j' + Sigma_j)
    E_LtPsiL <- matrix(0, K, K)
    for (j in seq_len(p)) {
        mu_j <- mu_lambda[j, , drop = FALSE]  # 1 x K
        E_LtPsiL <- E_LtPsiL + E_psi_inv[j] *
            (crossprod(mu_j) + sigma_lambda[[j]])
    }

    # Sigma_eta = (I_K + E[Lambda' Psi^{-1} Lambda])^{-1}
    sigma_eta <- solve(diag(K) + E_LtPsiL)

    # Residual after removing Gamma contribution
    if (!is.null(mu_gamma)) {
        Y_adj <- Y - X %*% t(mu_gamma)  # n x p
    } else {
        Y_adj <- Y
    }

    # mu_eta_i = Sigma_eta (E[Lambda]' E[Psi^{-1}] y_adj_i + E[B]' x_i)
    # Vectorised: mu_eta = (Y_adj diag(psi_inv) Lambda + X B) Sigma_eta'
    # since Sigma_eta is symmetric: Sigma_eta' = Sigma_eta
    mu_eta <- (Y_adj %*% (E_psi_inv * mu_lambda) + X %*% mu_B) %*% sigma_eta

    list(mu_eta = mu_eta, sigma_eta = sigma_eta)
}


#' @noRd
.cavi_update_lambda <- function(Y, X, mu_eta, sigma_eta, mu_gamma,
                                E_psi_inv, E_tau, n, p, K) {
    # E[h' h] = sum_i (mu_i mu_i' + Sigma_eta) = H'H + n * Sigma_eta
    E_HtH <- crossprod(mu_eta) + n * sigma_eta  # K x K

    # Residual: Y_adj = Y - Gamma X (if Gamma present)
    if (!is.null(mu_gamma)) {
        Y_adj <- Y - X %*% t(mu_gamma)
    } else {
        Y_adj <- Y
    }

    # E[h' Y_adj] -> K x p
    HtY <- crossprod(mu_eta, Y_adj)  # K x n . n x p = K x p

    mu_lambda <- matrix(0, p, K)
    sigma_lambda <- vector("list", p)

    for (j in seq_len(p)) {
        # Sigma_lambda_j = (psi_j^{-1} E[h'h] + diag(tau_j))^{-1}
        # use_ard enters through E_tau: it is re-estimated each sweep when ARD
        # is on and held at its fixed initial value when it is off, so the
        # update itself is the same expression either way.
        prec <- E_psi_inv[j] * E_HtH + diag(E_tau[j, ], nrow = K)
        sigma_lambda[[j]] <- solve(prec)

        # mu_lambda_j = psi_j^{-1} Sigma_lambda_j E[h]' y_adj_j
        mu_lambda[j, ] <- E_psi_inv[j] * sigma_lambda[[j]] %*% HtY[, j]
    }

    list(mu_lambda = mu_lambda, sigma_lambda = sigma_lambda)
}


#' @noRd
.cavi_update_gamma <- function(Y, X, mu_eta, mu_lambda, E_psi_inv,
                               E_xi, n, p, q, XtX, anchor_idx = integer(0)) {
    # Residual: Y - H Lambda'
    Y_adj <- Y - mu_eta %*% t(mu_lambda)  # n x p

    # XtY_adj = X' Y_adj -> q x p
    XtY <- crossprod(X, Y_adj)

    mu_gamma <- matrix(0, p, q)
    sigma_gamma <- vector("list", p)

    zero_q <- matrix(0, q, q)
    is_anchor <- logical(p)
    is_anchor[anchor_idx] <- TRUE

    for (j in seq_len(p)) {
        if (is_anchor[j]) {
            # Anchored indicator: gamma_j is fixed at exactly zero, so q(gamma_j)
            # is a point mass at 0 (no update, no posterior spread).
            sigma_gamma[[j]] <- zero_q
            next
        }
        # Sigma_gamma_j = (psi_j^{-1} X'X + diag(xi_j))^{-1}
        prec <- E_psi_inv[j] * XtX + diag(E_xi[j, ], nrow = q)
        sigma_gamma[[j]] <- solve(prec)

        # mu_gamma_j = psi_j^{-1} Sigma_gamma_j X' y_adj_j
        mu_gamma[j, ] <- E_psi_inv[j] * sigma_gamma[[j]] %*% XtY[, j]
    }

    list(mu_gamma = mu_gamma, sigma_gamma = sigma_gamma)
}


#' @noRd
.cavi_update_B <- function(X, mu_eta, XtX, q, K, sigma_B_sq = 100) {
    # q(B) = Matrix-Normal(mu_B, Sigma_{B,row}, I_K) with shared row covariance
    # Sigma_{B,row} = (X'X + sigma_B^{-2} I_q)^{-1} and mean mu_B = Sigma_B X' E[H].
    sigma_B <- solve(XtX + (1.0 / sigma_B_sq) * diag(q))
    mu_B    <- sigma_B %*% crossprod(X, mu_eta)
    list(mu_B = mu_B, sigma_B = sigma_B)
}


#' @noRd
.cavi_update_tau <- function(mu_lambda, sigma_lambda, a_tau, b_tau, p, K) {
    # q(tau_jk) = Gamma(a_hat, b_hat)
    # a_hat = a_tau + 0.5
    # b_hat = b_tau + 0.5 * (mu_jk^2 + Sigma_j_kk)
    # E[tau_jk] = a_hat / b_hat
    a_hat <- a_tau + 0.5
    E_tau <- matrix(0, p, K)
    for (j in seq_len(p)) {
        for (k in seq_len(K)) {
            b_hat <- b_tau + 0.5 * (mu_lambda[j, k]^2 +
                                    sigma_lambda[[j]][k, k])
            E_tau[j, k] <- a_hat / b_hat
        }
    }
    E_tau
}


#' @noRd
.cavi_update_xi <- function(mu_gamma, sigma_gamma, a_xi, b_xi, p, q,
                            anchor_idx = integer(0)) {
    a_hat <- a_xi + 0.5
    E_xi <- matrix(0, p, q)
    is_anchor <- logical(p)
    is_anchor[anchor_idx] <- TRUE
    for (j in seq_len(p)) {
        if (is_anchor[j]) {
            # gamma_j is fixed at zero, so its ARD precision carries no
            # information; hold it at the prior mean rather than letting
            # E[gamma^2] = 0 drive b_hat to b_xi and the precision to a/b_xi.
            E_xi[j, ] <- a_xi / b_xi
            next
        }
        for (l in seq_len(q)) {
            b_hat <- b_xi + 0.5 * (mu_gamma[j, l]^2 +
                                   sigma_gamma[[j]][l, l])
            E_xi[j, l] <- a_hat / b_hat
        }
    }
    E_xi
}


#' @noRd
.cavi_update_psi <- function(Y, X, mu_eta, sigma_eta, mu_lambda,
                             sigma_lambda, mu_gamma, sigma_gamma,
                             mu_B, a_psi, b_psi, n, p, K, q) {
    # q(psi_j^{-1}) = Gamma(a_hat, b_hat)
    # a_hat = a_psi + n/2
    # b_hat = b_psi + 0.5 * sum_i E[(y_ij - lambda_j' h_i - gamma_j' x_i)^2]
    a_hat <- a_psi + n / 2.0

    # Precompute E[h_i h_i'] contribution
    E_HtH <- crossprod(mu_eta) + n * sigma_eta  # K x K

    # Predicted values
    Y_hat <- mu_eta %*% t(mu_lambda)
    if (!is.null(mu_gamma)) {
        Y_hat <- Y_hat + X %*% t(mu_gamma)
    }
    resid <- Y - Y_hat  # n x p

    E_psi_inv <- numeric(p)
    XtX_full <- if (!is.null(mu_gamma)) crossprod(X) else NULL
    for (j in seq_len(p)) {
        mlj <- mu_lambda[j, ]
        # E[(y_j - lambda_j' h - gamma_j' x)^2] summed over i =
        #   sum_i resid_ij^2                          [deterministic residual]
        # + tr(Sigma_lambda_j E[h'h])                 [loading uncertainty]
        # + n * mu_lambda_j' Sigma_eta mu_lambda_j    [factor uncertainty via mean]
        # + tr(Sigma_gamma_j X'X)                     [direct-effect uncertainty]
        # where tr(Sigma_lambda_j E[h'h]) already carries n*tr(Sigma_eta Sigma_lambda_j).
        b_hat <- b_psi + 0.5 * (
            sum(resid[, j]^2) +
            sum(sigma_lambda[[j]] * E_HtH) +
            n * as.numeric(crossprod(mlj, sigma_eta %*% mlj))
        )
        if (!is.null(sigma_gamma)) {
            # full trace tr(Sigma_gamma_j X'X), not just its diagonal
            b_hat <- b_hat + 0.5 * sum(sigma_gamma[[j]] * XtX_full)
        }
        E_psi_inv[j] <- a_hat / b_hat
    }

    E_psi_inv
}


# ===========================================================================
# ELBO Computation
# ===========================================================================

#' Evidence Lower Bound (full, exact for the conjugate MIMIC model)
#'
#' Computes the complete ELBO L = E_q[log p(Y, h, Lambda, Gamma, B, psi, tau)]
#' - E_q[log q]. Uses E[log psi^{-1}] = digamma(a_tilde) - log(b_tilde) (NOT
#' log E[psi^{-1}]); the exact residual second moment (with the
#' n mu_lambda' Sigma_eta mu_lambda term and the full tr(Sigma_gamma X'X)); the
#' q(Lambda) and q(Gamma) entropies; the -KL of q(psi), q(tau_lambda),
#' q(tau_gamma) and the Matrix-Normal q(B). For the conjugate updates here this
#' increases monotonically.
#' @noRd
.compute_elbo <- function(Y, X, mu_eta, sigma_eta, mu_lambda, sigma_lambda,
                          mu_gamma, sigma_gamma, mu_B, sigma_B, E_psi_inv,
                          E_tau, E_xi,
                          a_tau, b_tau, a_psi, b_psi, a_xi, b_xi, sigma_B_sq,
                          n, p, q, K, use_ard, use_gamma,
                          anchor_idx = integer(0)) {

    a_hat_psi <- a_psi + n / 2.0
    a_hat_tau <- a_tau + 0.5
    a_hat_xi  <- a_xi + 0.5
    log2pi    <- log(2 * pi)
    E_HtH     <- crossprod(mu_eta) + n * sigma_eta
    XtX       <- if (use_gamma && q > 0L) crossprod(X) else NULL

    Y_hat <- mu_eta %*% t(mu_lambda)
    if (!is.null(mu_gamma)) Y_hat <- Y_hat + X %*% t(mu_gamma)
    resid <- Y - Y_hat

    b_hat_psi     <- a_hat_psi / E_psi_inv          # length p
    E_log_psi_inv <- digamma(a_hat_psi) - log(b_hat_psi)

    L <- 0
    for (j in seq_len(p)) {
        mlj <- mu_lambda[j, ]
        sq_resid <- sum(resid[, j]^2) +
                    sum(sigma_lambda[[j]] * E_HtH) +
                    n * as.numeric(crossprod(mlj, sigma_eta %*% mlj))
        if (use_gamma && !is.null(sigma_gamma)) {
            sq_resid <- sq_resid + sum(sigma_gamma[[j]] * XtX)
        }
        # E[log p(y_.j | .)]  with  E[log psi^-1] = digamma(a~) - log(b~)
        L <- L + 0.5 * n * (E_log_psi_inv[j] - log2pi) -
             0.5 * E_psi_inv[j] * sq_resid
        # -KL(q(psi_j) || p(psi_j))
        L <- L - .kl_gamma_rate(a_hat_psi, b_hat_psi[j], a_psi, b_psi)
        # q(lambda_j) differential entropy: 0.5 log det(2 pi e Sigma_lambda_j)
        L <- L + 0.5 * (K * (1 + log2pi) +
                        log(det(sigma_lambda[[j]]) + 1e-300))
        for (k in seq_len(K)) {
            E_lam_sq <- mu_lambda[j, k]^2 + sigma_lambda[[j]][k, k]
            if (use_ard) {
                b_hat_t <- a_hat_tau / E_tau[j, k]
                E_log_t <- digamma(a_hat_tau) - log(b_hat_t)
                L <- L + 0.5 * E_log_t - 0.5 * log2pi -
                     0.5 * E_tau[j, k] * E_lam_sq -
                     .kl_gamma_rate(a_hat_tau, b_hat_t, a_tau, b_tau)
            } else {
                L <- L + 0.5 * log(E_tau[j, k] + 1e-300) - 0.5 * log2pi -
                     0.5 * E_tau[j, k] * E_lam_sq
            }
        }
    }

    # --- Gamma direct-effect block (q(Gamma) entropy + Gamma/xi prior+KL) ---
    if (use_gamma && !is.null(mu_gamma)) {
        is_anchor <- logical(p)
        is_anchor[anchor_idx] <- TRUE
        for (j in seq_len(p)) {
            # Anchored rows are fixed at zero, not estimated: they contribute
            # no q(gamma_j) entropy and no gamma/xi prior or KL term.
            if (is_anchor[j]) next
            L <- L + 0.5 * (q * (1 + log2pi) +
                            log(det(sigma_gamma[[j]]) + 1e-300))
            for (l in seq_len(q)) {
                E_g_sq  <- mu_gamma[j, l]^2 + sigma_gamma[[j]][l, l]
                b_hat_x <- a_hat_xi / E_xi[j, l]
                E_log_x <- digamma(a_hat_xi) - log(b_hat_x)
                L <- L + 0.5 * E_log_x - 0.5 * log2pi -
                     0.5 * E_xi[j, l] * E_g_sq -
                     .kl_gamma_rate(a_hat_xi, b_hat_x, a_xi, b_xi)
            }
        }
    }

    # --- Latent factors: -KL(q(h) || p(h|B)) ---
    predicted_eta <- X %*% mu_B
    diff_eta <- mu_eta - predicted_eta
    # extra B-uncertainty contribution K * tr(Sigma_B X'X) from E[(B'x)'(B'x)]
    bxx <- if (!is.null(sigma_B)) K * sum(sigma_B * crossprod(X)) else 0
    kl_eta <- 0.5 * (sum(diff_eta^2) + bxx +
                     n * (sum(diag(sigma_eta)) - K -
                          log(det(sigma_eta) + 1e-300)))
    L <- L - kl_eta

    # --- B block: -KL(q(B) || p(B)) for the Matrix-Normal q(B) ---
    if (!is.null(mu_B)) {
        if (!is.null(sigma_B)) {
            q_dim <- nrow(mu_B)
            kl_B <- 0.5 * K * (sum(diag(sigma_B)) / sigma_B_sq - q_dim +
                               q_dim * log(sigma_B_sq) -
                               log(det(sigma_B) + 1e-300)) +
                    0.5 * sum(mu_B^2) / sigma_B_sq
            L <- L - kl_B
        } else {
            L <- L - 0.5 * sum(mu_B^2) / sigma_B_sq
        }
    }

    L
}


# ===========================================================================
# Helper: extract posterior standard deviations
# ===========================================================================

#' @noRd
.extract_posterior_sd <- function(sigma_list, p, K) {
    sd_mat <- matrix(0, p, K)
    for (j in seq_len(p)) {
        sd_mat[j, ] <- sqrt(pmax(diag(sigma_list[[j]]), 0))
    }
    sd_mat
}

#' @noRd
.extract_posterior_sd_rect <- function(sigma_list, p, d) {
    sd_mat <- matrix(0, p, d)
    for (j in seq_len(p)) {
        sd_mat[j, ] <- sqrt(pmax(diag(sigma_list[[j]]), 0))
    }
    sd_mat
}


# ===========================================================================
# Anchor selection
# ===========================================================================

#' Choose Identification Anchors for the Direct-Effect Model
#'
#' Picks the indicators whose direct effects should be held at zero so that the
#' structural coefficients \code{B} and the direct effects \code{Gamma} are
#' identified. See the Details section of \code{\link{vb_ard}} for why the
#' unanchored model is not identified by the likelihood.
#'
#' Anchors are the top indicators by loading signal-to-noise ratio. That
#' ranking is legitimate for this purpose: \code{Lambda} is identified from the
#' conditional covariance \eqn{\Lambda\Lambda' + \Psi}, which does not involve
#' the confounded \code{B}/\code{Gamma} split, so the choice does not use the
#' quantity it is meant to identify.
#'
#' @param fit A \code{vb_ard} object from a first-pass fit.
#' @param n_anchor Integer, how many anchors to return. Defaults to \code{K},
#'   the minimum that identifies \code{B}.
#' @param pool Integer, how many top-SNR indicators form the candidate pool
#'   from which anchors are drawn. Defaults to \code{max(3 * n_anchor, 10)},
#'   capped at the number of indicators.
#' @param tol Numeric, condition-number tolerance for the linear-independence
#'   check on the selected loading rows (default 1e-8).
#'
#' @return Character vector of variable names, suitable for the
#'   \code{gamma_anchor} argument of \code{\link{vb_ard}}.
#'
#' @details
#' For \code{K >= 2} the anchors' loading vectors must be linearly independent,
#' otherwise \eqn{\Lambda_A} is rank deficient and \code{B} stays unidentified.
#' The function walks the candidate pool greedily, keeping an indicator only
#' if it raises the rank of the anchor block, and errors if it cannot assemble
#' \code{n_anchor} independent rows.
#'
#' Within the pool, candidates are ordered by \emph{ascending} direct-effect
#' SNR, so anchors are the strong-loading indicators with the least evidence of
#' a direct effect. This matters: anchoring an indicator that truly has a
#' direct effect forces \eqn{\gamma_j = 0} where it should not be, biasing
#' \code{B}. Ranking by loading SNR alone can and does pick such an indicator,
#' because a strong loading and a direct effect are not mutually exclusive.
#'
#' The tie-break uses the unidentified \code{Gamma} posterior, which is
#' admittedly circular in principle. It is used only to order candidates that
#' already qualify on loading strength, never to estimate anything, and the
#' alternative -- ignoring the available evidence about which indicators look
#' impure -- is worse. Residual risk remains: refit with the next-best anchors
#' and compare \code{B} if the choice is consequential.
#'
#' @examples
#' set.seed(1)
#' d <- simulate_mimic(n = 150, p_total = 40, p_signal = 15, K = 1,
#'                     has_direct = TRUE)
#' pass1 <- vb_ard(d$Y, d$X, K = 1)
#' anchors <- select_gamma_anchor(pass1)
#' fit <- vb_ard(d$Y, d$X, K = 1, gamma_anchor = anchors)
#' fit$gamma_anchor
#'
#' @seealso \code{\link{vb_ard}}, \code{\link{biomimic}}
#' @export
select_gamma_anchor <- function(fit, n_anchor = NULL, pool = NULL,
                                tol = 1e-8) {
    stopifnot(inherits(fit, "vb_ard"))
    K <- fit$K
    if (is.null(n_anchor)) n_anchor <- K
    n_anchor <- as.integer(n_anchor)
    if (n_anchor < K) {
        warning(sprintf(
            "n_anchor = %d is below K = %d; B will remain unidentified.",
            n_anchor, K), call. = FALSE)
    }
    if (is.null(pool)) pool <- max(3L * n_anchor, 10L)
    pool <- min(as.integer(pool), length(fit$var_names))

    ranking <- rank_variables(fit, type = "lambda")
    ord     <- match(ranking$variable[seq_len(pool)], fit$var_names)
    Lam     <- fit$mu_lambda

    # Prefer, among these strong-loading candidates, the ones with the weakest
    # direct-effect evidence: anchoring a genuinely impure indicator biases B.
    if (!is.null(fit$snr_gamma)) {
        cov_cols <- setdiff(colnames(fit$snr_gamma),
                            c("intercept", "(Intercept)"))
        if (length(cov_cols)) {
            g <- apply(fit$snr_gamma[ord, cov_cols, drop = FALSE], 1,
                       function(r) max(r, na.rm = TRUE))
            ord <- ord[order(g)]
        }
    }

    keep <- integer(0)
    for (i in ord) {
        cand <- c(keep, i)
        block <- Lam[cand, , drop = FALSE]
        # rank via SVD: a new row must add a non-negligible singular value
        sv <- svd(block, nu = 0, nv = 0)$d
        if (min(sv) > tol * max(sv, 1)) keep <- cand
        if (length(keep) == n_anchor) break
    }

    if (length(keep) < n_anchor) {
        stop(sprintf(
            paste0("Could not find %d indicators with linearly independent ",
                   "loadings (found %d). Reduce K or n_anchor."),
            n_anchor, length(keep)), call. = FALSE)
    }

    .log_info("select_gamma_anchor: %s", paste(fit$var_names[keep],
                                               collapse = ", "))
    fit$var_names[keep]
}
