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
#' @param use_ard Logical, whether to use ARD priors on loadings (default TRUE).
#'   Recommended FALSE for K > 1.
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
#' @param verbose Logical, print progress (default FALSE).
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
                                       E_psi_inv, E_tau, n, p, K, use_ard)
        mu_lambda    <- res_lam$mu_lambda
        sigma_lambda <- res_lam$sigma_lambda

        # --- Update q(Gamma_j) if requested ---
        if (use_gamma) {
            res_gam <- .cavi_update_gamma(Y, X, mu_eta, mu_lambda, E_psi_inv,
                                          E_xi, n, p, q, XtX)
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
            E_xi <- .cavi_update_xi(mu_gamma, sigma_gamma, a_xi, b_xi, p, q)
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
                                      n, p, q, K, use_ard, use_gamma)

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
                                E_psi_inv, E_tau, n, p, K, use_ard) {
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
        if (use_ard) {
            prec <- E_psi_inv[j] * E_HtH + diag(E_tau[j, ], nrow = K)
        } else {
            prec <- E_psi_inv[j] * E_HtH + diag(E_tau[j, ], nrow = K)
        }
        sigma_lambda[[j]] <- solve(prec)

        # mu_lambda_j = psi_j^{-1} Sigma_lambda_j E[h]' y_adj_j
        mu_lambda[j, ] <- E_psi_inv[j] * sigma_lambda[[j]] %*% HtY[, j]
    }

    list(mu_lambda = mu_lambda, sigma_lambda = sigma_lambda)
}


#' @noRd
.cavi_update_gamma <- function(Y, X, mu_eta, mu_lambda, E_psi_inv,
                               E_xi, n, p, q, XtX) {
    # Residual: Y - H Lambda'
    Y_adj <- Y - mu_eta %*% t(mu_lambda)  # n x p

    # XtY_adj = X' Y_adj -> q x p
    XtY <- crossprod(X, Y_adj)

    mu_gamma <- matrix(0, p, q)
    sigma_gamma <- vector("list", p)

    for (j in seq_len(p)) {
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
.cavi_update_xi <- function(mu_gamma, sigma_gamma, a_xi, b_xi, p, q) {
    a_hat <- a_xi + 0.5
    E_xi <- matrix(0, p, q)
    for (j in seq_len(p)) {
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
                          n, p, q, K, use_ard, use_gamma) {

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
        for (j in seq_len(p)) {
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
