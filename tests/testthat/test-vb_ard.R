# --- Input validation tests ---

test_that("vb_ard rejects non-numeric input", {
    expect_error(vb_ard(data.frame(a = letters[1:5]), matrix(1, 5, 2)))
})

test_that("vb_ard rejects mismatched dimensions", {
    Y <- matrix(rnorm(50), 10, 5)
    X <- matrix(rnorm(18), 9, 2)  # wrong n
    expect_error(vb_ard(Y, X))
})

test_that("rank_variables requires vb_ard object", {
    expect_error(rank_variables(list()), "vb_ard")
})

test_that("top_k requires positive k", {
    expect_error(top_k(structure(list(), class = "vb_ard"), k = 0))
})

test_that("build_lavaan_model rejects invalid type", {
    expect_error(build_lavaan_model(
        structure(list(), class = "biomimic_selection"),
        type = "invalid"
    ))
})


# --- VB-ARD K=1 integration test ---

test_that("vb_ard K=1 recovers signal variables", {
    set.seed(123)
    n <- 200; p <- 30; q <- 2; K <- 1

    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    B <- c(0, 1.0)
    Lambda_true <- c(rep(1, 10), rep(0, 20))
    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    fit <- vb_ard(Y, X, K = 1, max_iter = 200, verbose = FALSE)

    # Should converge
    expect_true(fit$converged)

    # Should return correct dimensions
    expect_equal(dim(fit$mu_lambda), c(p, K))
    expect_equal(dim(fit$mu_B), c(q, K))
    expect_equal(dim(fit$mu_eta), c(n, K))
    expect_equal(length(fit$sigma_lambda), p)
    expect_equal(fit$K, K)
    expect_equal(fit$n, n)
    expect_equal(fit$p, p)
    expect_equal(fit$q, q)

    # SNR should rank signal variables higher than noise
    ranking <- rank_variables(fit, type = "lambda")
    expect_equal(nrow(ranking), p)

    # Top 10 should contain mostly signal variables (V1-V10)
    top10 <- ranking$variable[1:10]
    signal_vars <- paste0("V", 1:10)
    recall <- sum(top10 %in% signal_vars) / 10
    expect_gte(recall, 0.7)  # at least 7/10 recovered
})


# --- VB-ARD K=1 without Gamma ---

test_that("vb_ard works without Gamma", {
    set.seed(42)
    n <- 100; p <- 20
    X <- cbind(1, rnorm(n))
    B <- c(0, 0.5)
    Lambda_true <- c(rep(0.8, 8), rep(0, 12))
    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = 0.7), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    fit <- vb_ard(Y, X, K = 1, use_gamma = FALSE, max_iter = 200)

    expect_true(fit$converged)
    expect_null(fit$mu_gamma)
    expect_null(fit$snr_gamma)

    # Should still rank correctly
    ranking <- rank_variables(fit)
    expect_equal(nrow(ranking), p)
})


# --- VB-ARD with Gamma extension ---

test_that("vb_ard Gamma extension detects direct effects", {
    set.seed(456)
    n <- 300; p <- 25; q <- 2

    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    B <- c(0, 1.0)
    Lambda_true <- c(rep(1, 10), rep(0, 15))

    # Direct effect on variables 11-13 (not loaded on factor)
    Gamma_true <- matrix(0, p, q)
    Gamma_true[11:13, 2] <- 1.5  # strong direct effect of group

    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + X %*% t(Gamma_true) +
         matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    fit <- vb_ard(Y, X, K = 1, use_gamma = TRUE, max_iter = 300)

    expect_true(fit$converged)
    expect_false(is.null(fit$snr_gamma))
    expect_equal(dim(fit$snr_gamma), c(p, q))

    # Gamma ranking should flag V11-V13
    gamma_rank <- rank_variables(fit, type = "gamma")
    top5_gamma <- gamma_rank$variable[1:5]
    direct_vars <- paste0("V", 11:13)
    expect_gte(sum(top5_gamma %in% direct_vars), 2)
})


# --- Top-k selection ---

test_that("top_k returns correct structure", {
    set.seed(789)
    n <- 150; p <- 30
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    B <- c(0, 1)
    Lambda_true <- c(rep(1, 15), rep(0, 15))
    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = 0.7), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    fit <- vb_ard(Y, X, K = 1, max_iter = 200)
    sel <- top_k(fit, k = 10)

    expect_s3_class(sel, "biomimic_selection")
    expect_equal(length(sel$variables), 10)
    expect_true(is.data.frame(sel$data))
    expect_equal(nrow(sel$data), n)
    # data should have 10 Y columns + covariates (minus intercept)
    expect_true("group" %in% names(sel$data))
    expect_equal(sel$K, 1)
    expect_null(sel$indicator_groups)  # K=1 has no groups
})


# --- ELBO monotonicity ---

test_that("ELBO stabilises (convergence)", {
    set.seed(101)
    n <- 100; p <- 15
    X <- cbind(1, rnorm(n))
    Lambda_true <- c(rep(1, 8), rep(0, 7))
    h <- X %*% c(0, 0.5) + rnorm(n)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = 0.7), n, p)

    fit <- vb_ard(Y, X, K = 1, max_iter = 300, use_gamma = FALSE)

    # ELBO should stabilise (converge)
    expect_true(fit$converged)
    # Relative change in last iterations should be small
    elbo <- fit$elbo
    last_change <- abs(elbo[length(elbo)] - elbo[length(elbo) - 1]) /
                   (abs(elbo[length(elbo)]) + 1e-10)
    expect_lt(last_change, 1e-4)
})


# --- K=2 test ---

test_that("vb_ard K=2 runs and produces correct dimensions", {
    set.seed(202)
    n <- 200; p <- 20; K <- 2
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")

    B <- matrix(c(0, 1, 0, 0.5), nrow = 2, ncol = 2)
    Lambda_true <- matrix(0, p, K)
    Lambda_true[1:10, 1] <- 1
    Lambda_true[11:20, 2] <- 0.8
    h <- X %*% B + matrix(rnorm(n * K), n, K)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = 0.5), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    fit <- vb_ard(Y, X, K = 2, use_ard = FALSE, max_iter = 200)

    expect_equal(dim(fit$mu_lambda), c(p, K))
    expect_equal(dim(fit$mu_B), c(2, K))
    expect_equal(dim(fit$snr_lambda), c(p, K))

    # With K=2, ranking should still work
    ranking <- rank_variables(fit)
    expect_true("snr_k1" %in% names(ranking))
    expect_true("snr_k2" %in% names(ranking))
})
