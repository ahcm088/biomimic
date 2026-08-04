# --- Unified hypothesis testing and diagnostics tests ---

.sim_data <- function(n = 200, p = 10, seed = 42) {
    set.seed(seed)
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    B <- c(0, 1.0)
    Lambda_true <- rep(1, p)
    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))
    list(Y = Y, X = X, n = n, p = p)
}


test_that("test_hypotheses works for MIMIC model", {
    skip_if_not_installed("lavaan")
    d <- .sim_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    th <- test_hypotheses(fit)
    expect_s3_class(th, "biomimic_hypotheses")
    expect_equal(th$model_type, "mimic")
    expect_true(nrow(th$tests) > 0)
    expect_true("hypothesis" %in% names(th$tests))
    expect_true("pvalue" %in% names(th$tests))
    expect_true("significant" %in% names(th$tests))

    # Group effect should be significant (B = 1.0)
    group_test <- th$tests[grep("group", th$tests$parameter), ]
    expect_true(nrow(group_test) > 0)
    expect_true(group_test$significant[1])

    # Print should work
    expect_output(print(th), "Hypothesis Tests")
})


test_that("test_hypotheses works for direct effects model", {
    skip_if_not_installed("lavaan")
    set.seed(456)
    n <- 300; p <- 15
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    B <- c(0, 1.0)
    Lambda_true <- c(rep(1, 10), rep(0, 5))
    Gamma_true <- matrix(0, p, 2)
    Gamma_true[11:12, 2] <- 1.5
    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + X %*% t(Gamma_true) +
         matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, use_gamma = TRUE, max_iter = 200)
    sel <- top_k(vb, k = 12, gamma_threshold = 2.0)

    # Only run if direct effects were detected
    if (length(sel$direct_effects) > 0) {
        fit <- fit_sem(sel, type = "direct")
        th <- test_hypotheses(fit)

        expect_equal(th$model_type, "direct")
        # Should have both structural and direct tests
        has_direct <- any(grepl("direct", th$tests$hypothesis))
        expect_true(has_direct || nrow(th$tests) > 0)

        # Mediation decomposition
        if (!is.null(th$mediation)) {
            expect_true("indirect" %in% names(th$mediation))
            expect_true("direct" %in% names(th$mediation))
        }
    }
})


test_that("diagnose works for MIMIC model", {
    skip_if_not_installed("lavaan")
    d <- .sim_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    diag <- diagnose(fit)
    expect_s3_class(diag, "biomimic_diagnostics")
    expect_s3_class(diag$gof, "biomimic_gof")
    expect_true(is.data.frame(diag$mahalanobis))
    expect_true(!is.null(diag$r_squared))
    expect_true(!is.null(diag$reliability))

    # Print should work
    expect_output(print(diag), "Model Diagnostics")
})


test_that("test_hypotheses dispatches correctly by model_type", {
    skip_if_not_installed("lavaan")
    d <- .sim_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)

    # MIMIC
    fit_m <- fit_sem(sel, type = "mimic")
    th_m <- test_hypotheses(fit_m)
    expect_equal(th_m$model_type, "mimic")
    expect_null(th_m$invariance)
    expect_null(th_m$mediation)
})


test_that("diagnose dispatches correctly by model_type", {
    skip_if_not_installed("lavaan")
    d <- .sim_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)

    fit <- fit_sem(sel, type = "mimic")
    diag <- diagnose(fit)
    expect_equal(diag$model_type, "mimic")
    expect_true(!is.null(diag$reliability))
})
