# --- Co-ordination diagnostic ---

test_that("coordination_diagnostic requires a vb_ard object", {
    expect_error(coordination_diagnostic(list()), "vb_ard")
})

test_that("coordination_diagnostic rejects a degenerate panel size", {
    sim <- simulate_mimic(n = 60L, p_total = 40L, p_signal = 10L, seed = 1)
    fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 50L)
    expect_error(coordination_diagnostic(fit, k = 2L))
})

test_that("a simulated factor is reported as co-ordinated", {
    # The signal variables share one latent factor by construction, so the
    # lambda channel should land well above the background of random sets.
    sim <- simulate_mimic(n = 150L, p_total = 200L, p_signal = 15L,
                          lambda_strength = 1.5, seed = 42)
    fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 300L)
    d <- coordination_diagnostic(fit, k = 15L, n_null = 50L)

    expect_s3_class(d, "biomimic_regime")
    expect_true(all(c("channel", "n_variables", "mean_abs_r", "background",
                      "ratio") %in% names(d$channels)))
    lam <- d$channels[d$channels$channel == "lambda", ]
    expect_equal(nrow(lam), 1L)
    expect_gt(lam$ratio, 2)
    expect_true(is.finite(lam$mean_abs_r) && lam$mean_abs_r <= 1)

    # Empirical p-value: bounded by 1 / (1 + n_null) and, for a panel this far
    # above the background, at that bound.
    expect_true("p_value" %in% names(d$channels))
    expect_true(all(d$channels$p_value > 0 & d$channels$p_value <= 1))
    expect_equal(lam$p_value, 1 / (1 + length(d$background_all$lambda)))
})

test_that("the reference argument restricts the samples used", {
    sim <- simulate_mimic(n = 120L, p_total = 100L, p_signal = 12L, seed = 7)
    fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 200L)
    ctl <- which(sim$X[, "group"] == 0)

    d_res <- coordination_diagnostic(fit, k = 10L, n_null = 30L)
    d_ctl <- coordination_diagnostic(fit, k = 10L, reference = ctl,
                                     n_null = 30L)

    expect_equal(d_res$reference_n, nrow(sim$Y))
    expect_equal(d_ctl$reference_n, length(ctl))
    expect_true(d_res$residualised)
    expect_false(d_ctl$residualised)
    # Same panels either way: the channels are ranked on the fit, not on the
    # subset, so only the correlation changes.
    expect_identical(d_res$panels$lambda, d_ctl$panels$lambda)
})

test_that("too few reference samples is an error, not a silent NA", {
    sim <- simulate_mimic(n = 60L, p_total = 50L, p_signal = 10L, seed = 3)
    fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 50L)
    expect_error(coordination_diagnostic(fit, k = 5L, reference = 1:3),
                 "at least 5")
})

test_that("the gamma channel is reported only when gamma was fitted", {
    sim <- simulate_mimic(n = 80L, p_total = 60L, p_signal = 10L, seed = 11)
    fit_no <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = FALSE, max_iter = 100L)
    d <- coordination_diagnostic(fit_no, k = 8L, n_null = 20L)
    expect_identical(d$channels$channel, "lambda")
})

test_that("print returns its argument invisibly", {
    sim <- simulate_mimic(n = 80L, p_total = 60L, p_signal = 10L, seed = 5)
    fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 100L)
    d <- coordination_diagnostic(fit, k = 8L, n_null = 20L)
    expect_output(out <- print(d), "co-ordination diagnostic")
    expect_identical(out, d)
})
