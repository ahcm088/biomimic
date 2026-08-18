# Tests for parametric-bootstrap calibration of the gamma threshold

make_fit <- function(seed = 11, n = 120, p = 25) {
    d <- simulate_mimic(n = n, p_total = p, p_signal = 10L, K = 1L,
                        B_group = 1.0, gen_dist = "normal",
                        has_direct = FALSE, seed = seed)
    list(fit = vb_ard(d$Y, d$X, K = 1L, max_iter = 200), d = d)
}

test_that("calibrate_gamma_threshold returns a sane calibration object", {
    f <- make_fit()
    sel <- top_k(f$fit, k = 10)
    cal <- calibrate_gamma_threshold(f$fit, sel$variables,
                                     n_boot = 25, n_cores = 1, seed = 7)
    expect_s3_class(cal, "biomimic_gamma_calibration")
    expect_true(is.finite(cal$threshold) && cal$threshold > 0)
    expect_equal(cal$n_boot, 25)
    expect_equal(cal$level, 0.95)
    expect_length(cal$boot_max, 25)
    # threshold is the level-quantile of the maxima
    expect_equal(cal$threshold,
                 as.numeric(quantile(cal$boot_max, 0.95)))
    expect_output(print(cal), "parametric bootstrap")
})

test_that("calibration is reproducible for a fixed seed", {
    f <- make_fit()
    sel <- top_k(f$fit, k = 10)
    c1 <- calibrate_gamma_threshold(f$fit, sel$variables,
                                    n_boot = 10, n_cores = 1, seed = 3)
    c2 <- calibrate_gamma_threshold(f$fit, sel$variables,
                                    n_boot = 10, n_cores = 1, seed = 3)
    expect_identical(c1$boot_max, c2$boot_max)
})

test_that("top_k accepts gamma_threshold = 'bootstrap'", {
    f <- make_fit()
    sel <- top_k(f$fit, k = 10, gamma_threshold = "bootstrap",
                 n_boot = 25, n_cores = 1)
    expect_s3_class(sel$gamma_calibration, "biomimic_gamma_calibration")
    expect_true(is.numeric(sel$gamma_threshold))
    expect_equal(sel$gamma_threshold, sel$gamma_calibration$threshold)
    # true gamma is 0 everywhere: an FWER-95 threshold should flag nothing
    # in most datasets; allow at most 1 (it is a 5% event by construction)
    expect_lte(length(sel$direct_effects), 1L)
})

test_that("numeric threshold path records the value and no calibration", {
    f <- make_fit()
    sel <- top_k(f$fit, k = 10, gamma_threshold = 2.5)
    expect_null(sel$gamma_calibration)
    expect_equal(sel$gamma_threshold, 2.5)
})

test_that("biomimic() passes gamma_threshold through to top_k", {
    f <- make_fit()
    d <- f$d
    bf <- biomimic(d$Y, d$X, K = 1, k = 10, gamma_anchor = NULL,
                   gamma_threshold = 99)
    expect_equal(bf$selection$gamma_threshold, 99)
    expect_length(bf$selection$direct_effects, 0L)
})

test_that("calibration refuses fits without gamma and unknown variables", {
    d <- simulate_mimic(n = 80, p_total = 15, p_signal = 8L, K = 1L,
                        seed = 5)
    f0 <- vb_ard(d$Y, d$X, K = 1L, use_gamma = FALSE, max_iter = 100)
    expect_error(calibrate_gamma_threshold(f0, colnames(d$Y)[1:5]),
                 "use_gamma")
    f1 <- vb_ard(d$Y, d$X, K = 1L, max_iter = 100)
    expect_error(calibrate_gamma_threshold(f1, c("V1", "NOPE")),
                 "not found")
})

test_that("anchored fit calibrates on non-anchored variables only", {
    f <- make_fit()
    a <- select_gamma_anchor(f$fit)
    fa <- vb_ard(f$d$Y, f$d$X, K = 1L, gamma_anchor = a, max_iter = 200)
    sel <- top_k(fa, k = 10)
    cal <- calibrate_gamma_threshold(fa, sel$variables,
                                     n_boot = 10, n_cores = 1)
    expect_false(a %in% cal$variables)
    expect_true(is.finite(cal$threshold))
})
