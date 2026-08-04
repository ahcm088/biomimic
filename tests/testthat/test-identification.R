# --- Identification checks tests ---

test_that("build_lavaan_model rejects mimic with < 3 indicators", {
    sel <- structure(
        list(
            variables = c("y1", "y2"),  # only 2
            indicator_groups = NULL,
            direct_effects = character(0),
            K = 1L,
            fit = list(X = matrix(0, 10, 2,
                                  dimnames = list(NULL, c("(Intercept)", "group"))))
        ),
        class = "biomimic_selection"
    )
    expect_error(build_lavaan_model(sel, type = "mimic"),
                 "at least 3 indicators")
})


test_that("build_lavaan_model warns for multi_mimic with 2 indicators per factor", {
    sel <- structure(
        list(
            variables = c("y1", "y2", "y3", "y4"),
            indicator_groups = list(1:2, 3:4),  # 2 per factor
            direct_effects = character(0),
            K = 2L,
            fit = list(X = matrix(0, 10, 2,
                                  dimnames = list(NULL, c("(Intercept)", "group"))))
        ),
        class = "biomimic_selection"
    )
    expect_warning(build_lavaan_model(sel, type = "multi_mimic"),
                   "only 2 indicators")
})


test_that("build_lavaan_model errors for multi_mimic with 1 indicator", {
    sel <- structure(
        list(
            variables = c("y1", "y2", "y3"),
            indicator_groups = list(1L, 2:3),  # factor 1 has only 1
            direct_effects = character(0),
            K = 2L,
            fit = list(X = matrix(0, 10, 2,
                                  dimnames = list(NULL, c("(Intercept)", "group"))))
        ),
        class = "biomimic_selection"
    )
    expect_error(build_lavaan_model(sel, type = "multi_mimic"),
                 "1 indicator")
})


test_that("build_lavaan_model passes for well-specified models", {
    sel <- structure(
        list(
            variables = c("y1", "y2", "y3", "y4", "y5"),
            indicator_groups = NULL,
            direct_effects = character(0),
            K = 1L,
            fit = list(X = matrix(0, 10, 2,
                                  dimnames = list(NULL, c("(Intercept)", "group"))))
        ),
        class = "biomimic_selection"
    )
    # Should not error or warn
    syntax <- build_lavaan_model(sel, type = "mimic")
    expect_true(grepl("F1 =~", syntax))
})


test_that("check_identification works for identified MIMIC model", {
    skip_if_not_installed("lavaan")

    set.seed(42)
    n <- 200; p <- 10
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    B <- c(0, 1.0)
    Lambda_true <- rep(1, p)
    h <- X %*% B + rnorm(n)
    Y <- h %*% t(Lambda_true) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    id <- check_identification(fit)
    expect_s3_class(id, "biomimic_identification")
    expect_true(id$identified)
    expect_true(id$df > 0)
    expect_true(id$converged)
    expect_true(id$n_free < id$n_equations)

    # Print should work
    expect_output(print(id), "IDENTIFIED")
})
