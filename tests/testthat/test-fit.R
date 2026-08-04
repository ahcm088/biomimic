# --- Stage 3 integration tests ---
# These tests require lavaan to be installed.

# Helper: generate simulated data for a K=1 MIMIC model
.sim_mimic_data <- function(n = 200, p = 10, seed = 42) {
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


test_that("fit_sem creates biomimic_fit with correct structure", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    expect_s3_class(fit, "biomimic_fit")
    expect_true(inherits(fit$lavaan_fit, "lavaan"))
    expect_equal(fit$model_type, "mimic")
    expect_equal(fit$K, 1L)
    expect_equal(length(fit$factor_labels), 1)
    expect_equal(length(fit$indicator_labels), 10)
})


test_that("summary.biomimic_fit produces output without errors", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    # summary should run without error and return biomimic_summary
    result <- expect_output(summary(fit), "biomimic SEM Results")
    expect_s3_class(result, "biomimic_summary")
    expect_true(nrow(result$loadings) > 0)
    expect_true(nrow(result$structural) > 0)
})


test_that("coef.biomimic_fit returns tidy data frame", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    # All coefficients
    cf <- coef(fit)
    expect_true(is.data.frame(cf))
    expect_true("estimate" %in% names(cf))
    expect_true("pvalue" %in% names(cf))
    expect_true("from_label" %in% names(cf))

    # Loadings only
    cf_load <- coef(fit, type = "loadings")
    expect_true(all(cf_load$type == "loading"))

    # Structural only
    cf_struct <- coef(fit, type = "structural")
    expect_true(all(cf_struct$type == "structural"))
})


test_that("confint.biomimic_fit returns intervals", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    ci <- confint(fit)
    expect_true(is.data.frame(ci))
    expect_true("ci_lower" %in% names(ci))
    expect_true("ci_upper" %in% names(ci))
    expect_true(all(ci$ci_lower <= ci$ci_upper))
})


test_that("goodness_of_fit returns valid indices", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    gof <- goodness_of_fit(fit)
    expect_s3_class(gof, "biomimic_gof")
    expect_true(gof$rmsea >= 0)
    expect_true(gof$cfi >= 0 && gof$cfi <= 1)
    expect_true(gof$srmr >= 0)
    expect_true(gof$df > 0)
})


test_that("residual_analysis returns residuals and summary", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    ra <- residual_analysis(fit)
    expect_true(is.matrix(ra$residuals))
    expect_true(is.data.frame(ra$summary))
    expect_true("abs_residual" %in% names(ra$summary))
})


test_that("mahalanobis_distance detects outliers", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    md <- mahalanobis_distance(fit)
    expect_true(is.data.frame(md))
    expect_equal(nrow(md), d$n)
    expect_true(all(md$mahal_dist >= 0))
    expect_true("outlier" %in% names(md))
})


test_that("predict_scores returns factor scores", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    scores <- predict_scores(fit, append_data = FALSE)
    expect_equal(nrow(scores), d$n)
    expect_equal(ncol(scores), 1)  # K=1

    # With appended data
    scores_full <- predict_scores(fit, append_data = TRUE)
    expect_true(ncol(scores_full) > 1)
    expect_equal(nrow(scores_full), d$n)
})


test_that("predict_observed returns predictions and residuals", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)
    fit <- fit_sem(sel, type = "mimic")

    po <- predict_observed(fit)
    expect_true(is.matrix(po$predicted))
    expect_true(is.matrix(po$residuals))
    expect_equal(dim(po$predicted), c(d$n, 10))
    expect_equal(dim(po$residuals), c(d$n, 10))
})


test_that("biomimic() runs the full pipeline", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data(n = 150, p = 20)
    fit <- biomimic(d$Y, d$X, K = 1, k = 10, max_iter = 100)

    expect_s3_class(fit, "biomimic_fit")
    expect_true(inherits(fit$lavaan_fit, "lavaan"))
    expect_equal(length(fit$selection$variables), 10)
})


test_that("biology labels are applied correctly", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data(n = 150, p = 10)
    vb <- vb_ard(d$Y, d$X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)

    labels <- setNames(paste0("Protein_", 1:10), paste0("V", 1:10))
    fit <- fit_sem(sel, type = "mimic", indicator_labels = labels,
                   factor_labels = "Immune activation")

    expect_equal(fit$factor_labels, "Immune activation")
    expect_equal(unname(fit$indicator_labels["V1"]), "Protein_1")

    cf <- coef(fit, type = "loadings")
    expect_true(all(grepl("Protein_", cf$from_label)))
})
