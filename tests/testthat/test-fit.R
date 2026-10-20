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


# --- Factor orientation (marker sign) ---

# Helper: B = +1 on a factor that every indicator but V1 loads positively on.
.sim_minority_marker <- function(n = 300, p = 10, seed = 7) {
    set.seed(seed)
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    h <- X %*% c(0, 1) + rnorm(n)
    lam <- c(-1, rep(1, p - 1))
    Y <- h %*% t(lam) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))
    list(Y = Y, X = X)
}

.b_group <- function(fit) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit)
    pe$est[pe$op == "~" & pe$lhs == "F1" & pe$rhs == "group"]
}

test_that("a minority-sign marker reverses B when orientation is off", {
    skip_if_not_installed("lavaan")
    d <- .sim_minority_marker()
    sel <- top_k(vb_ard(d$Y, d$X, K = 1, max_iter = 100), k = 10)
    sel$variables <- c("V1", setdiff(sel$variables, "V1"))

    fit <- fit_sem(sel, type = "mimic", orient = FALSE)
    expect_equal(fit$orientation$marker, "V1")
    expect_false(fit$orientation$aligned)
    expect_lt(.b_group(fit), 0)
})

test_that("fit_sem reorients a minority-sign marker by default", {
    skip_if_not_installed("lavaan")
    d <- .sim_minority_marker()
    sel <- top_k(vb_ard(d$Y, d$X, K = 1, max_iter = 100), k = 10)
    sel$variables <- c("V1", setdiff(sel$variables, "V1"))

    fit <- fit_sem(sel, type = "mimic")
    expect_true(fit$orientation$aligned)
    expect_false(fit$orientation$marker == "V1")
    expect_equal(fit$orientation$replaced, "V1")
    expect_gt(.b_group(fit), 0)
    expect_equal(fit$orientation$n_neg, 1L)
})

test_that("an aligned marker is left alone", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data()
    sel <- top_k(vb_ard(d$Y, d$X, K = 1, max_iter = 100), k = 10)
    fit <- fit_sem(sel, type = "mimic")
    expect_true(fit$orientation$aligned)
    expect_equal(fit$orientation$marker, sel$variables[1])
    expect_true(is.na(fit$orientation$replaced))
})

test_that("the direct model never reorients onto a direct-effect indicator", {
    skip_if_not_installed("lavaan")
    d <- .sim_minority_marker()
    sel <- top_k(vb_ard(d$Y, d$X, K = 1, max_iter = 100), k = 10)
    sel$variables <- c("V1", setdiff(sel$variables, "V1"))
    # Make every majority-sign indicator but one "direct": that one must win.
    sel$direct_effects <- setdiff(sel$variables, c("V1", "V5"))
    fit <- suppressWarnings(fit_sem(sel, type = "direct"))
    expect_equal(fit$orientation$marker, "V5")
})

test_that("select_gamma_anchor picks a majority-sign anchor for K = 1", {
    # V1 has the largest loading, so VB-ARD orients the factor on it, and it
    # is then the only positive loading among ten.
    set.seed(11)
    n <- 300; p <- 10
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    h <- X %*% c(0, 1) + rnorm(n)
    lam <- c(-3, rep(1, p - 1))
    Y <- h %*% t(lam) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, max_iter = 200)
    s <- sign(vb$mu_lambda[, 1])
    maj <- if (sum(s > 0) >= sum(s < 0)) 1 else -1
    a <- select_gamma_anchor(vb)
    expect_equal(unname(s[match(a, vb$var_names)]), maj)
})

test_that("lrt_structural tests every factor path and is marker invariant", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data(n = 300)
    sel <- top_k(vb_ard(d$Y, d$X, K = 1, max_iter = 100), k = 10)
    fit <- fit_sem(sel, type = "mimic")
    r <- lrt_structural(fit)
    expect_equal(r$covariate, "group")
    expect_true(all(c("p_wald", "p_lrt", "std_estimate") %in% names(r)))
    expect_lt(r$p_lrt, 0.001)            # B = 1 in the simulation

    # Same model, another marker: the LRT and standardised B do not move.
    sel2 <- sel
    sel2$variables <- c(sel$variables[5], setdiff(sel$variables,
                                                  sel$variables[5]))
    r2 <- lrt_structural(fit_sem(sel2, type = "mimic"))
    expect_equal(r2$chisq_diff, r$chisq_diff, tolerance = 1e-3)
    expect_equal(r2$std_estimate, r$std_estimate, tolerance = 1e-4)
})

test_that("lrt_structural restricts to the requested covariates", {
    skip_if_not_installed("lavaan")
    d <- .sim_mimic_data(n = 200)
    sel <- top_k(vb_ard(d$Y, d$X, K = 1, max_iter = 100), k = 10)
    fit <- fit_sem(sel, type = "mimic")
    expect_error(lrt_structural(fit, covariates = "nope"), "no structural")
})
