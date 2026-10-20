# Tests for the global-level diagnostics

make_global <- function(n = 200, p = 60, seed = 3) {
    set.seed(seed)
    g <- rbinom(n, 1, 0.5)
    level <- 0.8 * g + rnorm(n)                 # each sample's global level
    Y <- sapply(seq_len(p), function(j) 0.7 * level + rnorm(n))
    colnames(Y) <- paste0("v", seq_len(p))
    list(Y = Y, g = g, level = level)
}

test_that("global_level_check flags a panel that follows the global level", {
    d <- make_global()
    gl <- global_level_check(d$Y, paste0("v", 1:10), d$g)
    expect_s3_class(gl, "biomimic_global")
    expect_gt(gl$r_global, 0.8)
    # the group effect on the panel disappears once the global level is in
    expect_gt(gl$group_effect[["estimate"]], 0.3)
    expect_lt(abs(gl$group_effect_adjusted[["estimate"]]), 0.2)
})

test_that("a panel independent of the global level is not flagged", {
    d <- make_global()
    set.seed(4)
    own <- sapply(1:10, function(j) 0.8 * d$g + rnorm(nrow(d$Y)))
    colnames(own) <- paste0("s", 1:10)
    gl <- global_level_check(cbind(d$Y, own), colnames(own), d$g)
    expect_lt(gl$r_global, 0.6)
    expect_gt(gl$group_effect_adjusted[["estimate"]], 0.3)
})

test_that("panel_scores and external_check recover an external measure", {
    d <- make_global()
    sc <- panel_scores(d$Y, list(P = paste0("v", 1:10)))
    ex <- external_check(sc, data.frame(level = d$level))
    expect_equal(nrow(ex), 1L)
    expect_gt(ex$rho, 0.7)
})

test_that("normalise_arrays removes the per-array level", {
    d <- make_global()
    Ym <- normalise_arrays(d$Y, "median")
    expect_lt(abs(stats::cor(apply(Ym, 1, median), d$level)), 0.05)
    Yq <- normalise_arrays(d$Y, "quantile")
    expect_equal(dim(Yq), dim(d$Y))
    expect_lt(stats::sd(apply(Yq, 1, mean)), 1e-8)
})

test_that("coordination_diagnostic reports r_global", {
    sim <- simulate_mimic(n = 100L, p_total = 80L, p_signal = 12L, K = 1L, seed = 2)
    fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 150L, verbose = FALSE)
    cd <- coordination_diagnostic(fit, n_null = 20L)
    expect_true("r_global" %in% names(cd$channels))
})
