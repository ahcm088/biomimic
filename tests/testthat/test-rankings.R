# Tests for the limma and CATE rankings

make_ranked <- function(n = 120, p = 200, seed = 5) {
    set.seed(seed)
    g <- rbinom(n, 1, 0.5)
    h <- rnorm(n)                                   # hidden factor
    Y <- sapply(seq_len(p), function(j) 0.5 * h + rnorm(n))
    Y[, 1:10] <- Y[, 1:10] + 1.0 * g                # true group effects
    colnames(Y) <- paste0("v", seq_len(p))
    list(Y = Y, g = g)
}

test_that("rank_limma puts the true effects on top", {
    skip_if_not_installed("limma")
    d <- make_ranked()
    rk <- rank_limma(d$Y, data.frame(group = d$g))
    expect_named(rk, c("variable", "snr", "rank", "estimate", "t", "p", "adj_p"))
    expect_gte(sum(paste0("v", 1:10) %in% head(rk$variable, 15)), 9)
    expect_true(all(diff(rk$snr) <= 0))
})

test_that("rank_cate puts the true effects on top", {
    skip_if_not_installed("cate")
    skip_on_cran()
    d <- make_ranked()
    rk <- rank_cate(d$Y, data.frame(group = d$g), r = 1)
    expect_equal(attr(rk, "r"), 1)
    expect_gte(sum(paste0("v", 1:10) %in% head(rk$variable, 15)), 9)
})

test_that("a ranking feeds as_selection", {
    skip_if_not_installed("limma")
    d <- make_ranked()
    rk <- rank_limma(d$Y, data.frame(group = d$g))
    sel <- as_selection(d$Y, head(rk$variable, 6), covariates = data.frame(group = d$g))
    expect_s3_class(sel, "biomimic_selection")
})
