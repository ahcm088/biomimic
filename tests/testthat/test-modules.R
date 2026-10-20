# Tests for the module front ends and as_selection()

make_two_modules <- function(n = 200, seed = 1) {
    set.seed(seed)
    g <- rbinom(n, 1, 0.5)
    h1 <- 0.8 * g + rnorm(n)
    h2 <- rnorm(n)
    A <- sapply(1:8, function(j) 0.8 * h1 + rnorm(n, sd = 0.6))
    B <- sapply(1:8, function(j) 0.8 * h2 + rnorm(n, sd = 0.6))
    N <- matrix(rnorm(n * 24), n, 24)
    Y <- cbind(A, B, N)
    colnames(Y) <- c(paste0("a", 1:8), paste0("b", 1:8), paste0("n", 1:24))
    list(Y = Y, g = g)
}

test_that("modules_from_sets ranks genes and drops small sets", {
    d <- make_two_modules()
    mods <- modules_from_sets(list(A = paste0("a", 1:8), B = paste0("b", 1:8),
                                   tiny = c("a1", "zz")), d$Y)
    expect_s3_class(mods, "biomimic_modules")
    expect_named(mods$modules, c("A", "B"))
    expect_setequal(mods$modules$A$genes, paste0("a", 1:8))
    expect_true(all(diff(mods$modules$A$membership) <= 0))
    expect_gt(min(mods$modules$A$membership), 0.5)
})

test_that("annotate_modules finds the module with most markers", {
    d <- make_two_modules()
    mods <- modules_from_sets(list(X = paste0("a", 1:8), Y = paste0("b", 1:8)), d$Y)
    an <- annotate_modules(mods, list(first = c("a1", "a2", "a3"), second = c("b5", "b6")))
    expect_equal(an$module, c("X", "Y"))
    expect_equal(an$hits, c(3, 2))
})

test_that("as_selection builds a selection that fit_sem accepts (one and two factors)", {
    d <- make_two_modules()
    cov <- data.frame(group = d$g)
    s1 <- as_selection(d$Y, paste0("a", 1:6), covariates = cov)
    expect_s3_class(s1, "biomimic_selection")
    expect_null(s1$fit)
    f1 <- fit_sem(s1, type = "mimic")
    pe <- lavaan::parameterEstimates(f1$lavaan_fit)
    b <- pe$est[pe$op == "~" & pe$rhs == "group"]
    expect_gt(b, 0)

    s2 <- as_selection(d$Y, list(IFN = paste0("a", 1:6), LYM = paste0("b", 1:6)),
                       covariates = cov)
    expect_equal(s2$K, 2L)
    expect_length(s2$indicator_groups, 2L)
    f2 <- fit_sem(s2, type = "multi_mimic")
    expect_equal(f2$factor_labels, c("IFN", "LYM"))
    expect_true(lavaan::lavInspect(f2$lavaan_fit, "converged"))
})

test_that("as_selection sanitises names and keeps the original labels", {
    set.seed(2)
    Y <- matrix(rnorm(80 * 4), 80, 4, dimnames = list(NULL, c("HLA-A", "KRTAP13-2", "g3", "g4")))
    s <- as_selection(Y, c("HLA-A", "KRTAP13-2", "g3"))
    expect_equal(unname(s$labels), c("HLA-A", "KRTAP13-2", "g3"))
    expect_true(all(make.names(s$variables) == s$variables))
})

test_that("as_selection rejects overlapping or too small panels", {
    d <- make_two_modules()
    expect_error(as_selection(d$Y, list(A = c("a1", "a2"), B = c("a2", "b1"))), "overlap")
    expect_error(as_selection(d$Y, list(A = c("a1", "missing"))), "at least 2")
})

test_that("select_panels takes the top genes present in another cohort", {
    d <- make_two_modules()
    mods <- modules_from_sets(list(A = paste0("a", 1:8)), d$Y)
    Y2 <- d$Y[, setdiff(colnames(d$Y), "a1")]
    s <- select_panels(mods, "A", k = 5, Y = Y2)
    expect_length(s$variables, 5L)
    expect_false("a1" %in% s$labels)
})

test_that("modules_ard returns one module per factor", {
    d <- make_two_modules(n = 120)
    X <- cbind(intercept = 1, group = d$g)
    vb <- vb_ard(scale(d$Y), X, K = 2L, max_iter = 200L, verbose = FALSE)
    mods <- modules_ard(vb, n = 10)
    expect_length(mods$modules, 2L)
    expect_true(all(lengths(lapply(mods$modules, `[[`, "genes")) == 10L))
})

test_that("modules_wgcna finds the planted modules", {
    skip_if_not_installed("WGCNA")
    skip_on_cran()
    d <- make_two_modules(n = 150)
    mods <- suppressMessages(modules_wgcna(d$Y, n_genes = 40, min_size = 5,
                                           raw_var = apply(d$Y, 2, var)))
    an <- annotate_modules(mods, list(A = paste0("a", 1:8), B = paste0("b", 1:8)), top = 8)
    expect_gte(min(an$hits), 6)
    expect_false(an$module[1] == an$module[2])
})
