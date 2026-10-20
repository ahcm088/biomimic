# Tests for the structural-model builders

test_that("fit_bifactor separates a specific response from the global level", {
    set.seed(11)
    n <- 300; g <- rbinom(n, 1, 0.5)
    level <- rnorm(n)                             # global level, not shifted
    s <- 1.0 * g + rnorm(n)                       # specific response, shifted
    Y <- sapply(1:30, function(j) 0.7 * level + rnorm(n, sd = 0.7))
    S <- sapply(1:4, function(j) 0.6 * level + 0.7 * s + rnorm(n, sd = 0.5))
    Y <- cbind(Y, S); colnames(Y) <- c(paste0("v", 1:30), paste0("s", 1:4))
    bf <- fit_bifactor(Y, specific = list(S1 = paste0("s", 1:4)), general = paste0("v", 1:30),
                       covariates = data.frame(group = g), lrt = TRUE)
    expect_s3_class(bf, "biomimic_bifactor")
    expect_true(bf$converged)
    e <- bf$effects
    expect_gt(e$estimate[e$factor == "S1"], 0.5)
    expect_lt(e$p[e$factor == "S1"], 1e-4)
    expect_lt(abs(e$estimate[e$factor == "G"]), 0.3)
    expect_true(is.finite(e$p_lrt[e$factor == "S1"]))
})

test_that("fit_bifactor constrains two-indicator specific factors", {
    set.seed(12)
    n <- 300; g <- rbinom(n, 1, 0.5); level <- rnorm(n); s <- g + rnorm(n)
    Y <- cbind(sapply(1:20, function(j) 0.7 * level + rnorm(n, sd = 0.7)),
               sapply(1:2, function(j) 0.5 * level + 0.8 * s + rnorm(n, sd = 0.4)))
    colnames(Y) <- c(paste0("v", 1:20), "RO60", "SSB")
    bf <- fit_bifactor(Y, list(S_RoLa = c("RO60", "SSB")), general = 20L,
                       covariates = data.frame(group = g), lrt = FALSE)
    expect_true(bf$converged)
    expect_match(bf$syntax, "l_S_RoLa\\*RO60")
    expect_gt(bf$effects$estimate[bf$effects$factor == "S_RoLa"], 0.4)
})

test_that("fit_mediation recovers latent and observed indirect effects", {
    set.seed(13)
    n <- 1000; x <- rnorm(n)
    m1 <- 0.5 * x + rnorm(n)                      # latent mediator
    m2 <- -0.4 * x + rnorm(n)                     # observed mediator
    y <- 0.4 * m1 - 0.5 * m2 + 0.1 * x + rnorm(n)
    ind <- sapply(1:5, function(j) 0.8 * m1 + rnorm(n, sd = 0.6))
    d <- data.frame(ind, m2 = m2, x = x, y = y); names(d)[1:5] <- paste0("i", 1:5)
    md <- fit_mediation(d, exposure = "x", mediators = list(M1 = paste0("i", 1:5), C = "m2"),
                        outcome = "y")
    e <- md$effects
    expect_equal(e$mediator, c("M1", "C"))
    expect_gt(e$indirect[1], 0); expect_gt(e$indirect[2], 0)
    expect_lt(max(e$p_lrt), 1e-6)
    # the total effect equals the regression of the outcome on the exposure
    expect_equal(unname(md$total[["estimate"]]), unname(coef(lm(y ~ x))[2]), tolerance = 0.02)
})

test_that("decompose_effects separates factor paths from own effects", {
    set.seed(14)
    n <- 400; g <- rbinom(n, 1, 0.5)
    f1 <- 0.8 * g + rnorm(n); f2 <- -0.6 * g + rnorm(n)
    P1 <- sapply(1:6, function(j) 0.8 * f1 + rnorm(n, sd = 0.6))
    P2 <- sapply(1:6, function(j) 0.8 * f2 + rnorm(n, sd = 0.6))
    t_via <- 0.7 * f1 + rnorm(n, sd = 0.7)        # all through factor 1
    t_own <- 0.5 * g + rnorm(n)                   # own effect only
    Y <- cbind(P1, P2, t_via, t_own)
    colnames(Y) <- c(paste0("a", 1:6), paste0("b", 1:6), "t_via", "t_own")
    Y <- scale(Y)
    de <- decompose_effects(Y, list(IFN = paste0("a", 1:6), LYM = paste0("b", 1:6)),
                            targets = c("t_via", "t_own"), group = g, n_cores = 1L)
    expect_equal(de$feature, c("t_via", "t_own"))
    expect_true(all(c("via_IFN", "via_IFN_se", "via_LYM", "own", "own_se") %in% names(de)))
    expect_lt(abs(de$own[1]), 2.5 * de$own_se[1])
    expect_gt(de$via_IFN[1], 4 * de$via_IFN_se[1])
    expect_gt(de$own[2], 3 * de$own_se[2])
    expect_lt(abs(de$via_IFN[2]), 2.5 * de$via_IFN_se[2])
    expect_equal(de$total, de$via_IFN + de$via_LYM + de$own, tolerance = 1e-6)
})

test_that("within_between splits repeated measures", {
    d <- data.frame(id = rep(1:3, each = 2), x = c(1, 3, 2, 2, 5, 7))
    w <- within_between(d, "x", "id")
    expect_equal(w$x_b, c(2, 2, 2, 2, 6, 6))
    expect_equal(w$x_w, c(-1, 1, 0, 0, -1, 1))
})

test_that("fit_sem with cluster estimates a within-subject effect", {
    set.seed(17)
    N <- 80; T <- 4; id <- rep(seq_len(N), each = T)
    dose <- rnorm(N)[id] + rnorm(N * T)            # varies between and within
    level <- rnorm(N, sd = 0.8)[id]                # stable subject level
    f <- level - 0.5 * (dose - ave(dose, id)) + rnorm(N * T, sd = 0.5)
    Y <- sapply(1:5, function(j) 0.8 * f + rnorm(N * T, sd = 0.5))
    colnames(Y) <- paste0("y", 1:5)
    wb <- within_between(data.frame(id = id, dose = dose), "dose", "id")
    sel <- as_selection(Y, paste0("y", 1:5), covariates = wb[, c("dose_w", "dose_b")])
    fit <- fit_sem(sel, cluster = id)
    expect_true(lavaan::lavInspect(fit$lavaan_fit, "converged"))
    expect_equal(fit$lavaan_fit@Data@cluster, ".cluster")
    lr <- lrt_structural(fit)
    w <- lr[lr$covariate == "dose_w", ]
    expect_lt(w$estimate, -0.3)
    expect_lt(w$p_lrt, 1e-4)
    expect_gt(lr$p_lrt[lr$covariate == "dose_b"], 0.01)
})

test_that("fit_riclpm recovers a cross-lagged path", {
    set.seed(15)
    N <- 400; T <- 3
    riA <- rnorm(N); riB <- rnorm(N)
    wA <- matrix(0, N, T); wB <- matrix(0, N, T)
    wA[, 1] <- rnorm(N); wB[, 1] <- rnorm(N)
    for (t in 2:T) {
        wA[, t] <- 0.3 * wA[, t - 1] + rnorm(N, sd = 0.8)
        wB[, t] <- 0.2 * wB[, t - 1] + 0.4 * wA[, t - 1] + rnorm(N, sd = 0.8)
    }
    long <- data.frame(id = rep(seq_len(N), each = T), visit = rep(seq_len(T), N),
                       A = as.vector(t(riA + wA)), B = as.vector(t(riB + wB)))
    rc <- fit_riclpm(long, id = "id", time = "visit", constructs = c("A", "B"))
    p <- rc$paths
    expect_equal(rc$subjects, N)
    ab <- p$estimate[p$from == "A" & p$to == "B"]
    ba <- p$estimate[p$from == "B" & p$to == "A"]
    expect_gt(ab, 0.2)
    expect_lt(abs(ba), 0.15)
})

test_that("test_invariance names the intercept that differs between groups", {
    set.seed(16)
    n <- 300; grp <- rep(c("a", "b"), each = n / 2)
    f <- rnorm(n)
    Y <- sapply(1:6, function(j) 0.8 * f + rnorm(n, sd = 0.6))
    Y[grp == "b", 3] <- Y[grp == "b", 3] + 1.0      # intercept shift in one feature
    colnames(Y) <- paste0("y", 1:6)
    sel <- as_selection(Y, paste0("y", 1:6), covariates = data.frame(grp = grp))
    inv <- test_invariance(sel, group_var = "grp")
    expect_equal(inv$worst$intercepts$variable[1], "y3")
})
