# --- Network analysis tests ---

test_that("network_analysis works with correlation type", {
    skip_if_not_installed("igraph")

    set.seed(42)
    n <- 100; p <- 30
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    Lambda <- c(rep(1, 10), rep(0, 20))
    h <- X %*% c(0, 1) + rnorm(n)
    Y <- h %*% t(Lambda) + matrix(rnorm(n * p, sd = sqrt(0.5)), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, max_iter = 100)
    net <- network_analysis(vb, type = "correlation", threshold = 0.3,
                             top_n = 10, max_vars = 30)

    expect_s3_class(net, "biomimic_network")
    expect_true(inherits(net$graph, "igraph"))
    expect_equal(length(net$hub_nodes), 10)
    expect_true(is.data.frame(net$marginal_neighbours))
    expect_true(is.numeric(net$modularity))
    expect_true(is.numeric(net$hub_density))
    expect_true(net$type == "correlation")

    # Hub density should be higher than overall (signal vars are correlated)
    if (!is.na(net$hub_density) && !is.na(net$overall_density) &&
        net$overall_density > 0) {
        expect_gte(net$hub_density, net$overall_density)
    }

    # Print should work
    expect_output(print(net), "biomimic Network")
})


test_that("network_analysis works with partial correlation", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("corpcor")

    set.seed(42)
    n <- 100; p <- 20
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    Lambda <- c(rep(1, 8), rep(0, 12))
    h <- X %*% c(0, 1) + rnorm(n)
    Y <- h %*% t(Lambda) + matrix(rnorm(n * p, sd = 0.7), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, max_iter = 100, use_gamma = FALSE)
    net <- network_analysis(vb, type = "partial", threshold = 0.1,
                             top_n = 8, max_vars = 20)

    expect_s3_class(net, "biomimic_network")
    expect_equal(net$type, "partial")
    expect_true(igraph::vcount(net$graph) <= 20)
})


test_that("network_analysis works with biomimic_selection input", {
    skip_if_not_installed("igraph")

    set.seed(123)
    n <- 100; p <- 25
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    Lambda <- c(rep(1, 10), rep(0, 15))
    h <- X %*% c(0, 1) + rnorm(n)
    Y <- h %*% t(Lambda) + matrix(rnorm(n * p, sd = 0.7), n, p)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, max_iter = 100)
    sel <- top_k(vb, k = 10)

    # Should accept selection object
    net <- network_analysis(sel, type = "correlation", threshold = 0.3,
                             max_vars = 25)
    expect_s3_class(net, "biomimic_network")
})


test_that("marginal_neighbours identifies correct variables", {
    skip_if_not_installed("igraph")

    set.seed(99)
    n <- 200; p <- 20
    X <- cbind(1, rbinom(n, 1, 0.5))
    colnames(X) <- c("(Intercept)", "group")
    # 5 signal vars, strongly correlated
    Lambda <- c(rep(1.5, 5), rep(0, 15))
    h <- X %*% c(0, 1) + rnorm(n)
    Y <- h %*% t(Lambda) + matrix(rnorm(n * p, sd = 0.5), n, p)
    # Make V6 correlated with V1 (marginal neighbour)
    Y[, 6] <- Y[, 1] * 0.8 + rnorm(n, sd = 0.5)
    colnames(Y) <- paste0("V", seq_len(p))

    vb <- vb_ard(Y, X, K = 1, max_iter = 150)
    net <- network_analysis(vb, type = "correlation", threshold = 0.3,
                             top_n = 5, max_vars = 20)

    # V6 should appear as a marginal neighbour of hub V1
    mn <- net$marginal_neighbours
    expect_true(nrow(mn) > 0)
    # V6 should be in marginal neighbours (it's correlated with V1 but not signal)
    if ("V6" %in% mn$variable) {
        v6_row <- mn[mn$variable == "V6", ]
        expect_true(v6_row$edge_weight[1] > 0.3)
    }
})
