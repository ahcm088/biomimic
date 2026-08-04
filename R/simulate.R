# --------------------------------------------------------------------------
# Data Generation for Monte Carlo Simulations
# Exported for reproducibility in papers and user experiments
# --------------------------------------------------------------------------

#' Simulate Data from a MIMIC Model
#'
#' Generates data from a known MIMIC model with optional direct effects,
#' multiple factors, and non-normal latent distributions. Designed for
#' Monte Carlo evaluation of the biomimic pipeline.
#'
#' @section Model:
#' The data-generating process is:
#' \deqn{h_i = B' x_i + \zeta_i, \quad \zeta_i \sim F}
#' \deqn{y_i = \Lambda h_i + \Gamma x_i + \varepsilon_i,
#'       \quad \varepsilon_i \sim N(0, \psi I)}
#'
#' where \eqn{F} is Normal, Student-t, or Contaminated Normal depending
#' on \code{gen_dist}. The first \code{p_signal} variables have nonzero
#' loadings (signal); the remaining \code{p_total - p_signal} are pure noise.
#'
#' @param n Integer, sample size.
#' @param p_total Integer, total number of observed variables.
#' @param p_signal Integer, number of signal variables (default 15).
#' @param K Integer, number of latent factors (default 1).
#' @param B_group Numeric or vector, group effect on latent factor(s).
#'   For K=1: scalar. For K=2: vector of length 2 (default c(1.0, 0.5)).
#'   Set to 0 for null hypothesis (no group effect).
#' @param lambda_strength Numeric, magnitude of factor loadings for signal
#'   variables (default 1.0).
#' @param gen_dist Character, distribution of latent disturbances:
#'   \code{"normal"}, \code{"t"} (Student-t), or \code{"contaminated"}
#'   (contaminated normal). Default \code{"normal"}.
#' @param nu Numeric, degrees of freedom for Student-t (default 3).
#' @param xi Numeric, contamination proportion for contaminated normal
#'   (default 0.1).
#' @param contam_scale Numeric, scale inflation for contaminated component
#'   (default 3).
#' @param has_direct Logical, whether to include direct effects (default FALSE).
#' @param n_direct Integer, number of variables with direct effects (default 3).
#'   Direct effects are placed on the first \code{n_direct} noise variables
#'   (i.e., variables just outside the signal set).
#' @param gamma_strength Numeric, magnitude of direct effects (default 1.5).
#' @param psi Numeric, idiosyncratic variance (default 0.5).
#' @param group_prevalence Numeric, probability of group=1 (default 0.5).
#' @param seed Integer or NULL, random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y}{Numeric matrix (n x p_total)}
#'   \item{X}{Numeric matrix (n x 2) with intercept and group indicator}
#'   \item{Lambda_true}{Numeric matrix (p_total x K) of true loadings}
#'   \item{B_true}{Numeric matrix (2 x K) of true structural coefficients}
#'   \item{Gamma_true}{Numeric matrix (p_total x 2) of true direct effects}
#'   \item{Psi_true}{Numeric vector (p_total) of true idiosyncratic variances}
#'   \item{h_true}{Numeric matrix (n x K) of true latent scores}
#'   \item{signal_vars}{Character vector of signal variable names}
#'   \item{noise_vars}{Character vector of noise variable names}
#'   \item{direct_vars}{Character vector of variables with direct effects}
#'   \item{indicator_groups}{List of integer vectors (for K>=2)}
#'   \item{n}{Sample size}
#'   \item{p_total}{Total variables}
#'   \item{p_signal}{Signal variables}
#'   \item{K}{Number of factors}
#'   \item{gen_dist}{Distribution used}
#' }
#'
#' @examples
#' d <- simulate_mimic(n = 200, p_total = 50, K = 1, B_group = 1.0)
#' dim(d$Y)  # 200 x 50
#' d$signal_vars  # first 15 variables
#'
#' @export
simulate_mimic <- function(n = 200L, p_total = 50L, p_signal = 15L,
                            K = 1L,
                            B_group = 1.0,
                            lambda_strength = 1.0,
                            gen_dist = c("normal", "t", "contaminated"),
                            nu = 3, xi = 0.1, contam_scale = 3,
                            has_direct = FALSE,
                            n_direct = 3L,
                            gamma_strength = 1.5,
                            psi = 0.5,
                            group_prevalence = 0.5,
                            seed = NULL) {

    gen_dist <- match.arg(gen_dist)
    if (!is.null(seed)) set.seed(seed)

    stopifnot(p_signal <= p_total)
    stopifnot(K >= 1L)

    # --- Covariates: intercept + binary group ---
    group <- stats::rbinom(n, 1, group_prevalence)
    X <- cbind(1, group)
    colnames(X) <- c("(Intercept)", "group")
    q <- 2L

    # --- Structural coefficients B (q x K) ---
    B_true <- matrix(0, nrow = q, ncol = K)
    # Intercept = 0 for all factors
    if (K == 1L) {
        B_true[2, 1] <- B_group[1]
    } else {
        B_vec <- rep_len(B_group, K)
        B_true[2, ] <- B_vec
    }

    # --- Latent disturbances ---
    zeta <- .generate_latent(n, K, gen_dist, nu, xi, contam_scale)

    # --- Latent scores ---
    h_true <- X %*% B_true + zeta  # n x K

    # --- Loading matrix Lambda (p_total x K) ---
    Lambda_true <- matrix(0, nrow = p_total, ncol = K)

    if (K == 1L) {
        Lambda_true[seq_len(p_signal), 1] <- lambda_strength
    } else {
        # Simple structure: divide signal vars into K blocks
        block_size <- p_signal %/% K
        remainder <- p_signal %% K
        start <- 1L
        indicator_groups <- vector("list", K)
        for (k in seq_len(K)) {
            sz <- block_size + (if (k <= remainder) 1L else 0L)
            idx <- seq(start, start + sz - 1L)
            Lambda_true[idx, k] <- lambda_strength
            indicator_groups[[k]] <- idx
            start <- start + sz
        }
    }

    # --- Direct effects Gamma (p_total x q) ---
    Gamma_true <- matrix(0, nrow = p_total, ncol = q)
    direct_var_idx <- integer(0)

    if (has_direct && n_direct > 0) {
        # Place direct effects on the first n_direct SIGNAL proteins, so they
        # carry BOTH a loading and a direct covariate effect (partial
        # mediation), matching the main-text description ("direct effects for
        # signal proteins").
        direct_var_idx <- seq_len(min(n_direct, p_signal))
        Gamma_true[direct_var_idx, 2] <- gamma_strength  # group direct effect
    }

    # --- Observed data ---
    Psi_true <- rep(psi, p_total)
    E <- matrix(stats::rnorm(n * p_total, sd = sqrt(psi)), n, p_total)
    Y <- h_true %*% t(Lambda_true) + X %*% t(Gamma_true) + E

    # --- Variable names ---
    var_names <- paste0("V", seq_len(p_total))
    colnames(Y) <- var_names
    signal_vars <- var_names[seq_len(p_signal)]
    noise_vars <- var_names[(p_signal + 1L):p_total]
    direct_vars <- if (length(direct_var_idx) > 0) {
        var_names[direct_var_idx]
    } else {
        character(0)
    }

    # Indicator groups for K >= 2
    if (K == 1L) {
        indicator_groups <- list(seq_len(p_signal))
    }

    list(
        Y                = Y,
        X                = X,
        Lambda_true      = Lambda_true,
        B_true           = B_true,
        Gamma_true       = Gamma_true,
        Psi_true         = Psi_true,
        h_true           = h_true,
        signal_vars      = signal_vars,
        noise_vars       = noise_vars,
        direct_vars      = direct_vars,
        indicator_groups = indicator_groups,
        n                = n,
        p_total          = p_total,
        p_signal         = p_signal,
        K                = K,
        gen_dist         = gen_dist
    )
}


#' @noRd
.generate_latent <- function(n, K, gen_dist, nu, xi, contam_scale) {
    # The structural disturbance is N(0, I_K) under the model, so the
    # non-normal generators are standardised to UNIT variance: this isolates
    # the tail shape from a change of scale, keeping cross-distribution
    # bias/RMSE/coverage comparisons on the same footing.
    switch(gen_dist,
        normal = {
            matrix(stats::rnorm(n * K), n, K)
        },
        t = {
            sd_t <- if (nu > 2) sqrt(nu / (nu - 2)) else 1  # Var(t_nu)=nu/(nu-2)
            matrix(stats::rt(n * K, df = nu), n, K) / sd_t
        },
        contaminated = {
            z <- matrix(stats::rnorm(n * K), n, K)
            contam <- stats::rbinom(n, 1, xi)
            # Scale contaminated observations
            for (i in seq_len(n)) {
                if (contam[i] == 1L) {
                    z[i, ] <- z[i, ] * contam_scale
                }
            }
            # Mixture variance: (1-xi)*1 + xi*contam_scale^2  -> standardise
            z / sqrt((1 - xi) * 1 + xi * contam_scale^2)
        }
    )
}


#' Compute Simulation Metrics for Parameter Recovery
#'
#' Given true and estimated parameter values across replications, computes
#' bias, RMSE, and empirical coverage of confidence intervals.
#'
#' @param estimates Numeric vector of point estimates across replications.
#' @param true_value Numeric scalar, the true parameter value.
#' @param ci_lower Numeric vector of CI lower bounds (same length as estimates).
#' @param ci_upper Numeric vector of CI upper bounds.
#'
#' @return A named list with: bias, rmse, coverage, mean_estimate, sd_estimate.
#'
#' @export
sim_metrics <- function(estimates, true_value, ci_lower = NULL,
                         ci_upper = NULL) {
    valid <- !is.na(estimates)
    est <- estimates[valid]

    bias <- mean(est) - true_value
    rmse <- sqrt(mean((est - true_value)^2))

    coverage <- NA_real_
    if (!is.null(ci_lower) && !is.null(ci_upper)) {
        lo <- ci_lower[valid]
        hi <- ci_upper[valid]
        coverage <- mean(lo <= true_value & true_value <= hi)
    }

    list(
        bias          = bias,
        rmse          = rmse,
        coverage      = coverage,
        mean_estimate = mean(est),
        sd_estimate   = stats::sd(est),
        n_valid       = sum(valid),
        n_total       = length(estimates)
    )
}


#' Compute Screening Metrics
#'
#' Evaluates how well the VB-ARD SNR ranking recovers the true signal
#' variables.
#'
#' @param ranking A data.frame from \code{rank_variables()}.
#' @param true_signal Character vector of true signal variable names.
#' @param k Integer, selection cutoff (default 15).
#'
#' @return A named list with: auc, f1, sensitivity, specificity,
#'   precision, n_selected.
#'
#' @export
sim_screening_metrics <- function(ranking, true_signal, k = 15L) {
    selected <- ranking$variable[seq_len(min(k, nrow(ranking)))]

    tp <- sum(selected %in% true_signal)
    fp <- sum(!selected %in% true_signal)
    fn <- sum(!true_signal %in% selected)
    # TP+FP+TN+FN = N (all ranked variables); selected = TP+FP = length(selected)
    tn <- nrow(ranking) - length(selected) - fn

    sensitivity <- tp / max(tp + fn, 1)
    specificity <- tn / max(tn + fp, 1)
    precision   <- tp / max(tp + fp, 1)
    f1 <- if (precision + sensitivity > 0) {
        2 * precision * sensitivity / (precision + sensitivity)
    } else 0

    # Ranking quality: AUC (= normalised Mann-Whitney U), the probability that
    # a signal variable outranks a non-signal one by SNR. This is the correct
    # measure of how well a CONTINUOUS score separates a BINARY signal set
    # (a Spearman correlation between continuous SNR and a binary indicator is
    # hard-bounded well below 1 by the ties and is NOT a meaningful rank metric).
    is_signal <- as.integer(ranking$variable %in% true_signal)
    pos <- ranking$snr[is_signal == 1L]
    neg <- ranking$snr[is_signal == 0L]
    auc <- if (length(pos) > 0L && length(neg) > 0L) {
        r <- rank(c(pos, neg))
        (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
            (length(pos) * length(neg))
    } else NA_real_

    list(
        auc          = auc,
        f1           = f1,
        sensitivity  = sensitivity,
        specificity  = specificity,
        precision    = precision,
        tp           = tp,
        fp           = fp,
        tn           = tn,
        n_selected   = length(selected)
    )
}


#' Run a factorial MIMIC pipeline simulation study
#'
#' Orchestrates the end-to-end pipeline simulation: for every scenario in
#' \code{design} and every replication it generates data with
#' [simulate_mimic()], runs the full [biomimic()] pipeline, and scores
#' screening recovery ([sim_screening_metrics()]) and the coverage of the
#' structural group effect. Replications are independent, deterministically
#' self-seeded, and optionally parallelised. This is the function behind the
#' paper's pipeline simulation block.
#'
#' @param design A data frame with one row per scenario. Recognised columns
#'   (any omitted use the defaults): \code{n}, \code{p_total}, \code{K},
#'   \code{gen_dist}.
#' @param n_reps Replications per scenario.
#' @param B_true True structural group effect used to generate data and to
#'   assess interval coverage (default 1).
#' @param k Panel size selected per replication (default 15).
#' @param parallel Logical; distribute replications over a PSOCK cluster.
#' @param n_cores Number of workers (default: \code{BIOMIMIC_NCORES} or
#'   \code{detectCores() - 1}).
#' @param seed_base Integer offset added to every per-replication seed.
#' @param max_iter Passed to [vb_ard()].
#'
#' @return A list with:
#'   \describe{
#'     \item{metrics}{Data frame, one row per scenario \eqn{\times} replication:
#'       scenario columns, \code{rep}, \code{model_type}, \code{converged},
#'       \code{screening_f1} and \code{B_covered}.}
#'     \item{summary}{Data frame, one row per scenario: scenario columns,
#'       \code{n_valid}, \code{convergence_rate}, \code{mean_screening_f1} and
#'       \code{coverage}.}
#'   }
#'
#' @examples
#' \donttest{
#' design <- expand.grid(n = c(100, 200), p_total = 50, K = 1)
#' sim <- run_simulation(design, n_reps = 5, parallel = FALSE)
#' sim$summary
#' }
#' @seealso [simulate_mimic()], [sim_screening_metrics()], [biomimic()]
#' @export
run_simulation <- function(design, n_reps = 100L, B_true = 1.0, k = 15L,
                           parallel = TRUE, n_cores = NULL, seed_base = 50000L,
                           max_iter = 200L) {
    design <- as.data.frame(design, stringsAsFactors = FALSE)
    defaults <- list(n = 200L, p_total = 50L, K = 1L, gen_dist = "normal")
    for (nm in names(defaults))
        if (is.null(design[[nm]])) design[[nm]] <- defaults[[nm]]
    n_scn <- nrow(design)
    n_reps <- as.integer(n_reps)
    jobs <- expand.grid(rep = seq_len(n_reps), scenario = seq_len(n_scn))

    run_one <- function(j) {
        s <- jobs$scenario[j]; r <- jobs$rep[j]; d <- design[s, , drop = FALSE]
        mt <- if (d$K == 1L) "mimic" else "multi_mimic"
        seed <- seed_base + s * 1000L + r
        bg <- if (d$K == 1L) B_true else c(B_true, 0.5)
        sim <- simulate_mimic(n = d$n, p_total = d$p_total, p_signal = 15L,
                              K = d$K, B_group = bg, gen_dist = d$gen_dist,
                              seed = seed)
        base <- data.frame(scenario = s, d, rep = r, model_type = mt,
                           row.names = NULL, stringsAsFactors = FALSE)
        res <- tryCatch(biomimic(sim$Y, sim$X, K = d$K, k = k, model_type = mt,
                                 max_iter = max_iter, use_ard = (d$K == 1L)),
                        error = function(e) NULL)
        if (is.null(res))
            return(cbind(base, converged = FALSE, screening_f1 = NA_real_,
                         B_covered = NA))
        sm <- sim_screening_metrics(
            rank_variables(res$selection$fit, type = "lambda"),
            sim$signal_vars, k = k)
        conv <- isTRUE(lavaan::lavInspect(res$lavaan_fit, "converged"))
        pe <- lavaan::parameterEstimates(res$lavaan_fit, ci = TRUE)
        st <- pe[pe$op == "~" & grepl("^F", pe$lhs) & pe$rhs == "group", ]
        b_cov <- NA
        if (nrow(st) > 0L && conv)
            b_cov <- st$ci.lower[1] <= B_true && B_true <= st$ci.upper[1]
        cbind(base, converged = conv, screening_f1 = sm$f1, B_covered = b_cov)
    }

    idx <- seq_len(nrow(jobs))
    if (is.null(n_cores)) {
        env <- suppressWarnings(as.integer(Sys.getenv("BIOMIMIC_NCORES", "")))
        n_cores <- if (!is.na(env) && env >= 1L) env
                   else max(1L, tryCatch(parallel::detectCores() - 1L,
                                         error = function(e) 1L))
    }
    use_par <- isTRUE(parallel) && n_cores > 1L &&
        requireNamespace("parallel", quietly = TRUE)
    if (use_par) {
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterExport(cl,
            c("jobs", "design", "B_true", "k", "max_iter", "seed_base"),
            envir = environment())
        parallel::clusterEvalQ(cl, suppressMessages(library(biomimic)))
        res <- parallel::parLapply(cl, idx, run_one)
    } else {
        res <- lapply(idx, run_one)
    }
    metrics <- do.call(rbind, res)

    summary <- do.call(rbind, lapply(seq_len(n_scn), function(s) {
        sub <- metrics[metrics$scenario == s, , drop = FALSE]
        data.frame(scenario = s, design[s, , drop = FALSE],
                   n_valid = sum(sub$converged, na.rm = TRUE),
                   convergence_rate = mean(sub$converged, na.rm = TRUE),
                   mean_screening_f1 = mean(sub$screening_f1, na.rm = TRUE),
                   coverage = mean(sub$B_covered, na.rm = TRUE),
                   row.names = NULL, stringsAsFactors = FALSE)
    }))
    list(metrics = metrics, summary = summary)
}
