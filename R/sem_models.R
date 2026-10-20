# --------------------------------------------------------------------------
# Structural models beyond the single-factor MIMIC
# Bifactor models, mediation, per-feature decomposition, random-intercept
# cross-lagged panels and the within/between split, each written for panels
# of omics features and fitted in lavaan.
# --------------------------------------------------------------------------


#' Bifactor Model: Specific Responses Beyond a General Factor
#'
#' A general factor loads on every indicator and orthogonal specific factors
#' load on named panels; covariates act on all factors (MIMIC). In
#' autoantibody arrays the general factor absorbs each serum's global level
#' of reactivity, and the specific factors test whether a known response
#' (for example Ro60-La in SLE) is shifted beyond it.
#'
#' @param Y Numeric matrix, samples x features, with column names.
#' @param specific Named list of character vectors, one per specific factor.
#' @param general Character vector of further features that load only on the
#'   general factor, or a single integer: that many features drawn at random
#'   from the rest of \code{Y} (with \code{seed}). Default 30.
#' @param covariates Data frame, one row per sample (e.g. group, ethnicity).
#' @param equal_loadings Character vector of specific factors whose loadings
#'   are constrained equal; by default those with two indicators, which are
#'   otherwise not identified.
#' @param seed Integer seed for the random general indicators.
#' @param estimator lavaan estimator (default \code{"MLR"}).
#' @param lrt Logical, also test each factor's covariate paths by likelihood
#'   ratio (default \code{TRUE}).
#'
#' @return A list of class \code{biomimic_bifactor}: \code{lavaan_fit},
#'   \code{effects} (factor, covariate, estimate, se, Wald p, CI, LRT p),
#'   \code{loadings} (mean standardised loading of each specific factor and
#'   of the general factor on its targets), \code{ecv} (explained common
#'   variance of the general factor on each panel), \code{syntax}.
#'
#' @details With few samples, a specific factor with weakly correlated
#'   indicators may not be empirically identified and the fit may fail; in
#'   the thesis this happened with a two-cytokine Th17 factor in 72 sera.
#'
#' @export
fit_bifactor <- function(Y, specific, general = 30L, covariates = NULL,
                         equal_loadings = NULL, seed = 1L,
                         estimator = "MLR", lrt = TRUE) {
    stopifnot(is.list(specific), !is.null(names(specific)))
    Y <- as.matrix(Y)
    specific <- lapply(specific, function(p) intersect(p, colnames(Y)))
    if (any(lengths(specific) < 2L)) stop("Each specific factor needs >= 2 features in Y.", call. = FALSE)
    targets <- unlist(specific, use.names = FALSE)
    if (is.numeric(general) && length(general) == 1L) {
        pool <- setdiff(colnames(Y), targets)
        general <- .with_seed(seed, sample(pool, min(general, length(pool))))
    }
    general <- setdiff(intersect(general, colnames(Y)), targets)
    if (is.null(equal_loadings)) equal_loadings <- names(specific)[lengths(specific) == 2L]
    allv <- c(targets, general)
    safe <- stats::setNames(make.names(allv, unique = TRUE), allv)
    D <- stats::setNames(as.data.frame(scale(Y[, allv, drop = FALSE])), safe)
    cv <- character(0)
    if (!is.null(covariates)) { covariates <- as.data.frame(covariates); cv <- names(covariates); D <- cbind(D, covariates) }
    fs <- c("G", names(specific))
    lines <- paste("G =~", paste(safe, collapse = " + "))
    for (s in names(specific)) {
        v <- safe[specific[[s]]]
        lines <- c(lines, paste0(s, " =~ ", paste0(if (s %in% equal_loadings) paste0("l_", s, "*") else "",
                                                   v, collapse = " + ")))
    }
    lines <- c(lines, utils::combn(fs, 2, function(x) paste(x[1], "~~ 0*", x[2])))
    path <- function(f, cvs) if (length(cvs)) paste(f, "~", paste(paste0("b_", f, "_", cvs, "*", cvs), collapse = " + "))
    syntax <- paste(c(lines, unlist(lapply(fs, path, cvs = cv))), collapse = "\n")
    fit <- tryCatch(lavaan::sem(syntax, data = D, estimator = estimator, std.lv = TRUE),
                    error = function(e) NULL)
    if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) {
        warning("fit_bifactor: the model did not converge.", call. = FALSE)
        return(structure(list(lavaan_fit = fit, effects = NULL, syntax = syntax,
                              converged = FALSE), class = "biomimic_bifactor"))
    }
    pe <- lavaan::parameterEstimates(fit, ci = TRUE)
    ss <- lavaan::standardizedSolution(fit)
    eff <- pe[pe$op == "~" & pe$lhs %in% fs, c("lhs", "rhs", "est", "se", "pvalue", "ci.lower", "ci.upper")]
    names(eff) <- c("factor", "covariate", "estimate", "se", "p", "lo", "hi")
    eff$p_lrt <- NA_real_
    if (isTRUE(lrt) && nrow(eff)) {
        for (i in seq_len(nrow(eff))) {
            lab <- paste0("b_", eff$factor[i], "_", eff$covariate[i])
            reduced <- sub(paste0(lab, "\\*"), "0*", syntax)
            fr <- tryCatch(lavaan::sem(reduced, data = D, estimator = estimator, std.lv = TRUE),
                           error = function(e) NULL)
            if (!is.null(fr) && isTRUE(lavaan::lavInspect(fr, "converged")))
                eff$p_lrt[i] <- tryCatch(lavaan::lavTestLRT(fit, fr)[["Pr(>Chisq)"]][2],
                                         error = function(e) NA_real_)
        }
    }
    ld <- ss[ss$op == "=~", ]
    loadings <- data.frame(
        factor = names(specific),
        specific_loading = vapply(names(specific), function(s) mean(ld$est.std[ld$lhs == s]), numeric(1)),
        general_loading_on_targets = vapply(names(specific), function(s)
            mean(ld$est.std[ld$lhs == "G" & ld$rhs %in% safe[specific[[s]]]]), numeric(1)),
        row.names = NULL)
    ecv <- vapply(names(specific), function(s) {
        g2 <- sum(ld$est.std[ld$lhs == "G" & ld$rhs %in% safe[specific[[s]]]]^2)
        s2 <- sum(ld$est.std[ld$lhs == s]^2)
        g2 / (g2 + s2)
    }, numeric(1))
    rownames(eff) <- NULL
    structure(list(lavaan_fit = fit, effects = eff, loadings = loadings, ecv = ecv,
                   general = general, n_indicators = length(allv), syntax = syntax,
                   converged = TRUE),
              class = "biomimic_bifactor")
}

#' @export
print.biomimic_bifactor <- function(x, ...) {
    cat("biomimic bifactor model\n")
    if (!isTRUE(x$converged)) { cat("  did not converge\n"); return(invisible(x)) }
    cat(sprintf("  general factor on %d features; specific: %s\n",
                x$n_indicators, paste(x$loadings$factor, collapse = ", ")))
    e <- x$effects
    for (i in seq_len(nrow(e))) cat(sprintf("  %-10s ~ %-12s %+.3f [%+.3f, %+.3f]  p = %.2g%s\n",
        e$factor[i], e$covariate[i], e$estimate[i], e$lo[i], e$hi[i], e$p[i],
        if (is.finite(e$p_lrt[i])) sprintf(" (LRT %.2g)", e$p_lrt[i]) else ""))
    cat("  share of each panel's common variance on the general factor:",
        paste(sprintf("%s %.2f", names(x$ecv), x$ecv), collapse = ", "), "\n")
    invisible(x)
}


#' Mediation with Latent or Observed Mediators
#'
#' Fits exposure -> mediators -> outcome with a direct exposure -> outcome
#' path; mediator residuals are correlated. Mediators given as feature
#' panels are latent factors measured by those features; mediators given as
#' a single column name are observed. Each mediator's contribution is tested
#' by likelihood ratio (its mediator -> outcome path fixed to zero).
#'
#' @param data Data frame with the exposure, the outcome, observed
#'   mediators, covariates and the panel features.
#' @param exposure,outcome Character, column names.
#' @param mediators Named list: each element a character vector of panel
#'   features (latent mediator) or a single column name (observed).
#' @param covariates Optional character vector of columns acting on every
#'   endogenous variable.
#' @param estimator lavaan estimator (default \code{"MLR"}).
#'
#' @return A list of class \code{biomimic_mediation}: \code{lavaan_fit},
#'   \code{effects} (each mediator's a, b, indirect effect with CI and p,
#'   LRT p for its b path, and share of the total), \code{direct},
#'   \code{total}.
#'
#' @export
fit_mediation <- function(data, exposure, mediators, outcome, covariates = NULL,
                          estimator = "MLR") {
    stopifnot(is.list(mediators), !is.null(names(mediators)))
    data <- as.data.frame(data)
    cv <- if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""
    mname <- character(0); lines <- character(0)
    for (m in names(mediators)) {
        v <- mediators[[m]]
        if (length(v) == 1L && v %in% names(data)) {
            mname[m] <- v
        } else {
            v <- intersect(v, names(data))
            if (length(v) < 2L) stop("Latent mediator ", m, " needs >= 2 features.", call. = FALSE)
            nm <- make.names(m)
            lines <- c(lines, paste(nm, "=~", paste(v, collapse = " + ")))
            mname[m] <- nm
        }
    }
    idx <- seq_along(mname)
    build <- function(drop = NULL) {
        b <- paste0("b", idx); b[names(mname) %in% drop] <- "0"
        L <- c(lines,
               sprintf("%s ~ a%d*%s%s", mname, idx, exposure, cv),
               sprintf("%s ~ %s + cp*%s%s", outcome, paste(paste0(b, "*", mname), collapse = " + "),
                       exposure, cv))
        if (length(mname) >= 2L) L <- c(L, utils::combn(unname(mname), 2, function(x) paste(x[1], "~~", x[2])))
        if (is.null(drop)) L <- c(L, sprintf("ind%d := a%d*b%d", idx, idx, idx),
                                  paste("total :=", paste(sprintf("a%d*b%d", idx, idx), collapse = " + "), "+ cp"))
        paste(L, collapse = "\n")
    }
    fit <- lavaan::sem(build(), data = data, estimator = estimator)
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) warning("fit_mediation: did not converge.", call. = FALSE)
    pe <- lavaan::parameterEstimates(fit, ci = TRUE)
    g <- function(l) pe[pe$label == l, ][1, ]
    tot <- g("total")$est
    effects <- do.call(rbind, lapply(idx, function(i) {
        fr <- tryCatch(lavaan::sem(build(drop = names(mname)[i]), data = data, estimator = estimator),
                       error = function(e) NULL)
        plr <- if (is.null(fr)) NA_real_ else
            tryCatch(lavaan::lavTestLRT(fit, fr)[["Pr(>Chisq)"]][2], error = function(e) NA_real_)
        ind <- g(paste0("ind", i))
        data.frame(mediator = names(mname)[i], a = g(paste0("a", i))$est, b = g(paste0("b", i))$est,
                   indirect = ind$est, lo = ind$ci.lower, hi = ind$ci.upper, p = ind$pvalue,
                   p_lrt = plr, share = ind$est / tot, stringsAsFactors = FALSE)
    }))
    dct <- g("cp")
    structure(list(lavaan_fit = fit, effects = effects,
                   direct = c(estimate = dct$est, lo = dct$ci.lower, hi = dct$ci.upper, p = dct$pvalue),
                   total = c(estimate = tot, lo = g("total")$ci.lower, hi = g("total")$ci.upper,
                             p = g("total")$pvalue),
                   syntax = build()), class = "biomimic_mediation")
}

#' @export
print.biomimic_mediation <- function(x, ...) {
    cat("biomimic mediation\n")
    e <- x$effects
    for (i in seq_len(nrow(e))) cat(sprintf("  via %-14s %+.3f [%+.3f, %+.3f]  p = %.2g  LRT %.2g  (%.0f%% of total)\n",
        e$mediator[i], e$indirect[i], e$lo[i], e$hi[i], e$p[i], e$p_lrt[i], 100 * e$share[i]))
    cat(sprintf("  direct           %+.3f, p = %.2g; total %+.3f\n", x$direct[["estimate"]],
                x$direct[["p"]], x$total[["estimate"]]))
    invisible(x)
}


#' Decompose Each Feature's Group Effect Through Named Factors
#'
#' For every target feature, fits the factors measured by their panels (all
#' regressed on the group, residuals correlated) and the target regressed on
#' the factors and on the group. The group effect on the target splits into
#' a part through each factor and the target's own part, with delta-method
#' standard errors.
#'
#' @param Y Numeric matrix, samples x features (standardise first).
#' @param panels Named list of character vectors, one per factor.
#' @param targets Character vector of features to decompose (panel features
#'   are skipped).
#' @param group Numeric vector, the group indicator.
#' @param estimator lavaan estimator (default \code{"MLR"}).
#' @param n_cores Integer, parallel workers (default all but one; honours
#'   \code{BIOMIMIC_NCORES}).
#'
#' @return A data frame, one row per target: \code{total}, for each factor
#'   \code{via_<name>} and \code{via_<name>_se}, \code{own}, \code{own_se},
#'   \code{converged}.
#'
#' @details In SLE blood the part through the interferon factor replicated
#'   across five cohorts and the own part did not (thesis, Chapter 2).
#'   With coherent panels the result is close to least-squares mediation on
#'   panel means; the latent version protects against attributing a factor's
#'   effect to the target when panels are small or weakly coherent.
#'
#' @export
decompose_effects <- function(Y, panels, targets, group, estimator = "MLR",
                              n_cores = NULL) {
    stopifnot(is.list(panels), !is.null(names(panels)))
    Y <- as.matrix(Y)
    panels <- lapply(panels, function(p) intersect(p, colnames(Y)))
    targets <- setdiff(intersect(targets, colnames(Y)), unlist(panels))
    pv <- unlist(panels, use.names = FALSE)
    fn <- make.names(names(panels))
    one <- function(task) {
        nmf <- function(x) make.names(x)
        df <- as.data.frame(task$Y); colnames(df) <- nmf(colnames(df)); df$group <- task$g
        L <- c(vapply(seq_along(task$panels), function(k)
                   paste0(task$fn[k], " =~ ", paste(nmf(task$panels[[k]]), collapse = " + ")), ""),
               sprintf("%s ~ b%d*group", task$fn, seq_along(task$fn)))
        if (length(task$fn) >= 2L) L <- c(L, utils::combn(task$fn, 2, function(x) paste(x[1], "~~", x[2])))
        tg <- nmf(task$gene)
        L <- c(L, paste0(tg, " ~ ", paste(sprintf("l%d*%s", seq_along(task$fn), task$fn), collapse = " + "),
                         " + gd*group"),
               sprintf("via%d := l%d*b%d", seq_along(task$fn), seq_along(task$fn), seq_along(task$fn)),
               paste("total :=", paste(sprintf("l%d*b%d", seq_along(task$fn), seq_along(task$fn)),
                                       collapse = " + "), "+ gd"))
        f <- tryCatch(lavaan::sem(paste(L, collapse = "\n"), data = df, estimator = task$est),
                      error = function(e) NULL)
        k <- length(task$fn)
        if (is.null(f) || !isTRUE(lavaan::lavInspect(f, "converged")))
            return(c(total = NA, stats::setNames(rep(NA, 2 * k), c(rbind(paste0("via", 1:k), paste0("via", 1:k, "_se")))),
                     own = NA, own_se = NA, converged = 0))
        pe <- lavaan::parameterEstimates(f)
        gt <- function(l) pe[pe$label == l, c("est", "se")][1, ]
        via <- unlist(lapply(seq_len(k), function(i) { x <- gt(paste0("via", i)); c(x$est, x$se) }))
        names(via) <- c(rbind(paste0("via", 1:k), paste0("via", 1:k, "_se")))
        c(total = gt("total")$est, via, own = gt("gd")$est, own_se = gt("gd")$se, converged = 1)
    }
    tasks <- lapply(targets, function(j) list(Y = Y[, c(pv, j), drop = FALSE], g = group,
                                              panels = panels, fn = fn, gene = j, est = estimator))
    n_cores <- n_cores %||% .biomimic_n_cores()
    res <- if (n_cores > 1L && length(tasks) > 1L) {
        cl <- parallel::makePSOCKcluster(n_cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::parLapply(cl, tasks, one)
    } else lapply(tasks, one)
    out <- as.data.frame(do.call(rbind, res))
    nm <- names(out)
    for (i in seq_along(fn)) nm <- sub(paste0("^via", i, "(_se)?$"), paste0("via_", names(panels)[i], "\\1"), nm)
    names(out) <- nm
    cbind(feature = targets, out, stringsAsFactors = FALSE)
}


#' Within- and Between-Subject Parts of Repeated Measures
#'
#' Splits each variable into the subject's mean (between) and the deviation
#' from it at each visit (within), so that a model can estimate the effect of
#' a change within the same person apart from differences between people.
#'
#' @param data Data frame in long format (one row per visit).
#' @param vars Character vector of variables to split.
#' @param id Character, the subject identifier column.
#'
#' @return \code{data} with columns \code{<var>_w} and \code{<var>_b} added.
#'
#' @export
within_between <- function(data, vars, id) {
    data <- as.data.frame(data)
    for (v in vars) {
        m <- stats::ave(data[[v]], data[[id]], FUN = function(x) mean(x, na.rm = TRUE))
        data[[paste0(v, "_w")]] <- data[[v]] - m
        data[[paste0(v, "_b")]] <- m
    }
    data
}


#' Random-Intercept Cross-Lagged Panel Model
#'
#' For repeated measures of several constructs, separates stable differences
#' between subjects (random intercepts) from within-subject dynamics, and
#' estimates autoregressive and cross-lagged paths between consecutive
#' waves (Hamaker, Kuiper and Grasman 2015). Paths are constrained equal
#' across lags.
#'
#' @param data Data frame in long format.
#' @param id,time Character, subject and visit-order columns.
#' @param constructs Character vector of observed construct columns (for
#'   example panel scores and a serology measure).
#' @param waves Integer, number of waves used per subject: the first
#'   \code{waves} visits of subjects with at least that many (default 3).
#' @param estimator lavaan estimator (default \code{"MLR"}).
#'
#' @return A list of class \code{biomimic_riclpm}: \code{lavaan_fit},
#'   \code{paths} (from, to, estimate, se, p, type), \code{subjects},
#'   \code{within_share} (share of each construct's variance within
#'   subjects, over all visits).
#'
#' @export
fit_riclpm <- function(data, id, time, constructs, waves = 3L, estimator = "MLR") {
    data <- as.data.frame(data)
    stopifnot(all(c(id, time, constructs) %in% names(data)), waves >= 3L)
    o <- order(data[[id]], data[[time]])
    d <- data[o, ]
    d$.wave <- stats::ave(seq_len(nrow(d)), d[[id]], FUN = seq_along)
    keep <- names(which(table(d[[id]][d$.wave <= waves]) == waves))
    d <- d[d[[id]] %in% keep & d$.wave <= waves, ]
    cn <- make.names(constructs)
    W <- data.frame(id = keep)
    for (t in seq_len(waves)) {
        dt <- d[d$.wave == t, ]; dt <- dt[match(keep, dt[[id]]), ]
        for (k in seq_along(constructs)) W[[paste0(cn[k], "_", t)]] <- dt[[constructs[k]]]
    }
    L <- character(0)
    for (k in seq_along(cn)) {
        L <- c(L, sprintf("RI_%s =~ %s", cn[k], paste(sprintf("1*%s_%d", cn[k], seq_len(waves)), collapse = " + ")))
        for (t in seq_len(waves)) L <- c(L, sprintf("w%s_%d =~ 1*%s_%d", cn[k], t, cn[k], t),
                                         sprintf("%s_%d ~~ 0*%s_%d", cn[k], t, cn[k], t))
    }
    if (length(cn) >= 2L) L <- c(L, utils::combn(paste0("RI_", cn), 2, function(x) paste(x[1], "~~", x[2])))
    for (r in paste0("RI_", cn)) for (k in cn) L <- c(L, sprintf("%s ~~ 0*w%s_1", r, k))
    lab <- function(from, to) paste0("p_", from, "_", to)
    for (t in 2:waves) for (k in cn) {
        L <- c(L, sprintf("w%s_%d ~ %s", k, t, paste(sprintf("%s*w%s_%d", lab(cn, k), cn, t - 1), collapse = " + ")))
    }
    for (t in seq_len(waves)) if (length(cn) >= 2L)
        L <- c(L, utils::combn(sprintf("w%s_%d", cn, t), 2, function(x) paste(x[1], "~~", x[2])))
    fit <- lavaan::lavaan(paste(L, collapse = "\n"), data = W, estimator = estimator, missing = "fiml",
                          meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE,
                          auto.var = TRUE, auto.fix.first = FALSE, auto.cov.lv.x = FALSE)
    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) warning("fit_riclpm: did not converge.", call. = FALSE)
    pe <- lavaan::parameterEstimates(fit)
    paths <- do.call(rbind, lapply(cn, function(f) do.call(rbind, lapply(cn, function(t) {
        x <- pe[pe$label == lab(f, t), ][1, ]
        data.frame(from = constructs[cn == f], to = constructs[cn == t], estimate = x$est, se = x$se,
                   p = x$pvalue, type = if (f == t) "autoregressive" else "cross-lagged",
                   stringsAsFactors = FALSE)
    }))))
    ws <- vapply(constructs, function(v) {
        x <- data[[v]]; m <- stats::ave(x, data[[id]], FUN = function(z) mean(z, na.rm = TRUE))
        1 - stats::var(m, na.rm = TRUE) / stats::var(x, na.rm = TRUE)
    }, numeric(1))
    structure(list(lavaan_fit = fit, paths = paths, subjects = length(keep), within_share = ws),
              class = "biomimic_riclpm")
}

#' @export
print.biomimic_riclpm <- function(x, ...) {
    cat(sprintf("biomimic random-intercept cross-lagged panel model (%d subjects)\n", x$subjects))
    p <- x$paths
    for (i in seq_len(nrow(p))) cat(sprintf("  %-14s %-20s -> %-20s %+.3f (se %.3f, p = %.2g)\n",
        p$type[i], p$from[i], p$to[i], p$estimate[i], p$se[i], p$p[i]))
    cat("  within-subject share of variance:",
        paste(sprintf("%s %.2f", names(x$within_share), x$within_share), collapse = ", "), "\n")
    invisible(x)
}


#' @noRd
.with_seed <- function(seed, expr) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
        assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    force(expr)
}

