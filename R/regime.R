# --------------------------------------------------------------------------
# Which regime are the data in? A diagnostic for Stage 1
# --------------------------------------------------------------------------
# A latent-pathway model is the right tool when the response is co-ordinated
# across variables and the wrong one when it is subject-specific. Stage 1
# cannot tell the difference on its own: it will always return a factor, and
# ranking by loading SNR will always return the variables that load on
# whatever the dominant shared mode happens to be -- biology in one dataset,
# assay structure in another.
#
# The two cases are separable before any model is interpreted. Take the
# variables a channel selects and measure how strongly they still co-vary once
# the design is accounted for, against random variable sets of the same size
# in the same data. A co-ordinated response leaves them correlated well above
# that background; a subject-specific one does not.
# --------------------------------------------------------------------------


#' Co-ordination Diagnostic for a VB-ARD Screen
#'
#' Reports, for each screening channel, how strongly the top-\code{k} selected
#' variables co-vary beyond the background correlation of the assay. This
#' separates a genuinely co-ordinated response -- for which the latent-pathway
#' model in Stage 3 is appropriate -- from a subject-specific one, where the
#' group signal lives in the direct effects and a common factor is the wrong
#' summary.
#'
#' Correlations are computed on residuals after projecting out \code{fit$X}, so
#' that the shared structure measured is what remains once the design (group,
#' covariates) has been accounted for. Supplying \code{reference} instead
#' restricts the correlation to a subset of samples, typically the controls,
#' which is equivalent in intent and easier to explain in a manuscript.
#'
#' No threshold is imposed and no verdict is returned. The ratio is a property
#' of the assay as much as of the biology, and the useful comparison is
#' between channels within one dataset and against the reference values in the
#' examples below.
#'
#' Each channel also gets an empirical p-value, the share of random sets at
#' least as co-ordinated as the panel, \eqn{(1 + \#\{b \ge obs\}) / (1 +
#' n_{null})}. Read it with its null in mind. The lambda channel ranks
#' variables by their loading on a shared factor, so its panel will beat
#' random sets almost by construction whenever any shared mode exists; a small
#' p-value there says little. The gamma channel does not select on
#' co-ordination, and there the p-value is a fair test.
#'
#' A high ratio says a panel is coherent, not that it carries the effect of
#' the design. A nuisance axis -- cell composition, sample collection -- can be
#' more co-ordinated than the signal of interest and will then own the lambda
#' channel. Where a second cohort exists, the deciding check is which
#' channel's selection replicates there, tested at the right level: a lambda
#' panel by the transfer of its factor, a gamma panel indicator by indicator
#' (do its variables keep their direct effects in the other cohort?). A gamma
#' panel typically mixes variables that rise and fall with the design; a
#' common factor is a poor summary of it, and the Wald test of its \code{B}
#' then depends on the marker -- use \code{\link{lrt_structural}} if it is
#' fitted at all.
#'
#' @param fit A fitted \code{vb_ard} object.
#' @param k Integer, panel size to evaluate (default 15, matching
#'   \code{\link{top_k}}).
#' @param reference Optional logical or integer vector selecting the samples to
#'   compute correlations on. If \code{NULL} (default) all samples are used and
#'   \code{fit$X} is projected out first.
#' @param n_null Integer, number of random variable sets forming the background
#'   (default 200).
#'
#' @return An object of class \code{biomimic_regime}: a list with
#'   \describe{
#'     \item{channels}{data.frame, one row per channel, with \code{mean_abs_r},
#'       \code{background}, \code{ratio}, the empirical \code{p_value} and the
#'       panel size actually used}
#'     \item{panels}{named list of the variables evaluated in each channel}
#'     \item{background_all}{the background correlations drawn, per channel}
#'     \item{reference_n}{number of samples the correlations used}
#'   }
#'
#' @examples
#' \donttest{
#' sim <- simulate_mimic(n = 120L, p_total = 150L, p_signal = 15L,
#'                       K = 1L, seed = 1)
#' fit <- vb_ard(sim$Y, sim$X, K = 1L, use_gamma = TRUE, max_iter = 200L)
#'
#' # The simulated signal variables share one factor by construction, so the
#' # lambda channel sits well above the background of random variable sets.
#' coordination_diagnostic(fit, n_null = 50L)
#' }
#'
#' @seealso \code{\link{rank_variables}}, \code{\link{top_k}}
#' @export
coordination_diagnostic <- function(fit, k = 15L, reference = NULL,
                                    n_null = 200L) {
    stopifnot(inherits(fit, "vb_ard"))
    k <- as.integer(k)
    stopifnot(k >= 3L, n_null >= 1L)

    Y <- fit$Y
    if (is.null(colnames(Y))) colnames(Y) <- fit$var_names

    if (is.null(reference)) {
        # Residualise on the design: what co-variation survives the model.
        X <- fit$X
        qrX <- qr(X)
        if (qrX$rank < ncol(X)) X <- X[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
        R <- Y - X %*% qr.solve(X, Y)
        ref_n <- nrow(Y)
    } else {
        R <- Y[reference, , drop = FALSE]
        ref_n <- nrow(R)
    }
    if (ref_n < 5L) {
        stop("need at least 5 samples to estimate correlations; got ", ref_n)
    }

    mean_abs_cor <- function(cols) {
        m <- stats::cor(R[, cols, drop = FALSE], use = "pairwise.complete.obs")
        mean(abs(m[upper.tri(m)]), na.rm = TRUE)
    }

    channels <- c("lambda", if (!is.null(fit$snr_gamma)) "gamma")
    rows <- list(); panels <- list(); bg_all <- list()
    for (ch in channels) {
        rk <- tryCatch(rank_variables(fit, type = ch), error = function(e) NULL)
        if (is.null(rk)) next
        vars <- utils::head(rk$variable, k)
        vars <- intersect(vars, colnames(R))
        if (length(vars) < 3L) next
        obs <- mean_abs_cor(vars)
        bg <- vapply(seq_len(n_null), function(i)
            mean_abs_cor(sample(colnames(R), length(vars))), numeric(1))
        bg <- bg[is.finite(bg)]
        panels[[ch]] <- vars
        bg_all[[ch]] <- bg
        # correlation of the panel mean with each sample's whole-array mean:
        # a panel that tracks it measures the sample's global level
        Ys <- scale(Y)
        r_global <- stats::cor(rowMeans(Ys[, vars, drop = FALSE]),
                               rowMeans(Ys, na.rm = TRUE), use = "complete.obs")
        rows[[length(rows) + 1L]] <- data.frame(
            channel = ch, n_variables = length(vars), mean_abs_r = obs,
            background = stats::median(bg), ratio = obs / stats::median(bg),
            p_value = (1 + sum(bg >= obs)) / (1 + length(bg)),
            r_global = r_global, stringsAsFactors = FALSE)
    }
    if (!length(rows)) stop("no channel produced an evaluable panel")

    structure(list(channels = do.call(rbind, rows), panels = panels,
                   background_all = bg_all, reference_n = ref_n,
                   residualised = is.null(reference), k = k),
              class = "biomimic_regime")
}


#' @param x A \code{biomimic_regime} object.
#' @param ... Ignored.
#' @rdname coordination_diagnostic
#' @export
print.biomimic_regime <- function(x, ...) {
    cat("biomimic co-ordination diagnostic\n")
    cat(sprintf("  %d samples, panels of %d, %s\n", x$reference_n, x$k,
                if (x$residualised) "residualised on the design"
                else "restricted to the reference samples"))
    cat("\n")
    d <- x$channels
    has_p <- "p_value" %in% names(d)
    has_g <- "r_global" %in% names(d)
    cat(sprintf("  %-8s %10s %12s %8s%s%s\n", "channel", "mean |r|",
                "background", "ratio", if (has_p) "   p (empirical)" else "",
                if (has_g) "   r(global)" else ""))
    for (i in seq_len(nrow(d))) {
        cat(sprintf("  %-8s %10.3f %12.3f %8.2f%s%s\n", d$channel[i],
                    d$mean_abs_r[i], d$background[i], d$ratio[i],
                    if (has_p) sprintf("   %15.3g", d$p_value[i]) else "",
                    if (has_g) sprintf("   %10.2f", d$r_global[i]) else ""))
    }
    cat("\n")
    cat("  A ratio near 1 means the panel is no more co-ordinated than an\n")
    cat("  arbitrary set of the same size: the group signal, if any, is in the\n")
    cat("  direct effects rather than in a shared factor. A high ratio says the\n")
    cat("  panel is coherent, not what makes it so: in a quantile-normalised\n")
    cat("  autoantibody array a panel with ratio 5.0 turned out to mirror each\n")
    cat("  serum's global level, and the interferon module of SLE blood had\n")
    cat("  ratios of 1.7 and 4.2 in two cohorts.\n")
    cat("\n")
    cat("  r(global) is the correlation of the panel mean with each sample's\n")
    cat("  whole-array mean. A panel that tracks it measures the sample's global\n")
    cat("  level (intensity, total reactivity, cell composition); check it with\n")
    cat("  global_level_check() and against external measures with\n")
    cat("  external_check(). The lambda p-value is small almost by\n")
    cat("  construction; the gamma one is a fair test.\n")
    invisible(x)
}
