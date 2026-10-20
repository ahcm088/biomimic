# --------------------------------------------------------------------------
# Is a panel measuring biology or each sample's global level?
# --------------------------------------------------------------------------
# A factor fitted to array data takes the largest axis of shared variation,
# whatever it is. In the validation behind this package that axis was the
# intensity of each array, the overall reactivity of each serum, or the cell
# composition of blood, more often than a pathway. These functions make the
# checks that exposed it routine: the panel against each sample's whole-array
# mean, the panel against external measurements, and the analysis repeated
# after per-array normalisation.
# --------------------------------------------------------------------------


#' Panel Scores
#'
#' Standardised mean of each panel's features, per sample.
#'
#' @param Y Numeric matrix, samples x features, with column names.
#' @param panels A character vector or a named list of character vectors, or
#'   a \code{biomimic_selection} (its panels, in the original names).
#' @param standardise Logical, standardise each feature before averaging
#'   (default \code{TRUE}).
#'
#' @return A data frame, one column per panel.
#'
#' @export
panel_scores <- function(Y, panels, standardise = TRUE) {
    panels <- .as_panel_list(panels)
    Z <- if (standardise) scale(Y) else as.matrix(Y)
    out <- lapply(panels, function(p) {
        p <- intersect(p, colnames(Z))
        if (!length(p)) return(rep(NA_real_, nrow(Z)))
        as.numeric(scale(rowMeans(Z[, p, drop = FALSE], na.rm = TRUE)))
    })
    as.data.frame(out, check.names = FALSE)
}


#' Does a Panel Follow Each Sample's Global Level?
#'
#' Correlates a panel's mean with each sample's whole-array mean and
#' re-estimates the group effect on the panel mean with the whole-array mean
#' as a covariate.
#'
#' @param Y Numeric matrix, samples x features (all features of the array).
#' @param panel Character vector, the panel's features.
#' @param group Numeric or logical vector, the group indicator (1 = cases).
#' @param covariates Optional data frame of adjustment covariates.
#' @param reference Optional logical vector; the correlation is also reported
#'   within these samples (typically the controls).
#'
#' @return A list of class \code{biomimic_global}: \code{r_global} (all
#'   samples), \code{r_global_reference}, \code{group_effect} and
#'   \code{group_effect_adjusted} (estimate, 95\% CI, p), and the array-mean
#'   difference between groups.
#'
#' @details A strong correlation is a warning, not a verdict. Each sample's
#'   whole-array mean follows technical intensity, total antibody binding and,
#'   in blood, the proportions of cell types; a panel that tracks the measured
#'   lymphocyte fraction will also track the whole-array mean. Where external
#'   measurements exist, compare the panel with them (\code{\link{external_check}}).
#'
#' @export
global_level_check <- function(Y, panel, group, covariates = NULL,
                               reference = NULL) {
    Z <- scale(as.matrix(Y))
    panel <- intersect(panel, colnames(Z))
    stopifnot(length(panel) >= 2L, length(group) == nrow(Z))
    pm <- as.numeric(scale(rowMeans(Z[, panel, drop = FALSE])))
    am <- as.numeric(scale(rowMeans(Z, na.rm = TRUE)))
    d <- data.frame(pm = pm, am = am, g = as.numeric(group))
    if (!is.null(covariates)) d <- cbind(d, as.data.frame(covariates))
    rhs <- paste(c("g", names(covariates)), collapse = " + ")
    eff <- function(f) {
        cf <- summary(stats::lm(f, d))$coefficients["g", ]
        ci <- cf[1] + c(-1, 1) * stats::qnorm(0.975) * cf[2]
        c(estimate = cf[[1]], lo = ci[1], hi = ci[2], p = cf[[4]])
    }
    out <- list(
        n_panel = length(panel),
        r_global = stats::cor(pm, am),
        r_global_reference = if (is.null(reference)) NA_real_ else
            stats::cor(pm[reference], am[reference]),
        group_effect = eff(stats::as.formula(paste("pm ~", rhs))),
        group_effect_adjusted = eff(stats::as.formula(paste("pm ~ am +", rhs))),
        array_mean_difference = eff(stats::as.formula(paste("am ~", rhs))))
    structure(out, class = "biomimic_global")
}

#' @export
print.biomimic_global <- function(x, ...) {
    cat(sprintf("Global-level check (%d panel features)\n", x$n_panel))
    cat(sprintf("  r(panel mean, whole-array mean) = %.2f%s\n", x$r_global,
                if (is.finite(x$r_global_reference))
                    sprintf(" (%.2f in the reference samples)", x$r_global_reference) else ""))
    f <- function(v) sprintf("%+.3f [%+.3f, %+.3f], p = %.2g", v[["estimate"]],
                             v[["lo"]], v[["hi"]], v[["p"]])
    cat("  group effect on the panel mean:               ", f(x$group_effect), "\n")
    cat("  ... adjusted for the whole-array mean:        ", f(x$group_effect_adjusted), "\n")
    cat("  group difference in the whole-array mean:     ", f(x$array_mean_difference), "\n")
    invisible(x)
}


#' Panel Scores Against External Measurements
#'
#' Spearman correlation of each panel score with each external measurement
#' (blood counts, clinical scores, serology), the check that identified the
#' lymphopenia factor of SLE blood.
#'
#' @param scores Data frame or matrix of panel scores (e.g. from
#'   \code{\link{panel_scores}}).
#' @param measures Data frame of external measurements, one row per sample.
#' @param method Correlation method (default \code{"spearman"}).
#'
#' @return A data frame with one row per score and measure: \code{rho},
#'   \code{p} and the number of complete samples.
#'
#' @export
external_check <- function(scores, measures, method = "spearman") {
    scores <- as.data.frame(scores); measures <- as.data.frame(measures)
    stopifnot(nrow(scores) == nrow(measures))
    rows <- list()
    for (s in names(scores)) for (m in names(measures)) {
        ok <- is.finite(scores[[s]]) & is.finite(measures[[m]])
        if (sum(ok) < 5L) next
        ct <- suppressWarnings(stats::cor.test(scores[[s]][ok], measures[[m]][ok],
                                               method = method, exact = FALSE))
        rows[[length(rows) + 1L]] <- data.frame(score = s, measure = m, n = sum(ok),
                                                rho = unname(ct$estimate), p = ct$p.value,
                                                stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}


#' Per-Array Normalisation
#'
#' Centres each sample on its median, or quantile-normalises the samples, for
#' the re-analysis check: a group effect that vanishes after per-array
#' normalisation was carried by each array's overall level.
#'
#' @param Y Numeric matrix, samples x features (log scale).
#' @param method \code{"median"} (default) or \code{"quantile"}.
#'
#' @return The normalised matrix, same dimensions and names.
#'
#' @details Quantile normalisation forces every sample to the same
#'   distribution. If some features rise in a sample, the others must fall in
#'   rank, which can create a block of features that move together against
#'   the rest of the array; see the thesis, Chapter 3, on mirrored blocks.
#'
#' @export
normalise_arrays <- function(Y, method = c("median", "quantile")) {
    method <- match.arg(method)
    Y <- as.matrix(Y)
    if (method == "median") {
        return(sweep(Y, 1, apply(Y, 1, stats::median, na.rm = TRUE)))
    }
    if (anyNA(Y)) stop("quantile normalisation needs complete data", call. = FALSE)
    o <- t(apply(Y, 1, order))
    ref <- rowMeans(apply(Y, 1, sort))            # mean of the sorted samples
    out <- Y
    for (i in seq_len(nrow(Y))) out[i, o[i, ]] <- ref
    out
}


#' @noRd
.as_panel_list <- function(panels) {
    if (inherits(panels, "biomimic_selection")) {
        p <- panels$panels %||% list(F1 = panels$variables)
        lab <- panels$labels
        if (!is.null(lab)) p <- lapply(p, function(v) unname(lab[v]))
        return(p)
    }
    if (is.character(panels)) return(list(panel = panels))
    stopifnot(is.list(panels))
    if (is.null(names(panels))) names(panels) <- paste0("panel", seq_along(panels))
    panels
}
