# --------------------------------------------------------------------------
# SEM Fitting Wrapper
# Stage 3 of the biomimic pipeline — high-level fitting interface
# --------------------------------------------------------------------------

#' Fit a Structural Equation Model from biomimic Selection
#'
#' Wraps \code{lavaan::sem()} with sensible defaults for biological data.
#' Accepts the output of Stage 2 (\code{top_k()}) and returns a
#' \code{biomimic_fit} object with biology-friendly methods.
#'
#' @param selection A \code{biomimic_selection} object from \code{top_k()}.
#' @param type Character, model type: \code{"mimic"}, \code{"multigroup"},
#'   \code{"direct"}, or \code{"multi_mimic"}. See \code{build_lavaan_model()}.
#' @param estimator Character, lavaan estimator (default \code{"MLR"} for
#'   robust maximum likelihood with Huber-White standard errors).
#' @param group_var Character, grouping variable name (required for
#'   \code{type = "multigroup"}).
#' @param covariate_names Character vector of covariate names. If \code{NULL},
#'   all non-intercept columns of X are used.
#' @param se Character, standard error method (default \code{"robust"}).
#' @param indicator_labels Named character vector mapping variable names to
#'   biology-friendly labels (e.g., \code{c(V1 = "IL-6", V2 = "TNF-alpha")}).
#'   If \code{NULL}, original variable names are used.
#' @param factor_labels Character vector of names for the latent factors
#'   (e.g., \code{"Immune activation"}). Length must match K.
#' @param orient Logical (default \code{TRUE}). If a factor's marker loads
#'   against the majority of its indicators, refit with a majority-sign marker
#'   so that the sign of \code{B} describes the majority of the panel. See
#'   Details.
#' @param cluster Optional vector, one subject identifier per sample, for
#'   repeated measures: lavaan then reports cluster-robust standard errors
#'   and test statistics. Combine with \code{\link{within_between}} to
#'   separate within-subject from between-subject effects.
#' @param ... Additional arguments passed to \code{lavaan::sem()}.
#'
#' @details
#' lavaan scales each factor by fixing its first (marker) loading to +1, so
#' the factor points the way the marker does. When the marker anti-correlates
#' with most of the panel, every loading but the marker's comes out negative
#' and \code{B} changes sign: "lower in patients" and "higher in patients"
#' then depend on which indicator happened to be listed first. The fit itself
#' is unchanged by the reflection; only the reading of the sign is.
#'
#' With \code{orient = TRUE} the marker of any such factor is replaced by the
#' strongest (standardised) majority-sign indicator and the model is refitted.
#' For \code{type = "direct"} the replacement is drawn only from indicators
#' without a direct effect, which identification requires of the marker. The
#' orientation actually used is returned in \code{orientation}.
#'
#' Orientation fixes the sign of \code{B}, not its test. The Wald z that
#' lavaan reports for a marker-scaled coefficient depends on which indicator
#' is the marker; on a co-ordinated panel barely, on a panel whose loadings
#' split in sign enough to decide significance. \code{\link{lrt_structural}}
#' gives the marker-invariant likelihood-ratio test.
#'
#' @return A \code{biomimic_fit} object (S3 class) with components:
#' \describe{
#'   \item{lavaan_fit}{The underlying \code{lavaan} object (full access)}
#'   \item{selection}{The Stage 2 selection object}
#'   \item{model_type}{Character, which model was fitted}
#'   \item{model_syntax}{Character, the lavaan model syntax}
#'   \item{indicator_labels}{Named character vector of biology-friendly labels}
#'   \item{factor_labels}{Character vector of factor names}
#'   \item{K}{Number of latent factors}
#'   \item{orientation}{Data frame, one row per factor: \code{marker},
#'     counts of positive and negative loadings, whether the marker sides
#'     with the majority, and the marker it replaced (\code{NA} if none)}
#' }
#'
#' @examples
#' \dontrun{
#' sel <- top_k(vb_fit, k = 15)
#' fit <- fit_sem(sel, type = "mimic")
#' summary(fit)
#' plot(fit, type = "loadings")
#' }
#'
#' @export
fit_sem <- function(selection,
                    type = c("mimic", "multigroup", "direct", "multi_mimic"),
                    estimator = "MLR",
                    group_var = NULL,
                    covariate_names = NULL,
                    se = "robust",
                    indicator_labels = NULL,
                    factor_labels = NULL,
                    orient = TRUE,
                    cluster = NULL,
                    ...) {

    type <- match.arg(type)
    stopifnot(inherits(selection, "biomimic_selection"))
    if (!is.null(cluster)) {
        stopifnot(length(cluster) == nrow(selection$data))
        selection$data$.cluster <- cluster
        selection$cluster <- ".cluster"
    }

    # Determine K
    K <- selection$K
    if (type == "mimic" || type == "direct") K <- 1L

    # Default factor labels: panel names for selections built by
    # as_selection(), generic names otherwise
    if (is.null(factor_labels)) {
        factor_labels <- if (!is.null(selection$panels) &&
                             length(selection$panels) == K) {
            names(selection$panels)
        } else if (K == 1L) {
            "Latent pathway"
        } else {
            paste("Pathway", seq_len(K))
        }
    }
    stopifnot(length(factor_labels) == K)

    # Default indicator labels: the original feature names when the
    # selection renamed them for lavaan, the variable names otherwise
    if (is.null(indicator_labels)) {
        indicator_labels <- selection$labels %||%
            setNames(selection$variables, selection$variables)
    }

    # Fit with lavaan
    .log_info("fit_sem: fitting %s model (K=%d, estimator=%s, p=%d indicators)",
              type, K, estimator, length(selection$variables))

    extra_args <- list(...)
    fit_once <- function(sel) {
        syntax <- build_lavaan_model(sel, type = type, group_var = group_var,
                                     covariate_names = covariate_names)
        .log_debug("fit_sem: model syntax:\n%s", syntax)
        args <- c(list(model = syntax, data = sel$data,
                       estimator = estimator, se = se), extra_args)
        if (!is.null(group_var)) args$group <- group_var
        if (!is.null(sel$cluster)) args$cluster <- sel$cluster
        list(syntax = syntax, fit = do.call(lavaan::sem, args))
    }

    res <- fit_once(selection)
    orientation <- .loading_orientation(res$fit)
    if (!is.null(orientation)) orientation$replaced <- NA_character_

    if (isTRUE(orient) && !is.null(orientation) &&
        any(!orientation$aligned)) {
        re <- .reorient_selection(selection, res$fit, orientation, type)
        if (length(re$replaced)) {
            for (f in names(re$replaced)) {
                .log_info(paste0("fit_sem: %s marker %s loaded against the ",
                                 "majority of its indicators; refitting with ",
                                 "marker %s"),
                          f, re$replaced[[f]], re$new_marker[[f]])
            }
            selection <- re$selection
            res <- fit_once(selection)
            orientation <- .loading_orientation(res$fit)
            orientation$replaced <- unname(
                re$replaced[orientation$factor])
            orientation$replaced[!orientation$factor %in%
                                     names(re$replaced)] <- NA_character_
            if (any(!orientation$aligned)) {
                .log_warn("fit_sem: reorientation did not align factor(s) %s",
                          paste(orientation$factor[!orientation$aligned],
                                collapse = ", "))
            }
        }
    }
    if (!is.null(orientation) && any(!orientation$aligned)) {
        .log_warn(paste0("fit_sem: marker of %s loads against the majority of ",
                         "its indicators; the sign of B describes the ",
                         "minority. Read B with the loading signs."),
                  paste(orientation$factor[!orientation$aligned],
                        collapse = ", "))
    }

    model_syntax <- res$syntax
    lavaan_fit   <- res$fit

    converged <- lavaan::lavInspect(lavaan_fit, "converged")
    if (converged) {
        fm <- lavaan::fitMeasures(lavaan_fit, c("rmsea", "cfi", "srmr"))
        .log_info("fit_sem: converged. RMSEA=%.3f, CFI=%.3f, SRMR=%.3f",
                  fm["rmsea"], fm["cfi"], fm["srmr"])
    } else {
        .log_error("fit_sem: lavaan did NOT converge")
    }

    structure(
        list(
            lavaan_fit       = lavaan_fit,
            selection        = selection,
            model_type       = type,
            model_syntax     = model_syntax,
            indicator_labels = indicator_labels,
            factor_labels    = factor_labels,
            K                = K,
            orientation      = orientation
        ),
        class = "biomimic_fit"
    )
}


#' @noRd
.loading_orientation <- function(lavaan_fit) {
    # One row per factor: which indicator is the fixed marker, and whether it
    # sides with the majority of the factor's loadings (ties count as aligned).
    pt <- lavaan::parTable(lavaan_fit)
    ld <- pt[pt$op == "=~", , drop = FALSE]
    if (!nrow(ld)) return(NULL)
    if ("group" %in% names(ld)) ld <- ld[ld$group == min(ld$group), ,
                                         drop = FALSE]
    rows <- lapply(unique(ld$lhs), function(f) {
        r <- ld[ld$lhs == f, , drop = FALSE]
        fixed <- r$rhs[r$free == 0L]
        n_pos <- sum(r$est > 0)
        n_neg <- sum(r$est < 0)
        data.frame(factor = f,
                   marker = if (length(fixed)) fixed[1] else NA_character_,
                   n_pos = n_pos, n_neg = n_neg,
                   aligned = n_pos >= n_neg,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' @noRd
.reorient_selection <- function(selection, lavaan_fit, orientation, type) {
    # For each misaligned factor, move its strongest majority-sign indicator
    # (by |standardised loading|) to the marker position.
    ss <- lavaan::standardizedSolution(lavaan_fit)
    ss <- ss[ss$op == "=~", , drop = FALSE]
    if ("group" %in% names(ss)) ss <- ss[ss$group == min(ss$group), ,
                                         drop = FALSE]
    ineligible <- if (type == "direct") selection$direct_effects else
        character(0)

    replaced <- list()
    new_marker <- list()
    for (i in which(!orientation$aligned)) {
        f <- orientation$factor[i]
        r <- ss[ss$lhs == f & ss$est.std < 0 & !ss$rhs %in% ineligible, ,
                drop = FALSE]
        if (!nrow(r)) next
        nm <- r$rhs[which.max(abs(r$est.std))]

        if (type == "multi_mimic") {
            k <- as.integer(sub("^F", "", f))
            idx <- match(nm, selection$variables)
            g <- selection$indicator_groups[[k]]
            selection$indicator_groups[[k]] <- c(idx, setdiff(g, idx))
        } else {
            selection$variables <- c(nm, setdiff(selection$variables, nm))
        }
        replaced[[f]] <- orientation$marker[i]
        new_marker[[f]] <- nm
    }
    list(selection = selection, replaced = unlist(replaced),
         new_marker = new_marker)
}
