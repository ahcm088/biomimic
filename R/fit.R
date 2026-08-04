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
#' @param ... Additional arguments passed to \code{lavaan::sem()}.
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
                    ...) {

    type <- match.arg(type)
    stopifnot(inherits(selection, "biomimic_selection"))

    # Build model syntax
    model_syntax <- build_lavaan_model(selection, type = type,
                                        group_var = group_var,
                                        covariate_names = covariate_names)

    # Determine K
    K <- selection$K
    if (type == "mimic" || type == "direct") K <- 1L

    # Default factor labels
    if (is.null(factor_labels)) {
        factor_labels <- if (K == 1L) {
            "Latent pathway"
        } else {
            paste("Pathway", seq_len(K))
        }
    }
    stopifnot(length(factor_labels) == K)

    # Default indicator labels (use variable names)
    if (is.null(indicator_labels)) {
        indicator_labels <- setNames(selection$variables, selection$variables)
    }

    # Fit with lavaan
    .log_info("fit_sem: fitting %s model (K=%d, estimator=%s, p=%d indicators)",
              type, K, estimator, length(selection$variables))
    .log_debug("fit_sem: model syntax:\n%s", model_syntax)

    lavaan_args <- list(model = model_syntax, data = selection$data,
                        estimator = estimator, se = se, ...)
    if (!is.null(group_var)) {
        lavaan_args$group <- group_var
    }
    lavaan_fit <- do.call(lavaan::sem, lavaan_args)

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
            K                = K
        ),
        class = "biomimic_fit"
    )
}
