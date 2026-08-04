# --------------------------------------------------------------------------
# Predictions: Factor Scores and Observed-Level Predictions
# --------------------------------------------------------------------------

#' Predict Factor Scores
#'
#' Estimates latent factor scores for each observation. These represent the
#' estimated value of the latent biological pathway for each individual,
#' enabling downstream analyses (e.g., survival models, clustering,
#' visualisation of group separation).
#'
#' @param fit A \code{biomimic_fit} object.
#' @param method Character, scoring method:
#'   \code{"regression"} (Thurstone/Thomson, default) or \code{"Bartlett"}.
#' @param newdata Optional data frame for prediction on new data.
#'   If \code{NULL}, scores are computed for the training data.
#' @param append_data Logical, whether to append the factor scores to the
#'   original data (default \code{TRUE}).
#'
#' @return A data frame with factor score columns named after
#'   \code{factor_labels}, optionally appended to the original data.
#'
#' @export
predict_scores <- function(fit, method = c("regression", "Bartlett"),
                            newdata = NULL, append_data = TRUE) {
    method <- match.arg(method)
    stopifnot(inherits(fit, "biomimic_fit"))

    if (!is.null(newdata)) {
        scores <- lavaan::lavPredict(fit$lavaan_fit, newdata = newdata,
                                      method = method)
    } else {
        scores <- lavaan::lavPredict(fit$lavaan_fit, method = method)
    }

    scores <- as.data.frame(scores)

    # Rename columns to biology-friendly labels
    K <- fit$K
    if (ncol(scores) == K) {
        colnames(scores) <- make.names(fit$factor_labels)
    }

    if (append_data) {
        dat <- if (!is.null(newdata)) newdata else fit$selection$data
        scores <- cbind(dat, scores)
    }

    scores
}


#' Predict Observed Variables from the Model
#'
#' Computes model-implied predicted values for each observed variable.
#' The residual (observed - predicted) highlights which proteins are
#' poorly explained by the latent factor model for each individual.
#'
#' @param fit A \code{biomimic_fit} object.
#'
#' @return A list with:
#' \describe{
#'   \item{predicted}{Matrix of predicted values (n x p)}
#'   \item{residuals}{Matrix of residuals (n x p)}
#'   \item{data}{The original observed data}
#' }
#'
#' @export
predict_observed <- function(fit) {
    stopifnot(inherits(fit, "biomimic_fit"))

    vars <- fit$selection$variables
    observed <- as.matrix(fit$selection$data[, vars, drop = FALSE])

    # Get model-implied values via factor scores and loadings
    scores <- lavaan::lavPredict(fit$lavaan_fit, method = "regression")
    scores <- as.matrix(scores)

    # Extract estimated loadings and intercepts
    pe <- lavaan::parameterEstimates(fit$lavaan_fit)
    loadings <- pe[pe$op == "=~", , drop = FALSE]

    # Build loading matrix
    K <- fit$K
    p <- length(vars)
    Lambda <- matrix(0, p, K)
    rownames(Lambda) <- vars
    for (i in seq_len(nrow(loadings))) {
        j <- match(loadings$rhs[i], vars)
        k <- as.integer(sub("F", "", loadings$lhs[i]))
        if (!is.na(j) && !is.na(k) && k <= K) {
            Lambda[j, k] <- loadings$est[i]
        }
    }

    # Extract intercepts
    intercepts <- pe[pe$op == "~1" & pe$lhs %in% vars, , drop = FALSE]
    nu <- setNames(rep(0, p), vars)
    if (nrow(intercepts) > 0) {
        for (i in seq_len(nrow(intercepts))) {
            idx <- match(intercepts$lhs[i], vars)
            if (!is.na(idx)) nu[idx] <- intercepts$est[i]
        }
    }

    # Predicted = nu + Lambda %*% scores'
    predicted <- t(nu + Lambda %*% t(scores))
    colnames(predicted) <- vars

    # Add direct-effect paths (indicator ~ covariates) for the 'direct' model:
    # the model-implied mean of a direct indicator is nu_j + lambda_j' F_i +
    # gamma_j' x_i; without this term the residuals of direct indicators are
    # systematically wrong.
    reg <- pe[pe$op == "~" & pe$lhs %in% vars, , drop = FALSE]
    if (nrow(reg) > 0) {
        dat <- fit$selection$data
        for (i in seq_len(nrow(reg))) {
            jj <- match(reg$lhs[i], vars)
            cov_name <- reg$rhs[i]
            if (!is.na(jj) && cov_name %in% colnames(dat)) {
                predicted[, jj] <- predicted[, jj] + reg$est[i] * dat[[cov_name]]
            }
        }
    }

    residuals <- observed - predicted

    list(
        predicted = predicted,
        residuals = residuals,
        data      = observed
    )
}
