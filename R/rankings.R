# --------------------------------------------------------------------------
# Feature rankings from established methods
# Alternatives to the VB-ARD rankings for the first stage of the pipeline.
# --------------------------------------------------------------------------

#' Rank Features with limma
#'
#' Moderated t-tests (\pkg{limma}) of one design coefficient, returned in the
#' format of \code{\link{rank_variables}} so that any ranking can feed the
#' rest of the pipeline.
#'
#' @param Y Numeric matrix, samples x features, with column names.
#' @param X Design matrix (with intercept) or data frame of covariates.
#' @param coef Character, the column of \code{X} to test (default
#'   \code{"group"}).
#'
#' @return A data frame sorted by \code{|t|}: \code{variable}, \code{snr}
#'   (\code{|t|}), \code{rank}, \code{estimate}, \code{t}, \code{p},
#'   \code{adj_p} (Benjamini-Hochberg).
#'
#' @export
rank_limma <- function(Y, X, coef = "group") {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("rank_limma() needs the limma package (Bioconductor).", call. = FALSE)
    }
    X <- .as_design(X)
    stopifnot(coef %in% colnames(X))
    fit <- limma::eBayes(limma::lmFit(t(Y), X))
    tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
    .ranking_frame(colnames(Y), tt$logFC, tt$t, tt$P.Value)
}


#' Rank Features with CATE
#'
#' Confounder-adjusted testing and estimation (Wang et al. 2017): the effect
#' of one covariate on each feature after removing \code{r} latent factors,
#' estimated from the residuals with a robust-regression adjustment. Needs
#' the \pkg{cate} package, which is archived on CRAN and must be installed
#' from the archive (with its dependency \pkg{esaBcv}).
#'
#' @param Y Numeric matrix, samples x features.
#' @param covariates Data frame with the covariate of interest and any
#'   nuisance covariates.
#' @param coef Character, the covariate of interest (default \code{"group"}).
#' @param r Integer, number of latent factors, or \code{NULL} (default) to
#'   estimate it with \code{cate::est.confounder.num(method = "ed")}.
#'
#' @return A data frame in the format of \code{\link{rank_limma}}, with the
#'   number of factors used in attribute \code{"r"}.
#'
#' @details In the thesis validation, CATE with its own estimate of \code{r}
#'   replicated age effects in DNA methylation at least as well as the
#'   VB-ARD gamma ranking with three factors, and was much faster.
#'
#' @export
rank_cate <- function(Y, covariates, coef = "group", r = NULL) {
    if (!requireNamespace("cate", quietly = TRUE)) {
        stop("rank_cate() needs the cate package, archived on CRAN; install it ",
             "from https://cran.r-project.org/src/contrib/Archive/cate/ ",
             "(with esaBcv from the same archive).", call. = FALSE)
    }
    covariates <- as.data.frame(covariates)
    stopifnot(coef %in% names(covariates), nrow(covariates) == nrow(Y))
    nuis <- setdiff(names(covariates), coef)
    f <- stats::as.formula(paste("~", coef, "|",
                                 if (length(nuis)) paste(nuis, collapse = " + ") else "1"))
    if (is.null(r)) r <- cate::est.confounder.num(f, X.data = covariates, Y = Y, method = "ed")
    fit <- cate::cate(f, X.data = covariates, Y = Y, r = r, adj.method = "rr", calibrate = FALSE)
    out <- .ranking_frame(colnames(Y), as.vector(fit$beta), as.vector(fit$beta.t),
                          as.vector(fit$beta.p.value))
    attr(out, "r") <- r
    out
}


#' @noRd
.ranking_frame <- function(vars, est, t, p) {
    d <- data.frame(variable = vars, snr = abs(t), estimate = est, t = t, p = p,
                    adj_p = stats::p.adjust(p, "BH"), stringsAsFactors = FALSE)
    d <- d[order(-d$snr), ]
    d$rank <- seq_len(nrow(d))
    rownames(d) <- NULL
    d[, c("variable", "snr", "rank", "estimate", "t", "p", "adj_p")]
}

#' @noRd
.as_design <- function(X) {
    if (is.data.frame(X)) {
        X <- stats::model.matrix(~ ., X)
        colnames(X)[1] <- "intercept"
    }
    X <- as.matrix(X)
    if (!any(apply(X, 2, function(v) all(v == 1)))) X <- cbind(intercept = 1, X)
    X
}
