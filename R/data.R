#' Neurodegenerative antibodyome (subset)
#'
#' A subset of the autoantibody protein-microarray (ProtoArray) dataset
#' analysed in the biomimic paper, for the Alzheimer's disease (AD) versus
#' control comparison. It is ready for the VB-ARD screening pipeline: features
#' are centred, there are no missing values, and a covariate matrix is
#' supplied.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{proteins}{Numeric matrix, 137 samples \eqn{\times} 400 autoantibody
#'     features (centred reactivities). Column names are gene symbols; rows are
#'     samples. Includes the 15 features selected by VB-ARD in the paper plus a
#'     reproducible random background sample.}
#'   \item{group}{Integer vector of length 137: \code{1} = AD case (n = 46),
#'     \code{0} = control (n = 91).}
#'   \item{X}{Numeric matrix, 137 \eqn{\times} 4: a ready-to-use covariate
#'     matrix with an \code{intercept} column, the \code{group} indicator, and
#'     \code{age} and \code{sex} adjustment covariates, to pass directly as the
#'     \code{X} argument of \code{\link{vb_ard}} or \code{\link{biomimic}}.}
#' }
#'
#' @details This is a CRAN-sized subset (400 of the 8,323 features) intended
#'   for examples and for reproducing the package figures; it preserves the
#'   AD group signal and the age/sex covariates. The full dataset is available
#'   from the package's GitHub repository.
#'
#' @source Gene Expression Omnibus accession \code{GSE62283}.
#'
#' @examples
#' data(neuro_antibodyome)
#' str(neuro_antibodyome, max.level = 1)
#' \donttest{
#' # The data are model-ready (centred, no missing values); X is supplied.
#' fit <- vb_ard(neuro_antibodyome$proteins, neuro_antibodyome$X, K = 1)
#' }
#' @name neuro_antibodyome
#' @docType data
#' @keywords datasets
"neuro_antibodyome"
