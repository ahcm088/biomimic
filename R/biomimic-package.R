# --------------------------------------------------------------------------
# Package-level imports and documentation
# --------------------------------------------------------------------------

#' @keywords internal
#' @importFrom stats setNames varimax
#' @importFrom ggplot2 %+replace%
#' @importFrom utils head
"_PACKAGE"

# Suppress R CMD check NOTE about .data pronoun from ggplot2
# .data is used in aes() calls for tidy evaluation
if (getRversion() >= "2.15.1") utils::globalVariables(".data")
