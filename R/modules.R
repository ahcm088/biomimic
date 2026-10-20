# --------------------------------------------------------------------------
# Modules and selections
# Front ends of the pipeline: gene sets, WGCNA modules or VB-ARD factors,
# turned into the biomimic_selection objects that fit_sem() and the other
# Stage-3 functions accept.
# --------------------------------------------------------------------------

#' Modules from Curated Gene Sets
#'
#' Turns named gene sets into a \code{biomimic_modules} object. Within each
#' set, genes are ranked by their correlation with the set's first principal
#' component (module membership), so that \code{\link{select_panels}} can
#' take the most representative genes.
#'
#' @param sets Named list of character vectors (gene identifiers).
#' @param Y Numeric matrix, samples x features, with column names.
#' @param min_size Integer, sets with fewer genes present in \code{Y} are
#'   dropped with a message (default 3).
#'
#' @return A \code{biomimic_modules} object: a list with \code{modules}
#'   (one element per set, each with \code{name}, \code{genes} ranked by
#'   membership, \code{membership} and \code{score}) and \code{method}.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(50 * 20), 50, 20, dimnames = list(NULL, paste0("g", 1:20)))
#' mods <- modules_from_sets(list(A = paste0("g", 1:5), B = paste0("g", 6:10)), Y)
#' mods
#'
#' @export
modules_from_sets <- function(sets, Y, min_size = 3L) {
    stopifnot(is.list(sets), !is.null(names(sets)), is.matrix(Y),
              !is.null(colnames(Y)))
    mods <- list()
    for (nm in names(sets)) {
        g <- intersect(sets[[nm]], colnames(Y))
        if (length(g) < min_size) {
            .log_info("modules_from_sets: set %s has %d genes in Y, dropped",
                      nm, length(g))
            next
        }
        mods[[nm]] <- .module_from_genes(nm, g, Y)
    }
    .new_modules(mods, method = "sets")
}


#' Gene Modules by WGCNA
#'
#' Runs a signed weighted gene co-expression network (WGCNA) on the features
#' of highest variance and returns its modules, with genes ranked by module
#' membership (kME, correlation with the module eigengene). Requires the
#' \pkg{WGCNA} package.
#'
#' @param Y Numeric matrix, samples x features, with column names. Values may
#'   be standardised; pass the variance of the unstandardised data in
#'   \code{raw_var} so that the variance filter is meaningful.
#' @param n_genes Integer, number of highest-variance features entering the
#'   network (default 5000).
#' @param power Integer soft-thresholding power, or \code{NULL} (default) for
#'   the smallest power with scale-free fit \eqn{R^2 \ge 0.8}, 12 if none.
#' @param raw_var Optional named numeric vector, the variance of each feature
#'   before standardisation.
#' @param min_size Integer, minimum module size (default 30).
#' @param merge_cut Numeric, eigengene dissimilarity below which modules are
#'   merged (default 0.25).
#' @param ... Further arguments passed to \code{WGCNA::blockwiseModules()}.
#'
#' @return A \code{biomimic_modules} object (see \code{\link{modules_from_sets}});
#'   the soft power used is stored in \code{power}. The unassigned (grey)
#'   module is dropped.
#'
#' @details In the SLE blood validation of the thesis (amendment A14), this
#'   front end found the interferon module in four of four cohorts, against
#'   two for \code{\link{modules_ard}}.
#'
#' @export
modules_wgcna <- function(Y, n_genes = 5000L, power = NULL, raw_var = NULL,
                          min_size = 30L, merge_cut = 0.25, ...) {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
        stop("modules_wgcna() needs the WGCNA package: install.packages('WGCNA').",
             call. = FALSE)
    }
    stopifnot(is.matrix(Y), !is.null(colnames(Y)))
    v <- if (is.null(raw_var)) apply(Y, 2, stats::var) else raw_var[colnames(Y)]
    if (is.null(raw_var) && stats::sd(v) < 1e-8) {
        warning("All features have the same variance (standardised data?); ",
                "pass raw_var for the variance filter.", call. = FALSE)
    }
    G <- names(sort(v, decreasing = TRUE))[seq_len(min(n_genes, ncol(Y)))]
    Z <- scale(Y[, G])
    if (is.null(power)) {
        # pickSoftThreshold() prints its table whatever the verbosity
        utils::capture.output(sft <- suppressWarnings(WGCNA::pickSoftThreshold(
            Z, networkType = "signed", verbose = 0,
            powerVector = c(1:10, seq(12, 20, 2)))))
        fi <- sft$fitIndices
        ok <- which(-sign(fi$slope) * fi$SFT.R.sq >= 0.8)
        power <- if (length(ok)) fi$Power[min(ok)] else 12L
    }
    # blockwiseModules() looks up `cor` in the calling frame and needs
    # WGCNA's version, which takes arguments that stats::cor() does not
    cor <- WGCNA::cor
    net <- WGCNA::blockwiseModules(Z, power = power, networkType = "signed",
                                   TOMType = "signed", minModuleSize = min_size,
                                   mergeCutHeight = merge_cut,
                                   numericLabels = TRUE,
                                   maxBlockSize = max(5000L, length(G) + 1L),
                                   verbose = 0, ...)
    labs <- setdiff(sort(unique(net$colors)), 0)
    mods <- lapply(labs, function(m) .module_from_genes(paste0("M", m),
                                                      G[net$colors == m], Z))
    names(mods) <- paste0("M", labs)
    out <- .new_modules(mods, method = "WGCNA")
    out$power <- power
    out
}


#' Modules from the Factors of a VB-ARD Fit
#'
#' Each factor of a \code{\link{vb_ard}} fit becomes a module, with genes
#' ranked by their loading signal-to-noise ratio on that factor.
#'
#' @param fit A \code{vb_ard} object.
#' @param n Integer, number of top genes kept per module (default 50).
#'
#' @return A \code{biomimic_modules} object.
#'
#' @details This is the screening front end that biomimic introduced; it
#'   estimates the factors and each feature's own (direct) effect in one
#'   model. It is offered as one option among others and is not claimed to
#'   be better: in the thesis validation it matched CATE as an adjustment and
#'   found the SLE interferon module in two of four cohorts, against four
#'   for \code{\link{modules_wgcna}}.
#'
#' @export
modules_ard <- function(fit, n = 50L) {
    stopifnot(inherits(fit, "vb_ard"))
    mods <- lapply(seq_len(fit$K), function(k) {
        s <- fit$snr_lambda[, k]
        names(s) <- fit$var_names
        o <- order(s, decreasing = TRUE)[seq_len(min(n, length(s)))]
        sc <- fit$mu_eta[, k]
        g <- fit$var_names[o]
        # orient so that most top genes load positively
        lam <- fit$mu_lambda[o, k]
        if (sum(lam > 0) < sum(lam < 0)) sc <- -sc
        list(name = paste0("F", k), genes = g, membership = unname(s[o]),
             score = sc)
    })
    names(mods) <- paste0("F", seq_len(fit$K))
    .new_modules(mods, method = "VB-ARD")
}


#' Label Modules by Marker Content
#'
#' For each named marker set, finds the module with the most markers among
#' its top genes. This is how modules were assigned in the thesis validation
#' (interferon-stimulated genes, plasma-cell, granulopoiesis markers).
#'
#' @param modules A \code{biomimic_modules} object.
#' @param markers Named list of character vectors.
#' @param top Integer, how many top genes of each module are searched
#'   (default 50).
#'
#' @return A data frame with one row per marker set: the best module, its
#'   number of markers among the top genes, and the counts for every module.
#'
#' @export
annotate_modules <- function(modules, markers, top = 50L) {
    stopifnot(inherits(modules, "biomimic_modules"), is.list(markers),
              !is.null(names(markers)))
    hits <- vapply(modules$modules, function(m)
        vapply(markers, function(s) sum(s %in% utils::head(m$genes, top)),
               numeric(1)), numeric(length(markers)))
    hits <- matrix(hits, nrow = length(markers),
                   dimnames = list(names(markers), names(modules$modules)))
    best <- apply(hits, 1, which.max)
    data.frame(markers = names(markers),
               module = colnames(hits)[best],
               hits = hits[cbind(seq_along(best), best)],
               counts = apply(hits, 1, paste, collapse = ";"),
               stringsAsFactors = FALSE, row.names = NULL)
}


#' Build a Selection from Named Panels
#'
#' Creates a \code{biomimic_selection} from one or several panels of features
#' chosen outside the VB-ARD screening (curated sets, WGCNA modules, limma or
#' CATE rankings), so that every Stage-3 function can use it.
#'
#' @param Y Numeric matrix, samples x features, with column names.
#' @param panels A character vector (one factor) or a named list of
#'   character vectors (one factor each; panels must not overlap).
#' @param covariates Optional data frame with one row per sample (group,
#'   adjustment covariates); its columns become the structural covariates.
#'
#' @return A \code{biomimic_selection} with \code{fit = NULL}, the panels
#'   (in lavaan-safe names), the original labels and the covariate names.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(60 * 6), 60, 6, dimnames = list(NULL, paste0("g", 1:6)))
#' sel <- as_selection(Y, list(A = c("g1", "g2", "g3")),
#'                     covariates = data.frame(group = rbinom(60, 1, 0.5)))
#' sel
#'
#' @export
as_selection <- function(Y, panels, covariates = NULL) {
    stopifnot(is.matrix(Y) || is.data.frame(Y))
    Y <- as.matrix(Y)
    if (is.character(panels)) panels <- list(F1 = panels)
    stopifnot(is.list(panels), length(panels) >= 1L)
    if (is.null(names(panels))) names(panels) <- paste0("F", seq_along(panels))
    panels <- lapply(panels, function(p) intersect(p, colnames(Y)))
    if (any(lengths(panels) < 2L)) {
        stop("Every panel needs at least 2 features present in Y.", call. = FALSE)
    }
    vars <- unlist(panels, use.names = FALSE)
    if (anyDuplicated(vars)) {
        stop("Panels overlap; a feature can load on one panel only.", call. = FALSE)
    }
    safe <- make.names(vars, unique = TRUE)
    labels <- stats::setNames(vars, safe)
    data <- stats::setNames(as.data.frame(Y[, vars, drop = FALSE]), safe)
    cov_names <- character(0)
    if (!is.null(covariates)) {
        covariates <- as.data.frame(covariates)
        stopifnot(nrow(covariates) == nrow(Y))
        cov_names <- names(covariates)
        data <- cbind(data, covariates)
    }
    K <- length(panels)
    sizes <- lengths(panels)
    groups <- if (K >= 2L) split(seq_along(vars), rep(seq_len(K), sizes)) else NULL
    if (!is.null(groups)) names(groups) <- NULL
    structure(list(
        variables = safe, ranking = NULL, data = data, K = K,
        indicator_groups = groups, direct_effects = character(0),
        gamma_threshold = NULL, gamma_calibration = NULL, fit = NULL,
        panels = stats::setNames(split(safe, rep(names(panels), sizes))[names(panels)],
                                 names(panels)),
        labels = labels, covariates = cov_names),
        class = "biomimic_selection")
}


#' Select the Top Genes of Chosen Modules
#'
#' @param modules A \code{biomimic_modules} object.
#' @param which Character or integer, the modules to use (one factor each).
#' @param k Integer, genes taken from the top of each module (default 15).
#' @param Y Numeric matrix, samples x features; genes absent from it are
#'   skipped, so the same modules can be applied to another cohort.
#' @param covariates Optional data frame of covariates (see
#'   \code{\link{as_selection}}).
#'
#' @return A \code{biomimic_selection}.
#'
#' @export
select_panels <- function(modules, which, k = 15L, Y, covariates = NULL) {
    stopifnot(inherits(modules, "biomimic_modules"))
    mm <- modules$modules[which]
    if (any(vapply(mm, is.null, logical(1)))) stop("Unknown module in 'which'.", call. = FALSE)
    panels <- lapply(mm, function(m) utils::head(intersect(m$genes, colnames(Y)), k))
    as_selection(Y, panels, covariates)
}


#' @export
print.biomimic_modules <- function(x, ...) {
    cat(sprintf("biomimic modules (%s): %d modules\n", x$method, length(x$modules)))
    for (m in x$modules) {
        cat(sprintf("  %-8s %4d genes: %s%s\n", m$name, length(m$genes),
                    paste(utils::head(m$genes, 6), collapse = ", "),
                    if (length(m$genes) > 6) ", ..." else ""))
    }
    invisible(x)
}


#' @noRd
.module_from_genes <- function(name, genes, Y) {
    Z <- scale(Y[, genes, drop = FALSE])
    sv <- svd(Z, nu = 1, nv = 1)
    sc <- sv$u[, 1] * sv$d[1]
    if (stats::cor(sc, rowMeans(Z)) < 0) sc <- -sc
    kme <- as.numeric(stats::cor(Z, sc))
    o <- order(kme, decreasing = TRUE)
    list(name = name, genes = genes[o], membership = kme[o], score = sc)
}

#' @noRd
.new_modules <- function(mods, method) {
    structure(list(modules = mods, method = method), class = "biomimic_modules")
}
