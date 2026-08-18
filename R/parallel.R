# --------------------------------------------------------------------------
# Parallel execution helpers for biomimic
# Embarrassingly-parallel work (bootstrap calibration) is distributed across
# cores with a PSOCK cluster (works on Windows, where fork/mclapply is
# unavailable). Mirrors the biocca helper of the same shape.
# --------------------------------------------------------------------------

#' Number of worker cores to use
#'
#' Honours the \code{BIOMIMIC_NCORES} environment variable; otherwise uses
#' \code{detectCores() - 1} (leaving one core for the OS).
#' @noRd
.biomimic_n_cores <- function() {
    env <- suppressWarnings(as.integer(Sys.getenv("BIOMIMIC_NCORES", "")))
    if (!is.na(env) && env >= 1L) return(env)
    nc <- tryCatch(parallel::detectCores(logical = TRUE),
                   error = function(e) 1L)
    if (is.na(nc) || nc < 1L) nc <- 1L
    max(1L, nc - 1L)
}

#' Apply a function over a list/vector, in parallel over a PSOCK cluster
#'
#' Spins up a PSOCK cluster, makes biomimic available on each worker
#' (installed package preferred; source-tree fallback for development),
#' exports the named \code{export} objects, and runs \code{parLapply}.
#' Falls back to a plain \code{lapply} when only one core is requested,
#' so callers need no separate serial branch.
#'
#' Reproducibility: no parallel RNG stream is set; \code{fun} must seed
#' itself deterministically per element (e.g. \code{set.seed(seed + b)}),
#' which makes the parallel output identical to the serial one.
#'
#' @param x List/vector of items to map over.
#' @param fun Function applied to each element.
#' @param export Character vector of object names to export to the workers.
#' @param envir Environment in which \code{export} names are found.
#' @param n_cores Number of workers (default \code{.biomimic_n_cores()}).
#' @noRd
.biomimic_parallel <- function(x, fun, export = character(0),
                               envir = parent.frame(), n_cores = NULL) {
    if (is.null(n_cores)) n_cores <- .biomimic_n_cores()
    has_par <- requireNamespace("parallel", quietly = TRUE)
    if (n_cores <= 1L || !has_par) return(lapply(x, fun))

    pkg_dir <- if (dir.exists("biomimic/R")) "biomimic/R"
               else if (dir.exists("R")) "R" else NULL
    pkg_files <- if (!is.null(pkg_dir)) {
        normalizePath(list.files(pkg_dir, pattern = "\\.R$",
                                 full.names = TRUE))
    } else character(0)

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    have_pkg <- all(unlist(parallel::clusterEvalQ(cl,
        requireNamespace("biomimic", quietly = TRUE))))
    if (have_pkg) {
        parallel::clusterEvalQ(cl, suppressMessages(library(biomimic)))
    } else if (length(pkg_files)) {
        parallel::clusterExport(cl, "pkg_files", envir = environment())
        parallel::clusterEvalQ(cl,
            for (.f in pkg_files) source(.f, local = FALSE))
    } else {
        stop("biomimic is not installed and no source tree was found ",
             "for the workers.", call. = FALSE)
    }
    if (length(export) > 0L) {
        parallel::clusterExport(cl, export, envir = envir)
    }
    parallel::parLapply(cl, x, fun)
}
