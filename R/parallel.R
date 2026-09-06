# --------------------------------------------------------------------------
# Parallel execution helpers for biomimic
# Embarrassingly-parallel work (bootstrap calibration) is distributed across
# cores with a PSOCK cluster (works on Windows, where fork/mclapply is
# unavailable). Mirrors the biocca helper of the same shape.
# --------------------------------------------------------------------------

#' Number of worker cores to use
#'
#' Honours the \code{BIOMIMIC_NCORES} environment variable; otherwise uses
#' \code{detectCores() - 1} (leaving one core for the OS). Under
#' \code{R CMD check --as-cran} the result is capped at 2, which CRAN requires
#' and \code{makeCluster()} enforces with an error above that.
#' @noRd
.biomimic_n_cores <- function() {
    env <- suppressWarnings(as.integer(Sys.getenv("BIOMIMIC_NCORES", "")))
    n <- if (!is.na(env) && env >= 1L) {
        env
    } else {
        nc <- tryCatch(parallel::detectCores(logical = TRUE),
                       error = function(e) 1L)
        if (is.na(nc) || nc < 1L) nc <- 1L
        max(1L, nc - 1L)
    }
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && !identical(tolower(chk), "false")) n <- min(n, 2L)
    n
}

# Records which build mismatches have already been reported, so the warning
# below fires once per session rather than once per parallel call.
.biomimic_worker_check_seen <- new.env(parent = emptyenv())

#' Abort if the workers would run a different biomimic than this session
#'
#' A function shipped to a PSOCK worker carries a *reference* to the biomimic
#' namespace, not a copy of it, so the worker resolves biomimic's internal
#' calls against whatever biomimic is installed there -- regardless of which
#' copy the master is running, and regardless of whether the source tree was
#' also sourced into the worker's global environment. A stale worker install
#' therefore runs old code for every replicate. That fails loudly when the old
#' version is missing an object, but silently when the difference is only
#' numerical, which would quietly contaminate bootstrap calibration and
#' simulations. Compare the two copies and refuse to proceed on a version
#' mismatch; warn when the versions agree but the builds differ -- either a
#' source tree reinstalled without a version bump, or a \code{load_all()}
#' session whose edits have not reached the workers.
#'
#' Skipped when biomimic is not installed on the workers (the source-tree
#' fallback in \code{.biomimic_parallel} covers that case) or when
#' \code{BIOMIMIC_SKIP_VERSION_CHECK} is set to a non-empty value.
#'
#' @param cl A running PSOCK cluster.
#' @noRd
.biomimic_check_worker_version <- function(cl) {
    if (nzchar(Sys.getenv("BIOMIMIC_SKIP_VERSION_CHECK")))
        return(invisible(NULL))

    # An uninstalled source tree (devtools::load_all) has no Built field, so
    # the [1L] indexing turns that into NA rather than a zero-length entry.
    stamp <- function(d) {
        if (!is.list(d)) return(NULL)
        c(as.character(d$Version)[1L], as.character(d$Built)[1L])
    }
    here <- stamp(tryCatch(utils::packageDescription("biomimic"),
                           error = function(e) NULL))
    if (is.null(here)) return(invisible(NULL))

    there <- tryCatch(parallel::clusterEvalQ(cl, {
        .d <- tryCatch(utils::packageDescription("biomimic"),
                       error = function(e) NULL)
        if (is.list(.d)) c(as.character(.d$Version)[1L],
                           as.character(.d$Built)[1L]) else NULL
    }), error = function(e) NULL)
    if (!length(there) || any(vapply(there, is.null, logical(1))))
        return(invisible(NULL))

    versions <- unique(vapply(there, `[`, character(1), 1L))
    if (!identical(versions, here[1L])) {
        stop(sprintf(paste0(
            "biomimic version mismatch: this session runs %s but the ",
            "parallel workers load %s. The workers resolve biomimic's ",
            "internal functions against their own installed copy, so they ",
            "would compute the results with different code. Reinstall ",
            "biomimic so the two agree (devtools::install(\"biomimic\")), or ",
            "set BIOMIMIC_NCORES=1 to run serially."),
            here[1L], paste(versions, collapse = ", ")), call. = FALSE)
    }

    builds <- unique(vapply(there, `[`, character(1), 2L))
    if (identical(builds, here[2L])) return(invisible(NULL))

    # Same version, different build. Warn once per distinct pairing: under
    # load_all() this holds for every parallel call in the session, and one
    # warning per bootstrap replicate set would be pure noise.
    key <- paste(c(here[2L], builds), collapse = "|")
    if (!is.null(.biomimic_worker_check_seen[[key]])) return(invisible(NULL))
    assign(key, TRUE, envir = .biomimic_worker_check_seen)

    msg <- if (is.na(here[2L])) {
        sprintf(paste0(
            "This session runs biomimic %s from an uninstalled source tree ",
            "(devtools::load_all), but the parallel workers run the installed ",
            "build (%s). Workers resolve biomimic's internal functions ",
            "through their installed copy, so your edits do NOT reach them. ",
            "Reinstall biomimic, or set BIOMIMIC_NCORES=1 to run serially."),
            here[1L], paste(builds, collapse = "; "))
    } else {
        sprintf(paste0(
            "The parallel workers load a different build of biomimic %s than ",
            "this session (workers: %s; here: %s). The version numbers agree, ",
            "so this is usually a source tree reinstalled without a version ",
            "bump -- the workers may still be running the older code."),
            here[1L], paste(builds, collapse = "; "), here[2L])
    }
    warning(paste(msg, "(Warned once per session.)"), call. = FALSE)
    invisible(NULL)
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
    # A PSOCK worker starts with its own default .libPaths() and does not
    # inherit the master's. That matters whenever biomimic lives somewhere the
    # default paths do not cover -- above all the temporary library that
    # R CMD build/check installs into, where the workers would otherwise fail
    # to find biomimic at all. Point them at the master's libraries so they
    # load the very copy this session is running.
    lib_paths <- .libPaths()
    parallel::clusterExport(cl, "lib_paths", envir = environment())
    parallel::clusterEvalQ(cl, .libPaths(lib_paths))
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
    # Both branches above leave the workers resolving biomimic internals
    # through their own installed copy, so this check applies to either one.
    .biomimic_check_worker_version(cl)
    if (length(export) > 0L) {
        parallel::clusterExport(cl, export, envir = envir)
    }
    parallel::parLapply(cl, x, fun)
}
