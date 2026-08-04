# --------------------------------------------------------------------------
# Logging System for biomimic
# Uses `logger` package when available, falls back to message()
# --------------------------------------------------------------------------

# Internal environment for log state (fallback mode)
.biomimic_env <- new.env(parent = emptyenv())
.biomimic_env$log_active <- FALSE
.biomimic_env$log_file <- NULL
.biomimic_env$log_con <- NULL
.biomimic_env$log_level <- "info"

# Level hierarchy
.LOG_LEVELS <- c(debug = 1L, info = 2L, warn = 3L, error = 4L)


#' Start Logging for biomimic
#'
#' Enables logging of all pipeline operations. If the \code{logger} package
#' is installed, it is used for full-featured logging (multiple appenders,
#' custom layouts, namespace isolation). Otherwise, a lightweight built-in
#' logger writes to file and/or console.
#'
#' @param file Character, path to log file. If \code{NULL}, logs to console
#'   only.
#' @param level Character, minimum log level: \code{"debug"}, \code{"info"}
#'   (default), \code{"warn"}, or \code{"error"}.
#' @param console Logical, also print to console when logging to file
#'   (default \code{TRUE}).
#'
#' @examples
#' \dontrun{
#' biomimic_log_start("analysis.log", level = "info")
#' fit <- biomimic(Y, X, K = 1, k = 15)
#' biomimic_log_stop()
#' }
#'
#' @export
biomimic_log_start <- function(file = NULL, level = "info",
                                console = TRUE) {
    level <- match.arg(level, names(.LOG_LEVELS))

    if (.has_logger()) {
        # Configure logger namespace
        logger::log_threshold(
            switch(level,
                debug = logger::DEBUG,
                info  = logger::INFO,
                warn  = logger::WARN,
                error = logger::ERROR),
            namespace = "biomimic"
        )

        if (!is.null(file)) {
            if (console) {
                logger::log_appender(
                    logger::appender_tee(file),
                    namespace = "biomimic"
                )
            } else {
                logger::log_appender(
                    logger::appender_file(file),
                    namespace = "biomimic"
                )
            }
        } else {
            logger::log_appender(
                logger::appender_stdout,
                namespace = "biomimic"
            )
        }

        logger::log_layout(
            logger::layout_glue_generator(
                format = "[{time}] {level} {msg}"
            ),
            namespace = "biomimic"
        )

        .log_info("Logging started (logger package) at level '%s'",
                  level)
    } else {
        # Fallback: built-in logger
        .biomimic_env$log_active <- TRUE
        .biomimic_env$log_level <- level

        if (!is.null(file)) {
            .biomimic_env$log_file <- file
            .biomimic_env$log_con <- file(file, open = "a")
            .biomimic_env$log_console <- console
        } else {
            .biomimic_env$log_file <- NULL
            .biomimic_env$log_con <- NULL
            .biomimic_env$log_console <- TRUE
        }

        .log_info("Logging started (built-in) at level '%s'", level)
    }

    invisible(NULL)
}


#' Stop Logging for biomimic
#'
#' Disables logging and closes any open log file connections.
#'
#' @export
biomimic_log_stop <- function() {
    if (.has_logger()) {
        .log_info("Logging stopped")
        logger::log_threshold(logger::WARN, namespace = "biomimic")
        logger::log_appender(logger::appender_stdout, namespace = "biomimic")
    } else {
        .log_info("Logging stopped")
        if (!is.null(.biomimic_env$log_con)) {
            close(.biomimic_env$log_con)
        }
        .biomimic_env$log_active <- FALSE
        .biomimic_env$log_file <- NULL
        .biomimic_env$log_con <- NULL
    }
    invisible(NULL)
}


# ===========================================================================
# Internal log functions — used throughout the package
# ===========================================================================

#' @noRd
.has_logger <- function() {
    requireNamespace("logger", quietly = TRUE)
}

# Pre-format with sprintf so messages render identically whether they go to
# the logger package (which would otherwise glue-interpret them) or the
# built-in fallback. Callers use sprintf-style format strings (%s, %d).
#' @noRd
.log_fmt <- function(msg, ...) {
    if (...length() > 0L) sprintf(msg, ...) else msg
}

#' @noRd
.log_debug <- function(msg, ...) {
    if (.has_logger()) {
        logger::log_debug(.log_fmt(msg, ...), namespace = "biomimic")
    } else {
        .builtin_log("DEBUG", msg, ...)
    }
}

#' @noRd
.log_info <- function(msg, ...) {
    if (.has_logger()) {
        logger::log_info(.log_fmt(msg, ...), namespace = "biomimic")
    } else {
        .builtin_log("INFO ", msg, ...)
    }
}

#' @noRd
.log_warn <- function(msg, ...) {
    if (.has_logger()) {
        logger::log_warn(.log_fmt(msg, ...), namespace = "biomimic")
    } else {
        .builtin_log("WARN ", msg, ...)
    }
}

#' @noRd
.log_error <- function(msg, ...) {
    if (.has_logger()) {
        logger::log_error(.log_fmt(msg, ...), namespace = "biomimic")
    } else {
        .builtin_log("ERROR", msg, ...)
    }
}


#' @noRd
.builtin_log <- function(level, msg, ...) {
    if (!.biomimic_env$log_active) return(invisible(NULL))

    # Check level threshold
    level_num <- .LOG_LEVELS[tolower(trimws(level))]
    threshold <- .LOG_LEVELS[.biomimic_env$log_level]
    if (is.na(level_num) || level_num < threshold) return(invisible(NULL))

    # Format message (use sprintf for built-in, since logger uses glue)
    formatted <- tryCatch(
        sprintf(msg, ...),
        error = function(e) msg
    )

    line <- sprintf("[%s] %s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    level, formatted)

    # Write to file
    if (!is.null(.biomimic_env$log_con)) {
        cat(line, file = .biomimic_env$log_con)
        flush(.biomimic_env$log_con)
    }

    # Write to console
    if (.biomimic_env$log_console) {
        message(trimws(line))
    }

    invisible(NULL)
}
