# --------------------------------------------------------------------------
# Publication-quality figure functions for the SEM-estimation figures
# (structural coefficients and the indirect/direct mediation decomposition),
# so the paper estimation figures are generatable from a fitted model.
# --------------------------------------------------------------------------

#' Structural-coefficient forest plot
#'
#' Forest plot of the MIMIC structural coefficients (the effect of each
#' covariate on the latent pathway) with 95% confidence intervals, coloured by
#' significance. Reproduces the structural panel of the paper's estimation
#' figure.
#'
#' @param fit A fitted [biomimic()] model (any MIMIC variant).
#' @param labels Optional named character vector mapping covariate names (the
#'   right-hand sides of the structural equation) to display labels.
#' @param sig_level Significance threshold for colouring (default 0.05).
#' @param title Plot title.
#'
#' @return A \pkg{ggplot} object.
#'
#' @examples
#' data(neuro_antibodyome)
#' \donttest{
#' fit <- biomimic(neuro_antibodyome$proteins, neuro_antibodyome$X, K = 1)
#' plot_structural_effects(fit,
#'     labels = c(group = "AD vs control", age = "Age", sex = "Sex"))
#' }
#' @seealso [plot_mediation()], [test_hypotheses()]
#' @export
plot_structural_effects <- function(fit, labels = NULL, sig_level = 0.05,
                                    title = "Structural coefficients") {
    stopifnot(inherits(fit, "biomimic_fit"))
    pe <- lavaan::parameterEstimates(fit$lavaan_fit)
    st <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]
    if (!nrow(st)) stop("No structural coefficients found in this fit.")
    term <- if (!is.null(labels) && all(st$rhs %in% names(labels)))
        unname(labels[st$rhs]) else st$rhs
    sig_lab <- sprintf("p < %g", sig_level)
    d <- data.frame(
        term = factor(term, levels = rev(unique(term))),
        estimate = st$est, lo = st$ci.lower, hi = st$ci.upper,
        sig = factor(ifelse(st$pvalue < sig_level, sig_lab, "n.s."),
                     levels = c(sig_lab, "n.s.")))
    pal <- stats::setNames(c("#E69F00", "grey60"), c(sig_lab, "n.s."))
    ggplot2::ggplot(d, ggplot2::aes(.data$estimate, .data$term,
                                    colour = .data$sig)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            colour = "grey55") +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                                height = 0.25, linewidth = 0.6) +
        ggplot2::geom_point(size = 2.8) +
        ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
        ggplot2::labs(title = title, y = NULL,
                      x = "Effect on latent pathway (95% CI)") +
        theme_biomimic() +
        ggplot2::theme(legend.position = "bottom")
}


#' Indirect-versus-direct mediation decomposition
#'
#' For a MIMIC-with-direct-effects model, decomposes the effect of a covariate
#' on each protein into the part mediated by the latent pathway (indirect) and
#' the residual direct effect, as paired bars. Reproduces the mediation panel
#' of the paper's estimation figure.
#'
#' @param fit A fitted [biomimic()] model with \code{model_type =
#'   "direct"} (so a mediation decomposition is available).
#' @param covariate Name of the covariate whose effect is decomposed
#'   (default \code{"group"}).
#' @param title Plot title.
#'
#' @return A \pkg{ggplot} object.
#'
#' @examples
#' data(neuro_antibodyome)
#' \donttest{
#' fit <- biomimic(neuro_antibodyome$proteins, neuro_antibodyome$X, K = 1,
#'                 model_type = "direct")
#' plot_mediation(fit)
#' }
#' @seealso [plot_structural_effects()], [test_hypotheses()]
#' @export
plot_mediation <- function(fit, covariate = "group",
                           title = "Group effect: indirect vs direct") {
    stopifnot(inherits(fit, "biomimic_fit"))
    med <- test_hypotheses(fit)$mediation
    if (is.null(med) || !nrow(med))
        stop("No mediation decomposition available; fit with ",
             "model_type = 'direct'.")
    if (!is.null(med$covariate))
        med <- med[med$covariate == covariate, , drop = FALSE]
    if (!nrow(med))
        stop(sprintf("No mediation rows for covariate '%s'.", covariate))
    med$protein <- sub("_1$", "", med$protein)
    dd <- data.frame(
        protein = factor(rep(med$protein, 2),
                         levels = rev(unique(med$protein))),
        type = factor(rep(c("Indirect (via pathway)", "Direct"),
                          each = nrow(med)),
                      levels = c("Indirect (via pathway)", "Direct")),
        effect = c(med$indirect, med$direct))
    pal <- c("Indirect (via pathway)" = "#0072B2", "Direct" = "#D55E00")
    ggplot2::ggplot(dd, ggplot2::aes(.data$effect, .data$protein,
                                     fill = .data$type)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
        ggplot2::geom_col(width = 0.65,
                          position = ggplot2::position_dodge(0.7)) +
        ggplot2::scale_fill_manual(values = pal, name = NULL) +
        ggplot2::labs(title = title, y = NULL,
                      x = "Effect on protein (case vs control)") +
        theme_biomimic() +
        ggplot2::theme(legend.position = "bottom")
}


#' Cross-condition recurrence of selected variables
#'
#' Tile plot of which selected variables recur across several conditions (e.g.
#' diseases), coloured by their selection rank, highlighting the shared core
#' versus condition-specific selections. Reproduces the cross-disease
#' recurrence figure of the paper.
#'
#' @param selections A named list, one element per condition. Each element is
#'   either a character vector of selected variable names (best first) or a
#'   data frame with a \code{variable} column and an optional \code{rank}
#'   column. The list names are the condition labels.
#' @param min_recurrence Keep variables selected in at least this many
#'   conditions (default 2).
#'
#' @return A \pkg{ggplot} object.
#'
#' @examples
#' sel <- list(
#'     AD      = c("SSR4", "CRP", "HESX1", "RAMP1", "TMED5"),
#'     PD      = c("CRP", "SSR4", "CLBA1", "HESX1", "QTRT1"),
#'     earlyPD = c("HESX1", "CRP", "SSR4", "TMED5", "ZDHHC7"),
#'     MS      = c("CLBA1", "QTRT1", "CRP", "SSR4", "MGC18216"))
#' plot_selection_recurrence(sel)
#' @export
plot_selection_recurrence <- function(selections, min_recurrence = 2L) {
    if (!is.list(selections) || is.null(names(selections)))
        stop("'selections' must be a named list, one element per condition.")
    long <- do.call(rbind, lapply(names(selections), function(k) {
        s <- selections[[k]]
        if (is.data.frame(s)) {
            v <- as.character(s$variable)
            r <- if (!is.null(s$rank)) s$rank else seq_along(v)
        } else {
            v <- as.character(s); r <- seq_along(v)
        }
        data.frame(variable = v, rank = r, condition = k,
                   stringsAsFactors = FALSE)
    }))
    cnt <- table(long$variable)
    recur <- names(cnt[cnt >= min_recurrence])
    if (!length(recur))
        stop("No variables recur in >= ", min_recurrence, " conditions.")
    n_specific <- sum(cnt == 1L)
    ord <- vapply(recur, function(p)
        c(cnt[[p]], min(long$rank[long$variable == p])), numeric(2))
    recur <- recur[order(-ord[1, ], ord[2, ])]
    conds <- names(selections)
    grid <- expand.grid(variable = recur, condition = conds,
                        stringsAsFactors = FALSE)
    grid <- merge(grid, long, by = c("variable", "condition"), all.x = TRUE)
    grid$variable <- factor(grid$variable, levels = rev(recur))
    grid$condition <- factor(grid$condition, levels = conds)
    ggplot2::ggplot(grid, ggplot2::aes(.data$condition, .data$variable,
                                       fill = .data$rank)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
        ggplot2::geom_text(ggplot2::aes(label = .data$rank), na.rm = TRUE,
                           size = 3, colour = "white", fontface = "bold") +
        ggplot2::scale_fill_gradient(low = "#08519C", high = "#9ECAE1",
                                     na.value = "grey93",
                                     name = "Selection\nrank (1 = top)") +
        ggplot2::labs(
            title = "Cross-condition recurrence of selections",
            subtitle = sprintf(paste("Variables selected in >= %d conditions",
                                     "(%d were condition-specific)."),
                               min_recurrence, n_specific),
            x = NULL, y = NULL) +
        theme_biomimic() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(face = "bold"),
                       plot.subtitle = ggplot2::element_text(size = 9,
                                                             colour = "grey35"))
}


#' Forest plot of group-effect estimates across conditions
#'
#' Forest plot of one structural group-effect estimate per condition (e.g. one
#' per disease), with optional 95% confidence intervals, coloured by
#' significance. Reproduces the cross-disease latent-effect (disease
#' progression) figure of the paper.
#'
#' @param effects A data frame with columns \code{label} (character),
#'   \code{estimate} (numeric) and \code{p} (numeric p-value); optionally
#'   \code{ci_lower} and \code{ci_upper} for confidence intervals. Rows are
#'   drawn top-to-bottom in the given order.
#' @param sig_level Significance threshold for colouring (default 0.05).
#' @param title,subtitle Plot title and optional subtitle.
#' @param xlab x-axis label (default \eqn{\hat{B}_{group}}).
#'
#' @return A \pkg{ggplot} object.
#'
#' @examples
#' effects <- data.frame(
#'     label    = c("early PD", "PD", "AD", "MS"),
#'     estimate = c(1.01, -0.50, -0.85, -0.20),
#'     ci_lower = c(0.55, -1.10, -1.23, -0.70),
#'     ci_upper = c(1.47, 0.10, -0.47, 0.30),
#'     p        = c(0.001, 0.08, 0.001, 0.40))
#' plot_group_effects(effects,
#'     subtitle = "Latent disease effect across comparisons")
#' @export
plot_group_effects <- function(effects, sig_level = 0.05,
                               title = "Group effect across conditions",
                               subtitle = NULL,
                               xlab = expression(hat(B)[group])) {
    if (!all(c("label", "estimate", "p") %in% names(effects)))
        stop("'effects' must have columns 'label', 'estimate' and 'p'.")
    e <- effects
    e$label <- factor(as.character(e$label),
                      levels = rev(as.character(e$label)))
    sig_lab <- sprintf("p < %g", sig_level)
    e$sig <- factor(ifelse(e$p < sig_level, sig_lab, "n.s."),
                    levels = c(sig_lab, "n.s."))
    pal <- stats::setNames(c("#B2182B", "grey60"), c(sig_lab, "n.s."))
    has_ci <- all(c("ci_lower", "ci_upper") %in% names(e))
    p <- ggplot2::ggplot(e, ggplot2::aes(.data$estimate, .data$label,
                                         colour = .data$sig)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            colour = "grey40")
    if (has_ci)
        p <- p + ggplot2::geom_errorbarh(
            ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
            height = 0.2, linewidth = 0.9)
    p +
        ggplot2::geom_point(size = 4) +
        ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
        ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
        theme_biomimic() +
        ggplot2::theme(legend.position = "top",
                       plot.subtitle = ggplot2::element_text(colour = "grey40",
                                                             size = 9))
}
