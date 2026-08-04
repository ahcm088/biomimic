# --------------------------------------------------------------------------
# Unified Plotting System for biomimic
# All plots use ggplot2 and return ggplot objects (customisable with +)
# --------------------------------------------------------------------------

# ===========================================================================
# Custom theme and palette
# ===========================================================================

#' biomimic ggplot2 Theme
#'
#' A publication-quality theme with white background and light grey grid.
#' Designed for bioinformatics figures in journals like BMC Bioinformatics.
#'
#' @param base_size Base font size (default 12).
#' @return A ggplot2 theme object.
#' @export
theme_biomimic <- function(base_size = 12) {
    ggplot2::theme_minimal(base_size = base_size) %+replace%
        ggplot2::theme(
            panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
            panel.grid.major  = ggplot2::element_line(colour = "grey92",
                                                       linewidth = 0.3),
            panel.grid.minor  = ggplot2::element_line(colour = "grey96",
                                                       linewidth = 0.15),
            plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
            plot.title        = ggplot2::element_text(size = base_size * 1.2,
                                                       face = "bold",
                                                       hjust = 0,
                                                       margin = ggplot2::margin(
                                                           b = base_size * 0.5)),
            plot.subtitle     = ggplot2::element_text(size = base_size * 0.9,
                                                       colour = "grey40",
                                                       hjust = 0),
            axis.title        = ggplot2::element_text(size = base_size * 0.95),
            axis.text         = ggplot2::element_text(size = base_size * 0.85),
            legend.background = ggplot2::element_rect(fill = "white", colour = NA),
            legend.key        = ggplot2::element_rect(fill = "white", colour = NA),
            strip.background  = ggplot2::element_rect(fill = "grey95", colour = NA),
            strip.text        = ggplot2::element_text(face = "bold",
                                                       size = base_size * 0.9)
        )
}

#' biomimic Colour Palette
#'
#' A qualitative colour palette with 15 distinct colours optimised for
#' colour-blind safety. Returns a character vector of hex codes.
#'
#' @param n Number of colours to return (default 15, max 20).
#' @return Character vector of hex colour codes.
#' @export
palette_biomimic <- function(n = 15L) {
    # 20-colour palette: Okabe-Ito extended + carefully chosen supplements
    full <- c(
        "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7", "#999999", "#882255", "#44AA99",
        "#332288", "#DDCC77", "#117733", "#88CCEE", "#AA4499",
        "#661100", "#6699CC", "#44BB99", "#EE6677", "#228833"
    )
    if (n > length(full)) n <- length(full)
    full[seq_len(n)]
}


# ===========================================================================
# Main plot dispatcher — S3 method
# ===========================================================================

#' Plot Method for biomimic_fit Objects
#'
#' Unified plotting interface. All plots return a \code{ggplot} object that
#' can be further customised with \code{+} (e.g., \code{+ labs(title = "...")}).
#'
#' @param x A \code{biomimic_fit} object.
#' @param type Character, the type of plot. See Details.
#' @param ... Additional arguments passed to the specific plot function.
#'   Common arguments: \code{group_var}, \code{factor} (which factor),
#'   \code{threshold}, \code{standardized}.
#'
#' @section Available plot types:
#' \describe{
#'   \item{\code{"loadings"}}{Factor loading bar plot (all models)}
#'   \item{\code{"volcano"}}{Loading magnitude vs. significance (all)}
#'   \item{\code{"rsquared"}}{R-squared per indicator (all)}
#'   \item{\code{"scores"}}{Factor scores by group — density (all)}
#'   \item{\code{"scores_box"}}{Factor scores by group — box/violin (all)}
#'   \item{\code{"scores_scatter"}}{Score scatter, factor 1 vs 2 (K>=2)}
#'   \item{\code{"residuals"}}{Residual correlation heatmap (all)}
#'   \item{\code{"qq"}}{QQ plot of residuals (all)}
#'   \item{\code{"mahalanobis"}}{Outlier detection (all)}
#'   \item{\code{"gof"}}{Goodness-of-fit summary (all)}
#'   \item{\code{"correlation"}}{Observed correlation heatmap (all)}
#'   \item{\code{"network"}}{Correlation network graph (all, needs igraph)}
#'   \item{\code{"forest"}}{Structural coefficient forest plot (mimic, multi_mimic, direct)}
#'   \item{\code{"loadings_heatmap"}}{Loading matrix heatmap (multi_mimic)}
#'   \item{\code{"factor_cor"}}{Factor correlation heatmap (multi_mimic)}
#'   \item{\code{"biplot"}}{Score biplot with loading arrows (multi_mimic)}
#'   \item{\code{"invariance"}}{Fit index trajectory (multigroup)}
#'   \item{\code{"loadings_group"}}{Per-group loading comparison (multigroup)}
#'   \item{\code{"delta_cfi"}}{Delta-CFI bar plot (multigroup)}
#'   \item{\code{"decomposition"}}{Indirect vs direct effect (direct)}
#'   \item{\code{"screening"}}{Stage 1 SNR ranking (all)}
#'   \item{\code{"scree"}}{Scree plot for K selection (all)}
#' }
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' fit <- fit_sem(sel, type = "mimic")
#' plot(fit, "loadings")
#' plot(fit, "scores", group_var = "group") +
#'     ggplot2::scale_fill_brewer(palette = "Set1")  # user customisation
#' }
#'
#' @export
plot.biomimic_fit <- function(x, type = "loadings", ...) {
    .check_ggplot2()

    type <- match.arg(type, choices = c(
        # Universal
        "loadings", "volcano", "rsquared", "scores", "scores_box",
        "scores_scatter", "residuals", "qq", "mahalanobis", "gof",
        "correlation", "network",
        # MIMIC / structural
        "forest",
        # Multi-factor
        "loadings_heatmap", "factor_cor", "biplot",
        # Multi-group
        "invariance", "loadings_group", "delta_cfi",
        # Direct effects
        "decomposition",
        # Stage 1-2
        "screening", "scree"
    ))

    switch(type,
        # Universal
        loadings         = .plot_loadings(x, ...),
        volcano          = .plot_volcano(x, ...),
        rsquared         = .plot_rsquared(x, ...),
        scores           = .plot_scores(x, plot_type = "density", ...),
        scores_box       = .plot_scores(x, plot_type = "violin", ...),
        scores_scatter   = .plot_scores_scatter(x, ...),
        residuals        = .plot_residuals(x, ...),
        qq               = .plot_qq(x, ...),
        mahalanobis      = .plot_mahalanobis(x, ...),
        gof              = .plot_gof(x, ...),
        correlation      = .plot_correlation(x, ...),
        network          = .plot_network(x, ...),
        # Structural
        forest           = .plot_forest(x, ...),
        # Multi-factor
        loadings_heatmap = .plot_loadings_heatmap(x, ...),
        factor_cor       = .plot_factor_cor(x, ...),
        biplot           = .plot_biplot(x, ...),
        # Multi-group
        invariance       = .plot_invariance(x, ...),
        loadings_group   = .plot_loadings_group(x, ...),
        delta_cfi        = .plot_delta_cfi(x, ...),
        # Direct
        decomposition    = .plot_decomposition(x, ...),
        # Stage 1-2
        screening        = .plot_screening(x, ...),
        scree            = .plot_scree(x, ...)
    )
}

#' @export
plot.biomimic_selection <- function(x, type = "screening", ...) {
    .check_ggplot2()
    type <- match.arg(type, c("screening", "scree"))
    switch(type,
        screening = .plot_screening_sel(x, ...),
        scree     = .plot_scree_sel(x, ...)
    )
}


# ===========================================================================
# Internal helper
# ===========================================================================

#' @noRd
.check_ggplot2 <- function() {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. Install with:\n",
             "  install.packages('ggplot2')")
    }
}

#' @noRd
.pal <- function(n = 15) palette_biomimic(n)


# ===========================================================================
# UNIVERSAL PLOTS
# ===========================================================================

#' @noRd
.plot_loadings <- function(fit, factor = 1L, standardized = TRUE,
                            threshold = 0.3, ...) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, standardized = TRUE)
    fname <- paste0("F", factor)
    ld <- pe[pe$op == "=~" & pe$lhs == fname, , drop = FALSE]

    est_col <- if (standardized && "std.all" %in% names(ld)) "std.all" else "est"

    df <- data.frame(
        variable = ld$rhs,
        loading  = ld[[est_col]],
        label    = .bio_label(ld$rhs, fit$indicator_labels),
        stringsAsFactors = FALSE
    )
    df$label <- stats::reorder(df$label, df$loading)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$loading, y = .data$label,
                                      fill = .data$loading > 0)) +
        ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
        ggplot2::geom_vline(xintercept = c(-threshold, threshold),
                            linetype = "dashed", colour = "grey60", linewidth = 0.3) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
        ggplot2::scale_fill_manual(values = c("TRUE" = .pal(2)[2],
                                               "FALSE" = .pal(2)[1])) +
        ggplot2::labs(
            x = if (standardized) "Standardized loading" else "Loading",
            y = NULL,
            title = paste("Factor Loadings:", fit$factor_labels[factor]),
            subtitle = paste0("Dashed lines at \u00b1", threshold)
        ) +
        theme_biomimic()
}


#' @noRd
.plot_volcano <- function(fit, factor = 1L, p_threshold = 0.05, ...) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, standardized = TRUE)
    fname <- paste0("F", factor)
    ld <- pe[pe$op == "=~" & pe$lhs == fname, , drop = FALSE]

    est_col <- if ("std.all" %in% names(ld)) "std.all" else "est"
    df <- data.frame(
        variable  = ld$rhs,
        loading   = ld[[est_col]],
        neglog10p = -log10(pmax(ld$pvalue, 1e-300)),
        pvalue    = ld$pvalue,
        label     = .bio_label(ld$rhs, fit$indicator_labels),
        significant = ld$pvalue < p_threshold,
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$loading, y = .data$neglog10p,
                                      colour = .data$significant)) +
        ggplot2::geom_point(size = 2.5, alpha = 0.8) +
        ggplot2::geom_hline(yintercept = -log10(p_threshold),
                            linetype = "dashed", colour = "grey50") +
        ggplot2::geom_text(data = df[df$significant, ],
                           ggplot2::aes(label = .data$label),
                           nudge_y = 0.15, size = 3, show.legend = FALSE) +
        ggplot2::scale_colour_manual(
            values = c("FALSE" = "grey70", "TRUE" = .pal(1)),
            labels = c("n.s.", paste0("p < ", p_threshold)),
            name = NULL) +
        ggplot2::labs(x = "Standardized loading", y = expression(-log[10](p)),
                      title = "Loading Volcano Plot",
                      subtitle = fit$factor_labels[factor]) +
        theme_biomimic()
}


#' @noRd
.plot_rsquared <- function(fit, ...) {
    r2 <- lavaan::lavInspect(fit$lavaan_fit, "rsquare")
    if (length(r2) == 0) stop("R-squared not available for this model.")

    df <- data.frame(
        variable = names(r2),
        rsq      = unname(r2),
        label    = .bio_label(names(r2), fit$indicator_labels),
        stringsAsFactors = FALSE
    )
    df$label <- stats::reorder(df$label, df$rsq)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$rsq, y = .data$label)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$rsq,
                                            yend = .data$label),
                              colour = "grey70", linewidth = 0.4) +
        ggplot2::geom_point(size = 3, colour = .pal(1)) +
        ggplot2::scale_x_continuous(limits = c(0, 1),
                                    labels = scales::percent_format(1)) +
        ggplot2::labs(x = "Variance explained (R\u00b2)", y = NULL,
                      title = "Variance Explained per Protein",
                      subtitle = "Proportion of protein variance captured by the latent model") +
        theme_biomimic()
}


#' @noRd
.plot_scores <- function(fit, group_var = NULL, factor = 1L,
                          plot_type = "density", group_labels = NULL, ...) {
    if (is.null(group_var)) {
        stop("group_var is required for score plots.")
    }
    scores_df <- predict_scores(fit, append_data = TRUE)
    score_col <- make.names(fit$factor_labels[factor])
    scores_df$score <- scores_df[[score_col]]
    scores_df$group <- as.factor(scores_df[[group_var]])
    if (!is.null(group_labels)) {
        levels(scores_df$group) <- group_labels[levels(scores_df$group)]
    }

    if (plot_type == "density") {
        ggplot2::ggplot(scores_df, ggplot2::aes(x = .data$score,
                                                  fill = .data$group)) +
            ggplot2::geom_density(alpha = 0.45, colour = NA) +
            ggplot2::geom_rug(ggplot2::aes(colour = .data$group),
                              alpha = 0.3, show.legend = FALSE) +
            ggplot2::scale_fill_manual(values = .pal(15), name = "Group") +
            ggplot2::scale_colour_manual(values = .pal(15)) +
            ggplot2::labs(x = paste("Estimated", fit$factor_labels[factor]),
                          y = "Density",
                          title = paste("Latent Pathway Activity:",
                                        fit$factor_labels[factor]),
                          subtitle = "Distribution of estimated factor scores by group") +
            theme_biomimic()
    } else {
        ggplot2::ggplot(scores_df, ggplot2::aes(x = .data$group,
                                                  y = .data$score,
                                                  fill = .data$group)) +
            ggplot2::geom_violin(alpha = 0.5, colour = NA) +
            ggplot2::geom_boxplot(width = 0.15, fill = "white",
                                  alpha = 0.85, outlier.alpha = 0.3) +
            ggplot2::scale_fill_manual(values = .pal(15), name = "Group") +
            ggplot2::labs(y = paste("Estimated", fit$factor_labels[factor]),
                          x = NULL,
                          title = paste("Latent Pathway Activity:",
                                        fit$factor_labels[factor])) +
            theme_biomimic()
    }
}


#' @noRd
.plot_scores_scatter <- function(fit, group_var = NULL, factor1 = 1L,
                                  factor2 = 2L, group_labels = NULL, ...) {
    if (fit$K < 2) stop("scores_scatter requires K >= 2.")
    if (is.null(group_var)) stop("group_var is required.")

    scores_df <- predict_scores(fit, append_data = TRUE)
    sc1 <- make.names(fit$factor_labels[factor1])
    sc2 <- make.names(fit$factor_labels[factor2])
    scores_df$f1 <- scores_df[[sc1]]
    scores_df$f2 <- scores_df[[sc2]]
    scores_df$group <- as.factor(scores_df[[group_var]])
    if (!is.null(group_labels)) {
        levels(scores_df$group) <- group_labels[levels(scores_df$group)]
    }

    ggplot2::ggplot(scores_df, ggplot2::aes(x = .data$f1, y = .data$f2,
                                             colour = .data$group)) +
        ggplot2::geom_point(alpha = 0.5, size = 1.5) +
        ggplot2::stat_ellipse(level = 0.68, linewidth = 0.6,
                               linetype = "solid") +
        ggplot2::stat_ellipse(level = 0.95, linewidth = 0.4,
                               linetype = "dashed") +
        ggplot2::scale_colour_manual(values = .pal(15), name = "Group") +
        ggplot2::labs(x = fit$factor_labels[factor1],
                      y = fit$factor_labels[factor2],
                      title = "Latent Score Scatter",
                      subtitle = "Ellipses: 68% and 95% concentration") +
        theme_biomimic()
}


#' @noRd
.plot_residuals <- function(fit, ...) {
    ra <- residual_analysis(fit)
    mat <- ra$residuals
    vnames <- rownames(mat)

    df <- expand.grid(var1 = vnames, var2 = vnames, stringsAsFactors = FALSE)
    df$residual <- as.vector(mat)
    df$label1 <- .bio_label(df$var1, fit$indicator_labels)
    df$label2 <- .bio_label(df$var2, fit$indicator_labels)
    df$label1 <- factor(df$label1, levels = .bio_label(vnames, fit$indicator_labels))
    df$label2 <- factor(df$label2,
                         levels = rev(.bio_label(vnames, fit$indicator_labels)))

    ggplot2::ggplot(df, ggplot2::aes(x = .data$label1, y = .data$label2,
                                      fill = .data$residual)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                      high = "#B2182B", midpoint = 0,
                                      name = "Residual\ncorrelation") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Standardized Residual Correlations",
                      subtitle = "Large values (|r| > 0.1) indicate local misfit") +
        theme_biomimic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' @noRd
.plot_qq <- function(fit, ...) {
    po <- predict_observed(fit)
    resid_mat <- po$residuals
    vars <- colnames(resid_mat)

    # Standardize residuals per variable
    df_list <- lapply(seq_along(vars), function(j) {
        r <- resid_mat[, j]
        r_std <- (r - mean(r)) / (stats::sd(r) + 1e-10)
        data.frame(
            variable = vars[j],
            label = .bio_label(vars[j], fit$indicator_labels),
            sample = sort(r_std),
            theoretical = stats::qnorm(stats::ppoints(length(r_std))),
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, df_list)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                              colour = "grey50") +
        ggplot2::geom_point(alpha = 0.3, size = 0.8, colour = .pal(1)) +
        ggplot2::facet_wrap(~ .data$label, scales = "free_y") +
        ggplot2::labs(x = "Theoretical quantiles", y = "Standardized residuals",
                      title = "QQ Plot of Model Residuals",
                      subtitle = "Departure from the diagonal indicates non-normality") +
        theme_biomimic()
}


#' @noRd
.plot_mahalanobis <- function(fit, p_threshold = 0.001, ...) {
    md <- mahalanobis_distance(fit, p_threshold = p_threshold)
    p <- length(fit$selection$variables)
    chi2_crit <- stats::qchisq(1 - p_threshold, df = p)

    ggplot2::ggplot(md, ggplot2::aes(x = .data$observation,
                                      y = .data$mahal_dist)) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$outlier),
                            alpha = 0.5, size = 1.5) +
        ggplot2::geom_hline(yintercept = chi2_crit, linetype = "dashed",
                            colour = "#CC3311", linewidth = 0.4) +
        ggplot2::geom_text(data = md[md$outlier, ],
                           ggplot2::aes(label = .data$observation),
                           nudge_y = max(md$mahal_dist) * 0.03,
                           size = 2.5, colour = "#CC3311") +
        ggplot2::scale_colour_manual(
            values = c("FALSE" = .pal(2)[2], "TRUE" = "#CC3311"),
            labels = c("Normal", "Outlier"), name = NULL) +
        ggplot2::labs(x = "Observation", y = "Mahalanobis distance",
                      title = "Multivariate Outlier Detection",
                      subtitle = sprintf("Threshold: \u03c7\u00b2(%d, p = %g)",
                                         p, p_threshold)) +
        theme_biomimic()
}


#' @noRd
.plot_gof <- function(fit, ...) {
    gof <- goodness_of_fit(fit)
    metrics <- data.frame(
        metric    = c("RMSEA", "CFI", "TLI", "SRMR"),
        value     = c(gof$rmsea, gof$cfi, gof$tli, gof$srmr),
        threshold = c(0.08, 0.90, 0.90, 0.08),
        excellent = c(0.05, 0.95, 0.95, 0.05),
        direction = c("lower", "higher", "higher", "lower"),
        stringsAsFactors = FALSE
    )
    metrics$metric <- factor(metrics$metric, levels = metrics$metric)
    metrics$good <- ifelse(metrics$direction == "lower",
                            metrics$value <= metrics$threshold,
                            metrics$value >= metrics$threshold)

    ggplot2::ggplot(metrics, ggplot2::aes(x = .data$metric, y = .data$value,
                                           fill = .data$good)) +
        ggplot2::geom_col(width = 0.6, alpha = 0.85) +
        ggplot2::geom_point(ggplot2::aes(y = .data$threshold),
                            shape = 4, size = 3, colour = "grey30") +
        ggplot2::geom_point(ggplot2::aes(y = .data$excellent),
                            shape = 1, size = 3, colour = "grey30") +
        ggplot2::scale_fill_manual(
            values = c("TRUE" = "#009E73", "FALSE" = "#D55E00"),
            labels = c("Below threshold", "Meets threshold"),
            name = NULL) +
        ggplot2::labs(x = NULL, y = "Value",
                      title = "Goodness-of-Fit Summary",
                      subtitle = "x = acceptable threshold, o = excellent threshold") +
        ggplot2::coord_flip() +
        theme_biomimic()
}


#' @noRd
.plot_correlation <- function(fit, method = "pearson", ...) {
    vars <- fit$selection$variables
    dat <- as.matrix(fit$selection$data[, vars, drop = FALSE])
    cor_mat <- stats::cor(dat, method = method)
    labels <- .bio_label(vars, fit$indicator_labels)
    rownames(cor_mat) <- colnames(cor_mat) <- labels

    # Hierarchical clustering for ordering
    hc <- stats::hclust(stats::as.dist(1 - abs(cor_mat)), method = "ward.D2")
    ord <- hc$order
    labels_ord <- labels[ord]

    df <- expand.grid(var1 = labels, var2 = labels, stringsAsFactors = FALSE)
    df$cor <- as.vector(cor_mat)
    df$var1 <- factor(df$var1, levels = labels_ord)
    df$var2 <- factor(df$var2, levels = rev(labels_ord))

    ggplot2::ggplot(df, ggplot2::aes(x = .data$var1, y = .data$var2,
                                      fill = .data$cor)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                      high = "#B2182B", midpoint = 0,
                                      limits = c(-1, 1),
                                      name = method) +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Observed Correlation Matrix",
                      subtitle = "Hierarchically clustered (Ward.D2)") +
        theme_biomimic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' @noRd
.plot_network <- function(fit, threshold = 0.3, method = "pearson",
                           layout = "fr", ...) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for network plots.\n",
             "  install.packages('igraph')")
    }
    if (!requireNamespace("ggraph", quietly = TRUE)) {
        stop("Package 'ggraph' is required for network plots.\n",
             "  install.packages('ggraph')")
    }

    vars <- fit$selection$variables
    dat <- as.matrix(fit$selection$data[, vars, drop = FALSE])
    cor_mat <- stats::cor(dat, method = method)
    labels <- .bio_label(vars, fit$indicator_labels)
    rownames(cor_mat) <- colnames(cor_mat) <- labels

    # Build adjacency
    adj <- cor_mat
    adj[abs(adj) < threshold] <- 0
    diag(adj) <- 0

    g <- igraph::graph_from_adjacency_matrix(abs(adj), mode = "undirected",
                                              weighted = TRUE)

    # Edge sign
    edges <- igraph::as_data_frame(g, what = "edges")
    if (nrow(edges) > 0) {
        for (i in seq_len(nrow(edges))) {
            r <- cor_mat[edges$from[i], edges$to[i]]
            edges$sign[i] <- ifelse(r >= 0, "positive", "negative")
        }
        igraph::E(g)$sign <- edges$sign
    }

    # Remove isolated nodes
    isolated <- which(igraph::degree(g) == 0)
    if (length(isolated) > 0 && length(isolated) < igraph::vcount(g)) {
        g <- igraph::delete_vertices(g, isolated)
    }

    ggraph::ggraph(g, layout = layout) +
        ggraph::geom_edge_link(
            ggplot2::aes(width = .data$weight, colour = .data$sign),
            alpha = 0.5) +
        ggraph::scale_edge_colour_manual(
            values = c("positive" = .pal(2)[2], "negative" = "#CC3311"),
            name = "Correlation") +
        ggraph::scale_edge_width(range = c(0.3, 2), guide = "none") +
        ggraph::geom_node_point(size = 5, colour = .pal(1), alpha = 0.9) +
        ggraph::geom_node_text(ggplot2::aes(label = .data$name),
                                repel = TRUE, size = 3) +
        ggplot2::labs(title = "Protein Correlation Network",
                      subtitle = sprintf("Edges: |r| > %.2f (%s)",
                                         threshold, method)) +
        ggraph::theme_graph(base_family = "") +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white",
                                                                colour = NA))
}


# ===========================================================================
# MIMIC / STRUCTURAL PLOTS
# ===========================================================================

#' @noRd
.plot_forest <- function(fit, ...) {
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, ci = TRUE, level = 0.95)
    struct <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), , drop = FALSE]
    if (nrow(struct) == 0) stop("No structural paths found.")

    struct$factor_label <- .bio_label(struct$lhs, .factor_label_map(fit))
    struct$predictor <- struct$rhs
    struct$significant <- struct$pvalue < 0.05
    struct$label <- paste(struct$predictor, "->", struct$factor_label)
    struct$label <- factor(struct$label, levels = rev(struct$label))

    ggplot2::ggplot(struct, ggplot2::aes(x = .data$est, y = .data$label,
                                          colour = .data$significant)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "solid",
                            colour = "grey70", linewidth = 0.3) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$ci.lower,
                                              xmax = .data$ci.upper),
                                height = 0.2, linewidth = 0.5) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_colour_manual(
            values = c("FALSE" = "grey60", "TRUE" = .pal(1)),
            labels = c("n.s.", "p < 0.05"), name = NULL) +
        ggplot2::labs(x = "Estimate (95% CI)", y = NULL,
                      title = "Structural Coefficients",
                      subtitle = "Effect of covariates on latent pathway(s)") +
        theme_biomimic()
}


# ===========================================================================
# MULTI-FACTOR PLOTS
# ===========================================================================

#' @noRd
.plot_loadings_heatmap <- function(fit, standardized = TRUE, ...) {
    if (fit$K < 2) stop("loadings_heatmap requires K >= 2.")

    pe <- lavaan::parameterEstimates(fit$lavaan_fit, standardized = TRUE)
    ld <- pe[pe$op == "=~", , drop = FALSE]
    est_col <- if (standardized && "std.all" %in% names(ld)) "std.all" else "est"

    df <- data.frame(
        indicator = .bio_label(ld$rhs, fit$indicator_labels),
        factor    = .bio_label(ld$lhs, .factor_label_map(fit)),
        loading   = ld[[est_col]],
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$factor, y = .data$indicator,
                                      fill = .data$loading)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$loading)),
                           size = 3) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                      high = "#B2182B", midpoint = 0,
                                      name = "Loading") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Factor Loading Matrix",
                      subtitle = "Block structure from VB-ARD Stage 2") +
        theme_biomimic()
}


#' @noRd
.plot_factor_cor <- function(fit, ...) {
    if (fit$K < 2) stop("factor_cor requires K >= 2.")

    # Extract factor (co)variances
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, standardized = TRUE)
    fvar <- pe[pe$op == "~~" & grepl("^F[0-9]+$", pe$lhs) &
               grepl("^F[0-9]+$", pe$rhs), , drop = FALSE]

    est_col <- if ("std.all" %in% names(fvar)) "std.all" else "est"
    K <- fit$K
    fnames <- fit$factor_labels

    mat <- matrix(0, K, K, dimnames = list(fnames, fnames))
    for (i in seq_len(nrow(fvar))) {
        r <- as.integer(sub("F", "", fvar$lhs[i]))
        c_idx <- as.integer(sub("F", "", fvar$rhs[i]))
        if (r <= K && c_idx <= K) {
            mat[r, c_idx] <- fvar[[est_col]][i]
            mat[c_idx, r] <- fvar[[est_col]][i]
        }
    }

    df <- expand.grid(f1 = fnames, f2 = fnames, stringsAsFactors = FALSE)
    df$cor <- as.vector(mat)
    df$f2 <- factor(df$f2, levels = rev(fnames))

    ggplot2::ggplot(df, ggplot2::aes(x = .data$f1, y = .data$f2,
                                      fill = .data$cor)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$cor)),
                           size = 4) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                      high = "#B2182B", midpoint = 0,
                                      limits = c(-1, 1), name = "Correlation") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Factor Correlations") +
        theme_biomimic()
}


#' @noRd
.plot_biplot <- function(fit, group_var = NULL, factor1 = 1L, factor2 = 2L,
                          group_labels = NULL, ...) {
    if (fit$K < 2) stop("biplot requires K >= 2.")

    # Scores
    scores_df <- predict_scores(fit, append_data = TRUE)
    sc1 <- make.names(fit$factor_labels[factor1])
    sc2 <- make.names(fit$factor_labels[factor2])
    scores_df$f1 <- scores_df[[sc1]]
    scores_df$f2 <- scores_df[[sc2]]

    # Loadings (for arrows)
    pe <- lavaan::parameterEstimates(fit$lavaan_fit, standardized = TRUE)
    ld <- pe[pe$op == "=~", , drop = FALSE]
    est_col <- if ("std.all" %in% names(ld)) "std.all" else "est"

    fn1 <- paste0("F", factor1); fn2 <- paste0("F", factor2)
    ld1 <- ld[ld$lhs == fn1, c("rhs", est_col)]
    ld2 <- ld[ld$lhs == fn2, c("rhs", est_col)]
    arrows <- merge(ld1, ld2, by = "rhs", suffixes = c("_1", "_2"))
    names(arrows) <- c("variable", "x", "y")
    arrows$label <- .bio_label(arrows$variable, fit$indicator_labels)

    # Scale arrows to fit in score space
    score_range <- max(abs(c(scores_df$f1, scores_df$f2))) * 0.8
    arrow_range <- max(abs(c(arrows$x, arrows$y)))
    scale_factor <- score_range / (arrow_range + 1e-10)
    arrows$x <- arrows$x * scale_factor
    arrows$y <- arrows$y * scale_factor

    g <- ggplot2::ggplot()

    if (!is.null(group_var)) {
        scores_df$group <- as.factor(scores_df[[group_var]])
        if (!is.null(group_labels)) {
            levels(scores_df$group) <- group_labels[levels(scores_df$group)]
        }
        g <- g + ggplot2::geom_point(data = scores_df,
                                      ggplot2::aes(x = .data$f1, y = .data$f2,
                                                   colour = .data$group),
                                      alpha = 0.35, size = 1.2) +
            ggplot2::scale_colour_manual(values = .pal(15), name = "Group")
    } else {
        g <- g + ggplot2::geom_point(data = scores_df,
                                      ggplot2::aes(x = .data$f1, y = .data$f2),
                                      alpha = 0.35, size = 1.2, colour = "grey60")
    }

    g + ggplot2::geom_segment(data = arrows,
                               ggplot2::aes(x = 0, y = 0,
                                            xend = .data$x, yend = .data$y),
                               arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                               colour = .pal(1), linewidth = 0.5, alpha = 0.8) +
        ggplot2::geom_text(data = arrows,
                           ggplot2::aes(x = .data$x, y = .data$y,
                                        label = .data$label),
                           size = 2.8, nudge_y = score_range * 0.04,
                           colour = "grey20") +
        ggplot2::labs(x = fit$factor_labels[factor1],
                      y = fit$factor_labels[factor2],
                      title = "Score Biplot",
                      subtitle = "Arrows: standardized loadings (scaled)") +
        theme_biomimic()
}


# ===========================================================================
# MULTI-GROUP CFA PLOTS
# ===========================================================================

#' @noRd
.plot_invariance <- function(fit, ...) {
    inv <- test_invariance(fit)
    fi <- inv$fit_indices

    # Reshape for plotting
    df <- data.frame(
        level = rep(fi$level, 3),
        metric = rep(c("CFI", "RMSEA", "SRMR"), each = nrow(fi)),
        value = c(fi$cfi, fi$rmsea, fi$srmr),
        stringsAsFactors = FALSE
    )
    df$level <- factor(df$level, levels = c("configural", "metric",
                                             "scalar", "strict"))

    thresholds <- data.frame(
        metric = c("CFI", "RMSEA", "SRMR"),
        threshold = c(0.90, 0.08, 0.08),
        excellent = c(0.95, 0.05, 0.05),
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$level, y = .data$value,
                                      group = 1)) +
        ggplot2::geom_line(colour = .pal(1), linewidth = 0.8) +
        ggplot2::geom_point(size = 3, colour = .pal(1)) +
        ggplot2::geom_hline(data = thresholds,
                            ggplot2::aes(yintercept = .data$threshold),
                            linetype = "dashed", colour = "grey50") +
        ggplot2::geom_hline(data = thresholds,
                            ggplot2::aes(yintercept = .data$excellent),
                            linetype = "dotted", colour = "#009E73") +
        ggplot2::facet_wrap(~ .data$metric, scales = "free_y") +
        ggplot2::labs(x = "Invariance Level", y = "Fit Index",
                      title = "Measurement Invariance Trajectory",
                      subtitle = "Dashed = acceptable; dotted = excellent") +
        theme_biomimic()
}


#' @noRd
.plot_loadings_group <- function(fit, factor = 1L, ...) {
    if (fit$model_type != "multigroup") {
        stop("loadings_group is only available for multigroup models.")
    }

    lam_list <- lavaan::lavInspect(fit$lavaan_fit, "est")
    if (!is.list(lam_list) || length(lam_list) < 2) {
        stop("Could not extract per-group parameters.")
    }

    groups <- names(lam_list)
    fname <- paste0("F", factor)

    df_list <- lapply(seq_along(groups), function(g) {
        lam <- lam_list[[g]]$lambda
        if (fname %in% colnames(lam)) {
            data.frame(
                variable = rownames(lam),
                loading  = lam[, fname],
                group    = groups[g],
                stringsAsFactors = FALSE
            )
        }
    })
    df <- do.call(rbind, df_list)
    df$label <- .bio_label(df$variable, fit$indicator_labels)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$loading, y = .data$label,
                                      colour = .data$group)) +
        ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.5)) +
        ggplot2::geom_segment(
            ggplot2::aes(x = 0, xend = .data$loading, yend = .data$label),
            position = ggplot2::position_dodge(0.5), linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
        ggplot2::scale_colour_manual(values = .pal(15), name = "Group") +
        ggplot2::labs(x = "Loading estimate", y = NULL,
                      title = "Per-Group Factor Loadings",
                      subtitle = "Large differences indicate differential item functioning (DIF)") +
        theme_biomimic()
}


#' @noRd
.plot_delta_cfi <- function(fit, ...) {
    inv <- test_invariance(fit)
    fi <- inv$fit_indices

    transitions <- data.frame(
        transition = c("Configural\n\u2192 Metric",
                        "Metric\n\u2192 Scalar",
                        "Scalar\n\u2192 Strict"),
        delta_cfi = diff(fi$cfi),
        stringsAsFactors = FALSE
    )
    transitions$transition <- factor(transitions$transition,
                                      levels = transitions$transition)
    transitions$pass <- abs(transitions$delta_cfi) < 0.01

    ggplot2::ggplot(transitions, ggplot2::aes(x = .data$transition,
                                               y = .data$delta_cfi,
                                               fill = .data$pass)) +
        ggplot2::geom_col(width = 0.6, alpha = 0.85) +
        ggplot2::geom_hline(yintercept = c(-0.01, 0.01), linetype = "dashed",
                            colour = "grey50") +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "#D55E00", "TRUE" = "#009E73"),
            labels = c("Fail (|\u0394CFI| \u2265 0.01)",
                        "Pass (|\u0394CFI| < 0.01)"),
            name = NULL) +
        ggplot2::labs(x = NULL, y = "\u0394CFI",
                      title = "Invariance Assessment: \u0394CFI Criterion",
                      subtitle = "Cheung & Rensvold (2002): |\u0394CFI| < 0.01") +
        theme_biomimic()
}


# ===========================================================================
# DIRECT EFFECTS PLOTS
# ===========================================================================

#' @noRd
.plot_decomposition <- function(fit, ...) {
    if (fit$model_type != "direct") {
        stop("decomposition is only available for MIMIC + direct effects models.")
    }

    th <- test_hypotheses(fit)
    med <- th$mediation
    if (is.null(med) || nrow(med) == 0) {
        stop("No mediation decomposition available.")
    }

    # Reshape for stacked bar
    df <- data.frame(
        protein = rep(med$protein, 2),
        type    = rep(c("Indirect (via pathway)", "Direct"), each = nrow(med)),
        effect  = c(med$indirect, med$direct),
        stringsAsFactors = FALSE
    )
    df$type <- factor(df$type, levels = c("Indirect (via pathway)", "Direct"))
    df$protein <- factor(df$protein, levels = rev(med$protein))

    ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$protein,
                                      fill = .data$type)) +
        ggplot2::geom_col(width = 0.6, alpha = 0.85,
                          position = ggplot2::position_stack()) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
        ggplot2::scale_fill_manual(values = c(.pal(2)[2], .pal(2)[1]),
                                    name = "Effect type") +
        ggplot2::labs(x = "Effect magnitude", y = NULL,
                      title = "Effect Decomposition: Indirect vs Direct",
                      subtitle = "How much of the group effect is mediated by the latent pathway?") +
        theme_biomimic()
}


# ===========================================================================
# STAGE 1-2 PLOTS
# ===========================================================================

#' @noRd
.plot_screening <- function(fit, top_n = 30L, ...) {
    .plot_screening_sel(fit$selection, top_n = top_n, ...)
}

#' @noRd
.plot_scree <- function(fit, ...) {
    .plot_scree_sel(fit$selection, ...)
}

#' @noRd
.plot_screening_sel <- function(selection, top_n = 30L, ...) {
    ranking <- selection$ranking
    n_show <- min(top_n, nrow(ranking))
    df <- ranking[seq_len(n_show), ]
    df$selected <- df$variable %in% selection$variables
    df$has_direct <- df$variable %in% selection$direct_effects
    df$variable <- factor(df$variable, levels = rev(df$variable))

    g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$snr, y = .data$variable,
                                           fill = .data$selected)) +
        ggplot2::geom_col(width = 0.7, alpha = 0.85) +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "grey75", "TRUE" = .pal(2)[2]),
            labels = c("Not selected", "Selected"), name = NULL)

    # Mark direct effects
    if (any(df$has_direct)) {
        de_df <- df[df$has_direct, ]
        g <- g + ggplot2::geom_point(data = de_df,
                                      ggplot2::aes(x = .data$snr + max(df$snr) * 0.02,
                                                   y = .data$variable),
                                      shape = 18, size = 3, colour = "#CC3311",
                                      inherit.aes = FALSE)
    }

    g + ggplot2::labs(x = "Signal-to-Noise Ratio (SNR)", y = NULL,
                      title = "VB-ARD Variable Screening (Stage 1)",
                      subtitle = sprintf("Top %d of %d variables. Red diamonds = direct effects.",
                                         n_show, nrow(ranking))) +
        theme_biomimic()
}


#' @noRd
.plot_scree_sel <- function(selection, ...) {
    # SVD of residual matrix from VB-ARD
    fit_vb <- selection$fit
    Y_resid <- fit_vb$Y - fit_vb$X %*% solve(crossprod(fit_vb$X),
                                                crossprod(fit_vb$X, fit_vb$Y))
    sv <- svd(Y_resid, nu = 0, nv = 0)
    d <- sv$d
    n_show <- min(15, length(d))
    var_expl <- (d^2) / sum(d^2)

    df <- data.frame(
        component = seq_len(n_show),
        variance  = var_expl[seq_len(n_show)],
        cumulative = cumsum(var_expl[seq_len(n_show)]),
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(df, ggplot2::aes(x = .data$component, y = .data$variance)) +
        ggplot2::geom_line(colour = .pal(1), linewidth = 0.8) +
        ggplot2::geom_point(size = 3, colour = .pal(1)) +
        ggplot2::geom_line(ggplot2::aes(y = .data$cumulative),
                           linetype = "dashed", colour = "grey50") +
        ggplot2::geom_point(ggplot2::aes(y = .data$cumulative),
                            size = 2, colour = "grey50", shape = 1) +
        ggplot2::scale_x_continuous(breaks = seq_len(n_show)) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(1)) +
        ggplot2::labs(x = "Component", y = "Proportion of variance",
                      title = "Scree Plot (Stage 2)",
                      subtitle = "Solid: individual; dashed: cumulative") +
        theme_biomimic()
}


# ===========================================================================
# Legacy wrappers (backwards compatibility)
# ===========================================================================

#' @rdname plot.biomimic_fit
#' @param fit A \code{biomimic_fit} object (used by legacy wrappers).
#' @export
plot_loadings <- function(fit, ...) plot.biomimic_fit(fit, type = "loadings", ...)

#' @rdname plot.biomimic_fit
#' @export
plot_scores <- function(fit, ...) plot.biomimic_fit(fit, type = "scores", ...)

#' @rdname plot.biomimic_fit
#' @export
plot_residuals <- function(fit, ...) plot.biomimic_fit(fit, type = "residuals", ...)

#' @rdname plot.biomimic_fit
#' @export
plot_mahalanobis <- function(fit, ...) plot.biomimic_fit(fit, type = "mahalanobis", ...)

#' @rdname plot.biomimic_fit
#' @export
plot_screening <- function(fit, ...) {
    if (inherits(fit, "biomimic_selection")) {
        plot.biomimic_selection(fit, type = "screening", ...)
    } else {
        plot.biomimic_fit(fit, type = "screening", ...)
    }
}

#' @rdname plot.biomimic_fit
#' @export
plot_path <- function(fit, ...) {
    stopifnot(inherits(fit, "biomimic_fit"))
    cat("=== Path Diagram ===\n\n")
    pe <- lavaan::parameterEstimates(fit$lavaan_fit)

    structural <- pe[pe$op == "~" & grepl("^F[0-9]+$", pe$lhs), ]
    if (nrow(structural) > 0) {
        cat("Structural (Covariates -> Latent Pathways):\n")
        for (i in seq_len(nrow(structural))) {
            sig <- .signif_stars(structural$pvalue[i])
            cat(sprintf("  %s --[%.3f%s]--> %s\n",
                        structural$rhs[i], structural$est[i], sig,
                        .bio_label(structural$lhs[i], .factor_label_map(fit))))
        }
        cat("\n")
    }

    loadings <- pe[pe$op == "=~", ]
    if (nrow(loadings) > 0) {
        cat("Measurement (Latent Pathways -> Proteins):\n")
        for (i in seq_len(nrow(loadings))) {
            sig <- .signif_stars(loadings$pvalue[i])
            cat(sprintf("  %s --[%.3f%s]--> %s\n",
                        .bio_label(loadings$lhs[i], .factor_label_map(fit)),
                        loadings$est[i], sig,
                        .bio_label(loadings$rhs[i], fit$indicator_labels)))
        }
        cat("\n")
    }

    direct <- pe[pe$op == "~" & !grepl("^F[0-9]+$", pe$lhs), ]
    if (nrow(direct) > 0) {
        cat("Direct Effects (Covariates -> Proteins):\n")
        for (i in seq_len(nrow(direct))) {
            sig <- .signif_stars(direct$pvalue[i])
            cat(sprintf("  %s --[%.3f%s]--> %s\n",
                        direct$rhs[i], direct$est[i], sig,
                        .bio_label(direct$lhs[i], fit$indicator_labels)))
        }
    }
    cat("\nFor graphical diagrams: semPlot::semPaths(fit$lavaan_fit)\n")
    invisible(fit)
}
