# --------------------------------------------------------------------------
# Network Analysis: SNR-annotated correlation and partial correlation graphs
# Connects VB-ARD screening with biological network interpretation
# --------------------------------------------------------------------------

#' Network Analysis of VB-ARD Screening Results
#'
#' Builds a correlation or partial correlation network from the observed
#' data, annotated with VB-ARD signal-to-noise ratios. Identifies the
#' neighbourhood of top-SNR proteins for pathway extension analysis.
#'
#' @section Network Types:
#' \describe{
#'   \item{\code{"correlation"}}{Pearson (or Spearman) correlations. Edges
#'     represent marginal associations. Fast and intuitive but may include
#'     indirect associations.}
#'   \item{\code{"partial"}}{Partial correlations via shrinkage estimation
#'     (Ledoit-Wolf or Schaefer-Strimmer). Edges represent direct
#'     associations, controlling for all other variables. Requires the
#'     \code{corpcor} package.}
#' }
#'
#' @param fit A \code{vb_ard} object (from Stage 1) or a
#'   \code{biomimic_selection} object (from Stage 2).
#' @param type Character, \code{"correlation"} (default) or \code{"partial"}.
#' @param threshold Numeric, minimum absolute (partial) correlation for an
#'   edge (default 0.3 for correlation, 0.15 for partial).
#' @param top_n Integer, number of top-SNR variables to highlight as "hub"
#'   nodes (default 15, matching typical top-k selection).
#' @param cor_method Character, correlation method for type="correlation":
#'   \code{"pearson"} (default) or \code{"spearman"}.
#' @param max_vars Integer, maximum number of variables to include in the
#'   network (default 100). Variables are selected by SNR ranking. Set to
#'   \code{NULL} to include all.
#'
#' @return A list of class \code{biomimic_network} with components:
#' \describe{
#'   \item{graph}{An \code{igraph} object with node attributes \code{snr},
#'     \code{selected} (top_n), and \code{label}; edge attribute \code{weight}
#'     (absolute correlation) and \code{sign}.}
#'   \item{adjacency}{The thresholded (partial) correlation matrix.}
#'   \item{full_cor}{The full (partial) correlation matrix before thresholding.}
#'   \item{hub_nodes}{Character vector of top-SNR variable names.}
#'   \item{marginal_neighbours}{Data frame of non-selected variables connected
#'     to hub nodes, with columns: \code{variable}, \code{snr},
#'     \code{hub_connected_to}, \code{edge_weight}.}
#'   \item{modularity}{Newman-Girvan modularity of detected communities.}
#'   \item{communities}{Community membership vector.}
#'   \item{hub_density}{Edge density within hub subgraph vs overall.}
#'   \item{type}{Network type used.}
#'   \item{threshold}{Threshold used.}
#' }
#'
#' @examples
#' \dontrun{
#' vb <- vb_ard(Y, X, K = 1)
#' net <- network_analysis(vb, type = "correlation", threshold = 0.3)
#' print(net)
#' plot(net)
#' }
#'
#' @export
network_analysis <- function(fit,
                              type = c("correlation", "partial"),
                              threshold = NULL,
                              top_n = 15L,
                              cor_method = "pearson",
                              max_vars = 100L) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for network analysis.\n",
             "  install.packages('igraph')")
    }

    type <- match.arg(type)

    # Default thresholds
    if (is.null(threshold)) {
        threshold <- if (type == "correlation") 0.3 else 0.15
    }

    # Extract data and SNR from fit object
    if (inherits(fit, "biomimic_selection")) {
        vb <- fit$fit
    } else if (inherits(fit, "vb_ard")) {
        vb <- fit
    } else {
        stop("fit must be a vb_ard or biomimic_selection object.")
    }

    Y <- vb$Y
    var_names <- vb$var_names
    snr_mat <- vb$snr_lambda
    # Max SNR across factors for each variable
    snr_vec <- apply(snr_mat, 1, max)
    names(snr_vec) <- var_names

    # Subset to top max_vars by SNR
    if (!is.null(max_vars) && length(var_names) > max_vars) {
        keep_idx <- order(snr_vec, decreasing = TRUE)[seq_len(max_vars)]
        Y <- Y[, keep_idx, drop = FALSE]
        var_names <- var_names[keep_idx]
        snr_vec <- snr_vec[keep_idx]
    }

    p <- length(var_names)

    .log_info("network_analysis: type=%s, threshold=%.3f, p=%d variables",
              type, threshold, p)

    # --- Compute correlation matrix ---
    if (type == "correlation") {
        cor_mat <- stats::cor(Y, method = cor_method)
    } else {
        if (!requireNamespace("corpcor", quietly = TRUE)) {
            stop("Package 'corpcor' is required for partial correlation.\n",
                 "  install.packages('corpcor')")
        }
        # Shrinkage-based partial correlation (Schaefer-Strimmer)
        cor_mat <- corpcor::cor2pcor(corpcor::cor.shrink(Y, verbose = FALSE))
        dimnames(cor_mat) <- list(var_names, var_names)
    }

    full_cor <- cor_mat

    # --- Threshold and build adjacency ---
    adj <- cor_mat
    adj[abs(adj) < threshold] <- 0
    diag(adj) <- 0

    # --- Build igraph ---
    g <- igraph::graph_from_adjacency_matrix(
        abs(adj), mode = "undirected", weighted = TRUE, diag = FALSE
    )

    # Edge sign
    edges <- igraph::as_data_frame(g, what = "edges")
    if (nrow(edges) > 0) {
        edge_signs <- vapply(seq_len(nrow(edges)), function(i) {
            r <- cor_mat[edges$from[i], edges$to[i]]
            ifelse(r >= 0, "positive", "negative")
        }, character(1))
        igraph::E(g)$sign <- edge_signs
    }

    # Node attributes
    igraph::V(g)$snr <- snr_vec[igraph::V(g)$name]
    hub_names <- names(sort(snr_vec, decreasing = TRUE))[seq_len(min(top_n, p))]
    igraph::V(g)$selected <- igraph::V(g)$name %in% hub_names
    igraph::V(g)$label <- igraph::V(g)$name

    # --- Community detection ---
    communities <- NULL
    modularity_val <- NA_real_
    if (igraph::ecount(g) > 0) {
        comm <- igraph::cluster_louvain(g)
        communities <- igraph::membership(comm)
        modularity_val <- igraph::modularity(comm)
        igraph::V(g)$community <- as.integer(communities[igraph::V(g)$name])
    }

    # --- Hub density ---
    hub_density <- NA_real_
    overall_density <- igraph::graph.density(g)
    hub_vertices <- which(igraph::V(g)$name %in% hub_names)
    if (length(hub_vertices) >= 2) {
        hub_subgraph <- igraph::induced_subgraph(g, hub_vertices)
        hub_density <- igraph::graph.density(hub_subgraph)
    }

    # --- Marginal neighbours ---
    non_hub <- setdiff(var_names, hub_names)
    neighbours_list <- list()

    for (hub in hub_names) {
        if (!hub %in% igraph::V(g)$name) next
        nb <- igraph::neighbors(g, hub)
        nb_names <- igraph::V(g)$name[nb]
        # Only non-hub neighbours
        marginal <- intersect(nb_names, non_hub)
        for (m in marginal) {
            w <- abs(cor_mat[hub, m])
            neighbours_list[[length(neighbours_list) + 1]] <- data.frame(
                variable = m,
                snr = snr_vec[m],
                hub_connected_to = hub,
                edge_weight = w,
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(neighbours_list) > 0) {
        marginal_neighbours <- do.call(rbind, neighbours_list)
        # Deduplicate: keep strongest connection per marginal variable
        marginal_neighbours <- marginal_neighbours[
            order(-marginal_neighbours$edge_weight), ]
        marginal_neighbours <- marginal_neighbours[
            !duplicated(marginal_neighbours$variable), ]
        rownames(marginal_neighbours) <- NULL
    } else {
        marginal_neighbours <- data.frame(
            variable = character(0), snr = numeric(0),
            hub_connected_to = character(0), edge_weight = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    .log_info("network_analysis: %d nodes, %d edges, modularity=%.3f, %d marginal neighbours",
              igraph::vcount(g), igraph::ecount(g), modularity_val,
              nrow(marginal_neighbours))

    result <- structure(
        list(
            graph                = g,
            adjacency            = adj,
            full_cor             = full_cor,
            hub_nodes            = hub_names,
            marginal_neighbours  = marginal_neighbours,
            modularity           = modularity_val,
            communities          = communities,
            hub_density          = hub_density,
            overall_density      = overall_density,
            snr                  = snr_vec,
            type                 = type,
            threshold            = threshold,
            top_n                = top_n
        ),
        class = "biomimic_network"
    )

    result
}


#' @export
print.biomimic_network <- function(x, ...) {
    cat("=== biomimic Network Analysis ===\n")
    cat(sprintf("  Type: %s (threshold = %.3f)\n", x$type, x$threshold))
    cat(sprintf("  Nodes: %d | Edges: %d\n",
                igraph::vcount(x$graph), igraph::ecount(x$graph)))
    cat(sprintf("  Hub nodes (top %d by SNR): %s\n",
                x$top_n, paste(head(x$hub_nodes, 5), collapse = ", ")))
    if (length(x$hub_nodes) > 5) cat(sprintf("    ... and %d more\n",
                                              length(x$hub_nodes) - 5))
    cat(sprintf("  Hub density: %.3f | Overall density: %.3f | Ratio: %.1fx\n",
                x$hub_density, x$overall_density,
                x$hub_density / max(x$overall_density, 1e-10)))
    cat(sprintf("  Modularity (Louvain): %.3f\n", x$modularity))
    cat(sprintf("  Marginal neighbours of hubs: %d variables\n",
                nrow(x$marginal_neighbours)))
    if (nrow(x$marginal_neighbours) > 0) {
        top_mn <- head(x$marginal_neighbours, 5)
        for (i in seq_len(nrow(top_mn))) {
            cat(sprintf("    %s (SNR=%.1f) <-> %s (|r|=%.3f)\n",
                        top_mn$variable[i], top_mn$snr[i],
                        top_mn$hub_connected_to[i], top_mn$edge_weight[i]))
        }
        if (nrow(x$marginal_neighbours) > 5)
            cat(sprintf("    ... and %d more\n",
                        nrow(x$marginal_neighbours) - 5))
    }
    invisible(x)
}


#' Plot Method for biomimic_network
#'
#' Produces a ggraph network plot with nodes sized and coloured by SNR.
#' Hub nodes (top-k) are highlighted. Edge colour indicates positive or
#' negative (partial) correlation.
#'
#' @param x A \code{biomimic_network} object.
#' @param layout Character, graph layout algorithm (default \code{"fr"}
#'   for Fruchterman-Reingold).
#' @param show_labels Logical, show node labels (default TRUE).
#' @param label_hubs_only Logical, only label hub nodes (default TRUE
#'   when > 30 nodes).
#' @param ... Ignored.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot.biomimic_network <- function(x, layout = "fr",
                                   show_labels = TRUE,
                                   label_hubs_only = NULL, ...) {
    .check_ggplot2()
    if (!requireNamespace("ggraph", quietly = TRUE)) {
        stop("Package 'ggraph' is required for network plots.\n",
             "  install.packages('ggraph')")
    }

    g <- x$graph
    if (igraph::vcount(g) == 0) {
        stop("Network has no nodes.")
    }

    # Auto: label only hubs if many nodes
    if (is.null(label_hubs_only)) {
        label_hubs_only <- igraph::vcount(g) > 30
    }

    # Prepare label column
    if (show_labels) {
        if (label_hubs_only) {
            igraph::V(g)$plot_label <- ifelse(
                igraph::V(g)$selected, igraph::V(g)$name, "")
        } else {
            igraph::V(g)$plot_label <- igraph::V(g)$name
        }
    } else {
        igraph::V(g)$plot_label <- ""
    }

    # Remove isolated nodes for cleaner plot
    isolated <- which(igraph::degree(g) == 0)
    if (length(isolated) > 0 && length(isolated) < igraph::vcount(g)) {
        g <- igraph::delete_vertices(g, isolated)
    }

    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        stop("No edges remain after filtering. Try a lower threshold.")
    }

    # SNR scaling for node size
    snr_vals <- igraph::V(g)$snr
    snr_scaled <- 2 + 6 * (snr_vals - min(snr_vals)) /
        max(max(snr_vals) - min(snr_vals), 1e-10)

    p <- ggraph::ggraph(g, layout = layout) +
        ggraph::geom_edge_link(
            ggplot2::aes(width = .data$weight, colour = .data$sign),
            alpha = 0.4) +
        ggraph::scale_edge_colour_manual(
            values = c("positive" = "#56B4E9", "negative" = "#D55E00"),
            name = paste(x$type, "sign")) +
        ggraph::scale_edge_width(range = c(0.2, 1.5), guide = "none") +
        ggplot2::geom_point(
            ggplot2::aes(x = .data$x, y = .data$y,
                         size = snr_scaled,
                         fill = .data$selected),
            shape = 21, colour = "grey30", stroke = 0.3) +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "grey80", "TRUE" = "#E69F00"),
            labels = c("Marginal", paste("Top", x$top_n, "(hub)")),
            name = "VB-ARD") +
        ggplot2::scale_size_identity()

    if (show_labels) {
        p <- p + ggplot2::geom_text(
            ggplot2::aes(x = .data$x, y = .data$y,
                         label = .data$plot_label),
            size = 2.5, repel = FALSE,
            nudge_y = 0.15, colour = "grey20")
    }

    p + ggplot2::labs(
            title = sprintf("Protein %s Network",
                            if (x$type == "partial") "Partial Correlation"
                            else "Correlation"),
            subtitle = sprintf(
                "|r| > %.2f | %d nodes, %d edges | modularity = %.3f | hub density %.1fx overall",
                x$threshold, igraph::vcount(g), igraph::ecount(g),
                x$modularity,
                x$hub_density / max(x$overall_density, 1e-10))
        ) +
        ggraph::theme_graph(base_family = "") +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            legend.position = "bottom"
        )
}
