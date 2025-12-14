#' @title Trace plot of probabilities parameter
#' @description It returns a trace plot of probabilities parameter after burn-in. Different colors represent different groups.
#'
#' @param Gibbs Gibbs sample from \code{BCFM} function.
#' @param burnin Number of burn-in period. If not specified, it uses the first tenths sample as burn-in.
#' @param main.title Title of the plot. Default is "Trace Plot: Cluster Probabilities"
#' @param x.label X-axis label. Default is "BCFM Iteration (post burn-in)"
#' @param y.label Y-axis label. Default is "Probability"
#' @param cluster_names Character vector of cluster names. If NULL, defaults to Cluster 1, Cluster 2, etc.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual theme_minimal labs theme ylim
#' @importFrom grDevices colorRampPalette
#'
#' @return A ggplot object showing trace plots
#' @export ggplot_probs.trace
#' @export

ggplot_probs.trace <- function(Gibbs,
                              burnin = NA,
                              main.title = "Trace Plot: Cluster Probabilities",
                              x.label = "BCFM Iteration (post burn-in)",
                              y.label = "Probability",
                              cluster_names = NULL){

  if(is.na(burnin)){burnin <- round(dim(Gibbs$probs)[1] / 10)}
  n.iter <- dim(Gibbs$probs)[1]
  G <- dim(Gibbs$probs)[2]
  turns <- seq(burnin + 1, n.iter)

  # Set cluster names
  if(is.null(cluster_names)){
    cluster_names <- paste0("Cluster ", 1:G)
  } else if(length(cluster_names) != G){
    warning("Length of cluster_names does not match number of clusters (G=", G,
            "). Using default labels.")
    cluster_names <- paste0("Cluster ", 1:G)
  }

  # Generate colors
  color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
  colors <- color_palette(G)
  names(colors) <- cluster_names

  # Prepare data
  probs_trace <- as.data.frame(Gibbs$probs[turns, ])
  colnames(probs_trace) <- cluster_names
  probs_trace$Iteration <- 1:nrow(probs_trace)

  # Convert to long format
  probs_trace <- pivot_longer(probs_trace,
                               -c(Iteration),
                               values_to = "Probability",
                               names_to = "Cluster")

  # Set factor levels
  probs_trace$Cluster <- factor(probs_trace$Cluster, levels = cluster_names)

  # Create plot
  p <- ggplot(probs_trace, aes(x = Iteration, y = Probability, color = Cluster)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_color_manual(values = colors) +
    ylim(0, 1) +
    theme_minimal() +
    labs(title = main.title,
         x = x.label,
         y = y.label,
         color = "Cluster") +
    theme(legend.position = "right")

  return(p)
}
