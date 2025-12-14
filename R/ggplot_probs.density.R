#' @title Density plot for posterior of probabilities
#' @description The function returns a density plot of the cluster assignment probability, p. Different color represent different Clusters.
#'
#' @param Gibbs MCMC sample simulated from \code{BCFMR} or \code{BCFMcpp} function
#' @param burnin Number of burn-in period. If not specified, it uses the first tenths as burn-in.
#' @param truep True values of probabilities. If not available, NA.
#' @param main.title Title of the plot. Default is "Posterior Densities of Cluster Probabilities"
#' @param x.label X-axis label. Default is "Probability"
#' @param y.label Y-axis label. Default is "Density"
#' @param cluster_names Character vector of cluster names. If NULL, defaults to Cluster 1, Cluster 2, etc.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes geom_density geom_vline theme_minimal scale_color_manual scale_fill_manual xlim labs theme
#' @importFrom grDevices colorRampPalette
#'
#' @return A ggplot object (grob from grid.arrange) with plot and legend
#' @export ggplot_probs.density
#' @export

ggplot_probs.density <- function(Gibbs,
                                  burnin = NA,
                                  truep = NA,
                                  main.title = "Posterior Densities of Cluster Probabilities",
                                  x.label = "Probability",
                                  y.label = "Density",
                                  cluster_names = NULL){

  n.all <- dim(Gibbs$probs)[1]
  G <- dim(Gibbs$probs)[2]

  if(is.na(burnin)){burnin <- round(n.all / 10)}
  turns <- seq(burnin + 1, n.all)

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

  # Prepare data (now using burnin!)
  Density <- data.frame()
  for(i in 1:G){
    temp_df <- data.frame(
      Cluster = cluster_names[i],
      Probability = Gibbs$probs[turns, i]
    )
    Density <- rbind(Density, temp_df)
  }

  Density$Cluster <- factor(Density$Cluster, levels = cluster_names)

  # Create main plot
  probsplot <- ggplot(Density, aes(x = Probability, color = Cluster)) +
    geom_density(linewidth = 1) +
    theme_minimal() +
    scale_color_manual(values = colors) +
    xlim(0, 1) +
    labs(title = main.title,
         x = x.label,
         y = y.label) +
    theme(legend.position = "none")

  # Add true values if available
  if(!is.na(truep[1])){
    true_data <- data.frame(
      xintercept = truep,
      Cluster = cluster_names
    )
    probsplot <- probsplot +
      geom_vline(data = true_data,
                 aes(xintercept = xintercept, color = Cluster),
                 linetype = "dashed",
                 linewidth = 0.5)
  }

  # Create legend plot
  probsplot.legend <- ggplot(Density, aes(x = Probability, fill = Cluster)) +
    geom_density() +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(legend.position = "right")

  # Extract legend
  legend <- get_legend(probsplot.legend)

  # Arrange plot with legend on the right
  p <- grid.arrange(probsplot, legend,
                    layout_matrix = rbind(c(1, 1, 1, 2)))

  return(p)
}

# Helper function to extract legend
get_legend <- function(plot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
