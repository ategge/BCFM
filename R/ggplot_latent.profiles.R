#' @title Plot Latent Factor Profiles by Cluster
#' @description Visualizes the latent factor profile means (mu) for each cluster, similar to Latent Profile Analysis (LPA) plots
#'
#' @param Gibbs Gibbs sample derived from \code{BCFM} function
#' @param burnin Number of burn-in period. If not specified, it uses the first tenth as burn-in period
#' @param factor_labels Character vector of factor names. If NULL, defaults to Factor 1, Factor 2, etc.
#' @param cluster_names Character vector of cluster names. If NULL, defaults to Cluster 1, Cluster 2, etc.
#' @param colors Named vector of colors for each cluster. If NULL, uses default color palette
#' @param title Plot title. Default is "Latent Factor Profiles by Cluster"
#' @param x_label X-axis label. Default is "Factor"
#' @param y_label Y-axis label. Default is "Posterior Mean"
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_discrete scale_color_manual theme_minimal labs theme element_text
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' ggplot_latent.profiles(Gibbs = SDresult$Result)
#'
#' # With custom labels
#' ggplot_latent.profiles(
#'   Gibbs = SDresult$Result,
#'   factor_labels = c("Cigarette", "E-cig", "Pipe", "Hookah", "Snus", "Cigar"),
#'   cluster_names = c("Traditionalists", "Experimenters", "Hybrid Users",
#'                     "Social Smokers", "Nicotine Enthusiasts")
#' )
#' }
#' @export ggplot_latent.profiles
#' @export
ggplot_latent.profiles <- function(Gibbs,
                                 burnin = NA,
                                 factor_labels = NULL,
                                 cluster_names = NULL,
                                 colors = NULL,
                                 title = "Latent Factor Profiles by Cluster",
                                 x_label = "Factor",
                                 y_label = "Posterior Mean") {

  # Packages loaded via package dependencies

  # --- Extract dimensions ---
  mu_draws <- Gibbs$mu
  n.iter <- dim(mu_draws)[1]
  G <- dim(mu_draws)[2]  # number of clusters
  K <- dim(mu_draws)[3]  # number of factors

  # --- Burn-in removal ---
  if(is.na(burnin)) {
    burnin <- round(n.iter / 10)
  }
  keep_iters <- seq(burnin + 1, n.iter)

  # --- Compute posterior means of mu for each cluster and factor ---
  mu_summary <- apply(mu_draws[keep_iters, , ], MARGIN = c(2, 3), FUN = mean)
  # mu_summary is [G x K] matrix: rows = clusters, cols = factors

  # --- Set default labels if not provided ---
  if(is.null(factor_labels)) {
    factor_labels <- paste0("Factor ", 1:K)
  } else if(length(factor_labels) != K) {
    warning("Length of factor_labels does not match number of factors (K=", K,
            "). Using default labels.")
    factor_labels <- paste0("Factor ", 1:K)
  }

  if(is.null(cluster_names)) {
    cluster_names <- paste0("Cluster ", 1:G)
  } else if(length(cluster_names) != G) {
    warning("Length of cluster_names does not match number of clusters (G=", G,
            "). Using default labels.")
    cluster_names <- paste0("Cluster ", 1:G)
  }

  # --- Put into dataframe ---
  cluster_factor_summary <- as.data.frame(mu_summary)
  colnames(cluster_factor_summary) <- factor_labels
  cluster_factor_summary$cluster_name <- cluster_names

  # --- Convert to long format ---
  cluster_factor_long <- cluster_factor_summary %>%
    pivot_longer(
      cols = all_of(factor_labels),
      names_to = "Factor",
      values_to = "Posterior_Mean"
    )

  # --- Set factor order ---
  cluster_factor_long$Factor <- factor(cluster_factor_long$Factor, levels = factor_labels)
  cluster_factor_long$cluster_name <- factor(cluster_factor_long$cluster_name, levels = cluster_names)

  # --- Define colors ---
  if(is.null(colors)) {
    # Use colorRampPalette to generate colors for any number of clusters
    color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
    default_colors <- color_palette(G)
    colors <- setNames(default_colors, cluster_names)
  } else if(length(colors) != G || !all(cluster_names %in% names(colors))) {
    warning("Colors vector does not match cluster names. Using default colors.")
    color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
    default_colors <- color_palette(G)
    colors <- setNames(default_colors, cluster_names)
  }

  # --- Create plot ---
  p <- ggplot(cluster_factor_long,
              aes(x = Factor, y = Posterior_Mean, group = cluster_name, color = cluster_name)) +
    geom_line(linewidth = 1, lineend = "round") +
    geom_point(size = 7, colour = "white", stroke = 0) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(title = title,
         x = x_label,
         y = y_label,
         color = "Cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  return(p)
}
