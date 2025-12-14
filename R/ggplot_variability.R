#' @title Variability explained by factors
#' @description Plots the proportion of variability explained by each factor based on eigenvalues of the correlation matrix. Useful for determining the number of factors.
#'
#' @param data The data matrix
#' @param nfactors Number of factors to display. Default is 5.
#' @param main.title Main title for the plot. Default is "Proportion of Variability Explained by Factors"
#' @param x.label X-axis label. Default is "Number of Factors"
#' @param y.label Y-axis label. Default is "Proportion of Variability"
#'
#' @importFrom stats cor
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal labs scale_x_continuous
#'
#' @return A ggplot object
#' @export ggplot_variability
#' @export

ggplot_variability <- function(data,
                              nfactors = 5,
                              main.title = "Proportion of Variability Explained by Factors",
                              x.label = "Number of Factors",
                              y.label = "Proportion of Variability"){

  # Calculate eigenvalues
  eigvals <- eigen(cor(data))$values

  # Ensure nfactors doesn't exceed available factors
  nfactors <- min(nfactors, length(eigvals))

  # Prepare data
  variability_data <- data.frame(
    Factor = 1:nfactors,
    Proportion = eigvals[1:nfactors] / ncol(data)
  )

  # Create plot
  p <- ggplot(variability_data, aes(x = Factor, y = Proportion)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    scale_x_continuous(breaks = 1:nfactors) +
    theme_minimal() +
    labs(title = main.title,
         x = x.label,
         y = y.label)

  return(p)
}
