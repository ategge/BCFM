#' @title The density plot of the diagonal of group covariance, Omega, with ggplot2
#' @description It returns multiple plots of the diagonal of group covariance, Omega using ggplot2 and gridExtra. It returns the result by each factor, different colors representing different factors
#'
#' @param Gibbs Gibbs sample from BCFM
#' @param group.select Group/cluster to plot. If not specified, the first group will be used.
#' @param true.val True values of Omega, if applicable.
#' @param burnin Number of burn-in period. If not specified, the first tenths is used as burn-in.
#' @param main.title Main title for the plot. Default is "Posterior Densities of Omega"
#' @param x.label X-axis label. Default is "Value"
#' @param y.label Y-axis label. Default is "Density"
#' @param factor_labels Character vector of factor names. If NULL, defaults to Factor 1, Factor 2, etc.
#' @param show.offdiag Show off-diagonal elements. Default is TRUE for any k.
#'
#' @importFrom ggplot2 aes geom_density labs theme_minimal scale_color_manual theme element_blank geom_vline
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices colorRampPalette
#'
#' @return A ggplot object showing densities of Omega elements
#' @export ggplot.omega.density
#' @export

ggplot.omega.density <- function(Gibbs,
                                  group.select = 1,
                                  true.val = NA,
                                  burnin = NA,
                                  main.title = "Posterior Densities of Omega",
                                  x.label = "Value",
                                  y.label = "Density",
                                  factor_labels = NULL,
                                  show.offdiag = TRUE){

  n.iter <- dim(Gibbs$Omega)[1]
  G <- dim(Gibbs$Omega)[2]
  k <- dim(Gibbs$Omega)[3]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  turns <- seq(burnin + 1, n.iter)

  # Set factor labels
  if(is.null(factor_labels)){
    factor_labels <- paste0("Factor ", 1:k)
  } else if(length(factor_labels) != k){
    warning("Length of factor_labels does not match number of factors (k=", k,
            "). Using default labels.")
    factor_labels <- paste0("Factor ", 1:k)
  }

  # Generate colors
  color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
  colors <- color_palette(k)
  names(colors) <- factor_labels

  Density <- Group <- NULL

  # Diagonal elements
  Omega.diag.density <- rep(NA, length(turns) * k)
  for(j in 1:k){
    Omega.diag.density[(turns + (j-1)*length(turns) - burnin)] <- Gibbs$Omega[turns, group.select, j, j]
  }
  Omega.diag.density <- as.data.frame(cbind(rep(factor_labels, each = length(turns)), Omega.diag.density))
  colnames(Omega.diag.density) <- c("Factor", "Density")
  Omega.diag.density$Density <- as.numeric(Omega.diag.density$Density)
  Omega.diag.density$Factor <- factor(Omega.diag.density$Factor, levels = factor_labels)

  p1 <- ggplot(Omega.diag.density, aes(x = Density, color = Factor)) +
    geom_density(linewidth = 1) +
    labs(title = "Diagonal Elements (Variances)",
         x = x.label,
         y = y.label) +
    theme_minimal() +
    scale_color_manual(values = colors) +
    theme(legend.position = "bottom")

  if(!is.na(true.val[1])){
    true_diag_data <- data.frame(
      xintercept = diag(true.val[group.select, , ]),
      Factor = factor_labels
    )
    p1 <- p1 + geom_vline(data = true_diag_data,
                          aes(xintercept = xintercept, color = Factor),
                          linetype = "dashed",
                          linewidth = 0.8)
  }

  # Off-diagonal elements (if k > 1 and show.offdiag = TRUE)
  if(k > 1 & show.offdiag){
    # Get all unique off-diagonal pairs
    offdiag.pairs <- which(upper.tri(matrix(1, k, k)), arr.ind = TRUE)
    n.offdiag <- nrow(offdiag.pairs)

    Omega.offdiag.density <- rep(NA, length(turns) * n.offdiag)
    pair_labels <- character(n.offdiag)

    for(j in 1:n.offdiag){
      row_idx <- offdiag.pairs[j, 1]
      col_idx <- offdiag.pairs[j, 2]
      Omega.offdiag.density[(turns + (j-1)*length(turns) - burnin)] <-
        Gibbs$Omega[turns, group.select, row_idx, col_idx]
      pair_labels[j] <- paste0(factor_labels[row_idx], "-", factor_labels[col_idx])
    }

    Omega.offdiag.density <- as.data.frame(cbind(rep(pair_labels, each = length(turns)),
                                                   Omega.offdiag.density))
    colnames(Omega.offdiag.density) <- c("Pair", "Density")
    Omega.offdiag.density$Density <- as.numeric(Omega.offdiag.density$Density)
    Omega.offdiag.density$Pair <- factor(Omega.offdiag.density$Pair, levels = pair_labels)

    # Generate colors for pairs
    color_palette_offdiag <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
    colors_offdiag <- color_palette_offdiag(n.offdiag)
    names(colors_offdiag) <- pair_labels

    p2 <- ggplot(Omega.offdiag.density, aes(x = Density, color = Pair)) +
      geom_density(linewidth = 1) +
      labs(title = "Off-Diagonal Elements (Covariances)",
           x = x.label,
           y = y.label) +
      theme_minimal() +
      scale_color_manual(values = colors_offdiag) +
      theme(legend.position = "bottom")

    if(!is.na(true.val[1])){
      true_offdiag_data <- data.frame(
        xintercept = numeric(n.offdiag),
        Pair = pair_labels
      )
      for(j in 1:n.offdiag){
        row_idx <- offdiag.pairs[j, 1]
        col_idx <- offdiag.pairs[j, 2]
        true_offdiag_data$xintercept[j] <- true.val[group.select, row_idx, col_idx]
      }
      p2 <- p2 + geom_vline(data = true_offdiag_data,
                            aes(xintercept = xintercept, color = Pair),
                            linetype = "dashed",
                            linewidth = 0.8)
    }

    # Arrange both plots
    p <- grid.arrange(p1, p2, ncol = 2, top = main.title)
  } else {
    # Only diagonal plot
    p <- grid.arrange(p1, ncol = 1, top = main.title)
  }

  return(p)
}
