#' @title Density of group means mu using ggplot2
#' @description The function returns multiple density plots of group mean parameter mu.
#'
#' @param Gibbs The Gibbs sample from BCFM function
#' @param true.val The true value of mu, if applicable
#' @param add.legend Add legend on extra pane
#' @param burnin Number of burn-in period. If not specified, it uses the first tenths as burn-in.
#' @param layout.dim Dimension of panes. If not specified, the plots are in one column.
#' @param main.title Main title for the entire plot. Default is "Posterior Densities of Group Means (mu)"
#' @param x.label X-axis label. Default is "mu"
#' @param y.label Y-axis label. Default is "Density"
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes geom_density geom_vline labs theme_minimal theme element_blank scale_color_manual
#' @importFrom grDevices colorRampPalette
#'
#' @return A ggplot object (grob from grid.arrange)
#' @export ggplot.mu.density
#' @export

ggplot.mu.density <- function(Gibbs,
                               true.val = NA,
                               add.legend = FALSE,
                               burnin = NA,
                               layout.dim = NA,
                               main.title = "Posterior Densities of Group Means (mu)",
                               x.label = "mu",
                               y.label = "Density"){

  n.iter <- dim(Gibbs$mu)[1]
  G <- dim(Gibbs$mu)[2]
  k <- dim(Gibbs$mu)[3]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  if(is.na(layout.dim[1])){
    if(add.legend){
      layout.dim <- c(k + 1, 1)  # Add extra row for legend
    } else {
      layout.dim <- c(k, 1)
    }
  }
  countOrder <- c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh")

  # Generate colors using colorRampPalette
  color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
  colors <- color_palette(G)

  turns <- seq(burnin + 1, n.iter)

  Density <- NULL
  mu.group <- rep(1:G, each = length(turns))
  mu.plots <- list()

  # Create layout matrix
  if(add.legend){
    layout.mat <- matrix(1:(k + 1), nrow = layout.dim[1], ncol = layout.dim[2], byrow = TRUE)
  } else {
    layout.mat <- matrix(1:k, nrow = layout.dim[1], ncol = layout.dim[2], byrow = TRUE)
  }

  for(i in 1:k){
    mu.all <- c()
    for(j in 1:G){
      mu.all <- c(mu.all, Gibbs$mu[turns, j, i])
    }
    mu.den <- as.data.frame(cbind(mu.all, mu.group))
    mu.den$mu.group <- as.factor(mu.den$mu.group)

    mu.plots[[i]] <- ggplot(mu.den, aes(x = mu.all, color = mu.group)) +
      geom_density(linewidth = 1) +
      labs(title = "",
           x = x.label,
           y = y.label,
           tag = paste("(", LETTERS[i], ")", sep = "")) +
      theme_minimal() +
      scale_color_manual(values = colors) +
      theme(legend.position = "none")

    if(!is.na(true.val[1])){
      mu.plots[[i]] <- mu.plots[[i]] +
        geom_vline(xintercept = true.val[, i],
                   color = colors[1:G],
                   linetype = "longdash")
    }
  }

  if(add.legend){
    levels(mu.den$mu.group) <- paste("Cluster", 1:G)
    colnames(mu.den)[2] <- "Cluster"
    mu.plot.legend <- ggplot(mu.den, aes(x = mu.all, fill = Cluster)) +
      geom_density() +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      theme(legend.position = "bottom")
    mu.plots[[(k + 1)]] <- get_legend(mu.plot.legend)
  }

  # Create the arranged plot with main title
  p <- grid.arrange(grobs = mu.plots,
                    layout_matrix = layout.mat,
                    top = main.title)

  return(p)
}

# Helper function to extract legend
get_legend <- function(plot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
