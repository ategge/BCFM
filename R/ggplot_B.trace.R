#' @title Trace plot for posterior of factor loadings
#' @description It returns a trace plot of factor loadings, B, showing MCMC convergence
#'
#' @param Gibbs MCMC sample simulated from \code{BCFMR} or \code{BCFMcpp} function
#' @param burnin Number of burn-in period. If not specified, no burn-in is removed.
#' @param permutation Permutation order vector, if applicable
#' @param true.val True values of factor loadings. If not available, NA.
#' @param factor.num The index of variable to plot (which variable's loadings across all factors to display). Use NULL to plot all variables.
#' @param var_labels Character vector of variable names. If NULL, uses Variable 1, Variable 2, etc.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual theme_minimal labs theme element_text
#' @importFrom grDevices colorRampPalette
#'
#' @return A ggplot object showing trace plots
#' @export ggplot_B.trace
#' @export
ggplot_B.trace <- function(Gibbs,
                           burnin = NA,
                           permutation = NA,
                           true.val = NA,
                           factor.num = 1,
                           var_labels = NULL){

  if(is.na(burnin)){burnin <- 0}
  n.iter <- dim(Gibbs$B)[1]
  r <- dim(Gibbs$B)[2]
  k <- dim(Gibbs$B)[3]

  repermute <- 1:r
  if(!is.na(permutation[1])){repermute <- order(permutation)}

  # Remove burn-in
  Gibbs$B <- Gibbs$B[(burnin + 1):n.iter, repermute, ]

  # Set variable labels
  if(is.null(var_labels)){
    var_labels <- paste0("Variable ", 1:r)
  } else if(length(var_labels) != r){
    warning("Length of var_labels does not match number of variables (r=", r,
            "). Using default labels.")
    var_labels <- paste0("Variable ", 1:r)
  }

  # Check if plotting all variables or just one
  if(is.null(factor.num)){
    # Plot all variables
    B.trace_list <- list()
    for(i in 1:r){
      B.temp <- Gibbs$B[, i, ]
      if(k == 1){
        B.temp <- matrix(B.temp, ncol = 1)
      }
      B.temp <- as.data.frame(B.temp)
      colnames(B.temp) <- paste0("Factor ", 1:k)
      B.temp$Iteration <- 1:nrow(B.temp)
      B.temp$Variable <- var_labels[i]
      B.trace_list[[i]] <- B.temp
    }
    B.trace <- do.call(rbind, B.trace_list)

    # Convert to long format
    B.trace <- pivot_longer(B.trace,
                            -c(Iteration, Variable),
                            values_to = "Value",
                            names_to = "Factor")

    # Set factor levels
    B.trace$Factor <- factor(B.trace$Factor, levels = paste0("Factor ", 1:k))
    B.trace$Variable <- factor(B.trace$Variable, levels = var_labels)

    # Generate colors
    color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
    colors <- color_palette(k)
    names(colors) <- paste0("Factor ", 1:k)

    # Create plot with facets
    p <- ggplot(data = B.trace, aes(x = Iteration, y = Value, color = Factor)) +
      geom_line(alpha = 0.8, linewidth = 0.8) +
      scale_color_manual(values = colors) +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      theme_minimal() +
      labs(title = "Trace Plots: All Variables",
           x = "BCFM Iteration (post burn-in)",
           y = "Loading Value",
           color = "Factor") +
      theme(legend.position = "right",
            strip.text = element_text(face = "bold"))

  } else {
    # Plot single variable
    B.trace <- Gibbs$B[, factor.num, ]

    # Handle case where k=1 (single factor)
    if(k == 1){
      B.trace <- matrix(B.trace, ncol = 1)
    }

    B.trace <- as.data.frame(B.trace)
    colnames(B.trace) <- paste0("Factor ", 1:k)
    B.trace$Iteration <- 1:nrow(B.trace)

    # Convert to long format
    B.trace <- pivot_longer(B.trace,
                            -c(Iteration),
                            values_to = "Value",
                            names_to = "Factor")

    # Set factor levels to maintain order
    B.trace$Factor <- factor(B.trace$Factor, levels = paste0("Factor ", 1:k))

    # Generate colors
    color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))
    colors <- color_palette(k)
    names(colors) <- paste0("Factor ", 1:k)

    # Create title
    var_name <- if(factor.num <= length(var_labels)){
      var_labels[factor.num]
    } else {
      paste("Variable", factor.num)
    }

    # Create plot
    p <- ggplot(data = B.trace, aes(x = Iteration, y = Value, color = Factor)) +
      geom_line(alpha = 0.8, linewidth = 0.8) +
      scale_color_manual(values = colors) +
      theme_minimal() +
      labs(title = paste("Trace Plot:", var_name),
           x = "BCFM Iteration (post burn-in)",
           y = "Loading Value",
           color = "Factor") +
      theme(legend.position = "right")
  }

  return(p)
}
