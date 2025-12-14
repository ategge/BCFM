#' @title Build factor loadings plot
#' @description The function builds a column-wise plots of factor loadings. The parameters fixed at 1 are displayed with red dashed vertical lines.
#'
#' @param Gibbs Result of Gibbs sampler from BCFM function
#' @param true.val True values of factor loadings. If not available, NA.
#' @param burnin Number of burn-in. If not set, it uses the first tenths as burn-in period.
#' @param permutation Permutation of variables. If not set, no permutation.
#' @param main.bool Add title of the plots. Default is TRUE.
#' @param var_labels Character vector of variable names. If NULL, defaults to Variable 1, Variable 2, etc.
#'
#' @importFrom stats quantile
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline geom_vline facet_wrap theme_minimal labs theme element_text
#'
#' @return A ggplot object
#' @export ggplot_B.CI
#' @export

ggplot_B.CI <- function(Gibbs,
                      true.val = NA,
                      burnin = NA,
                      permutation = NA,
                      main.bool = TRUE,
                      var_labels = NULL){

  # Setup
  if(is.na(burnin)){burnin <- round(dim(Gibbs$B)[1] / 10)}
  n.iter <- dim(Gibbs$B)[1]
  turns <- seq(burnin + 1, n.iter)
  r <- dim(Gibbs$B)[2]
  k <- dim(Gibbs$B)[3]
  countOrder <- c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh")
  repermute <- 1:r
  if(!is.na(permutation[1])){repermute <- order(permutation)}

  # Set variable labels
  if(is.null(var_labels)){
    var_labels <- paste0("Variable ", 1:r)
  } else if(length(var_labels) != r){
    warning("Length of var_labels does not match number of variables (r=", r,
            "). Using default labels.")
    var_labels <- paste0("Variable ", 1:r)
  }

  # Calculate credible intervals
  B.Gibbs.CI <- array(NA, dim = c(r, k, 3))
  B.Gibbs.CI[,,1] <- apply(Gibbs$B[turns, repermute,], c(2, 3), quantile, p = 0.025)
  B.Gibbs.CI[,,2] <- apply(Gibbs$B[turns, repermute,], c(2, 3), mean)
  B.Gibbs.CI[,,3] <- apply(Gibbs$B[turns, repermute,], c(2, 3), quantile, p = 0.975)

  # Prepare data for ggplot
  plot_data <- data.frame()
  true_data <- data.frame()
  vline_data <- data.frame()

  for(i in 1:k){
    factor_name <- if(main.bool) paste(countOrder[i], "factor") else paste("Factor", i)

    # Main CI data
    temp_df <- data.frame(
      Variable = 1:r,
      VarLabel = var_labels,
      Mean = B.Gibbs.CI[, i, 2],
      Lower = B.Gibbs.CI[, i, 1],
      Upper = B.Gibbs.CI[, i, 3],
      Factor = factor_name
    )
    plot_data <- rbind(plot_data, temp_df)

    # True values if available
    if(!is.na(true.val[1])){
      true_df <- data.frame(
        Variable = 1:r,
        VarLabel = var_labels,
        TrueValue = true.val[, i],
        Factor = factor_name
      )
      true_data <- rbind(true_data, true_df)
    }

    # Vertical lines for fixed parameters
    if(!is.na(permutation[1])){
      vline_df <- data.frame(
        xintercept = permutation[i],
        Factor = factor_name
      )
      vline_data <- rbind(vline_data, vline_df)
    }
  }

  # Create factor ordering
  factor_levels <- if(main.bool){
    paste(countOrder[1:k], "factor")
  } else {
    paste("Factor", 1:k)
  }
  plot_data$Factor <- factor(plot_data$Factor, levels = factor_levels)
  if(nrow(true_data) > 0){
    true_data$Factor <- factor(true_data$Factor, levels = factor_levels)
  }
  if(nrow(vline_data) > 0){
    vline_data$Factor <- factor(vline_data$Factor, levels = factor_levels)
  }

  # Build the plot
  p <- ggplot(plot_data, aes(x = Variable, y = Mean)) +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
    geom_point(size = 2.5, color = "black") +
    facet_wrap(~ Factor, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(x = "Variables",
         y = "Factor Loading",
         title = if(main.bool) "Factor Loadings with 95% Credible Intervals" else NULL) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  # Add true values if available
  if(nrow(true_data) > 0){
    p <- p + geom_point(data = true_data,
                        aes(x = Variable, y = TrueValue),
                        shape = 17, size = 3, color = "blue", alpha = 0.5)
  }

  # Add vertical lines for fixed parameters
  if(nrow(vline_data) > 0){
    p <- p + geom_vline(data = vline_data,
                        aes(xintercept = xintercept),
                        color = "red", linetype = "dashed", linewidth = 0.5)
  }

  return(p)
}
