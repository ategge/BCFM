#' @title A credible interval plot of posterior of sigma squared
#' @description It returns a credible interval plot of idiosyncratic variance, sigma squared. The lines are 95% intervals, while the circles are posterior mean.
#'
#' @param Gibbs Gibbs sample from \code{BCFM} function
#' @param burnin Number of burn-in period. If not specified, it uses the first tenths sample as burn-in period.
#' @param permutation Permutation vector, if applicable
#' @param main.bool Return main title. Default is TRUE.
#' @param main.title Main title for the plot. Default is expression for sigma squared.
#' @param x.label X-axis label. Default is "Variables"
#' @param y.label Y-axis label. Default is "Variance"
#' @param var_labels Character vector of variable names. If NULL, defaults to Variable 1, Variable 2, etc.
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar theme_minimal labs theme element_text
#'
#' @export ggplot.sigma2.CI
#' @export

ggplot.sigma2.CI <- function(Gibbs,
                           burnin = NA,
                           permutation = NA,
                           main.bool = TRUE,
                           main.title = NULL,
                           x.label = "Variables",
                           y.label = "Variance",
                           var_labels = NULL){

  n.iter <- dim(Gibbs$sigma2)[1]
  r <- dim(Gibbs$sigma2)[2]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  turns <- seq(burnin + 1, n.iter)
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
  sigma2.CI <- matrix(NA, 3, r)
  sigma2.CI[1,] <- apply(Gibbs$sigma2[turns, repermute], 2, quantile, 0.025)
  sigma2.CI[2,] <- apply(Gibbs$sigma2[turns, repermute], 2, mean)
  sigma2.CI[3,] <- apply(Gibbs$sigma2[turns, repermute], 2, quantile, 0.975)

  # Prepare data
  sigma2_data <- data.frame(
    Variable = 1:r,
    VarLabel = var_labels,
    Mean = sigma2.CI[2,],
    Lower = sigma2.CI[1,],
    Upper = sigma2.CI[3,]
  )

  # Set title
  if(is.null(main.title)){
    main.title <- if(main.bool) expression(paste(sigma^2, " 95% Credible Intervals")) else ""
  }

  # Create plot
  p <- ggplot(sigma2_data, aes(x = Variable, y = Mean)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
    geom_point(size = 3, color = "black") +
    theme_minimal() +
    labs(title = main.title,
         x = x.label,
         y = y.label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
