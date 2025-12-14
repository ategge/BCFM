#' @title A credible interval plot of posterior of factor loadings covariance, tau
#' @description It returns a credible interval plot of factor loadings covariance, tau. The lines are 95% intervals, while the circles are posterior mean.
#'
#' @param Gibbs Gibbs sample from \code{BCFM} function
#' @param burnin Number of burn-in period. If not specified, it uses the first tenths sample as burn-in period.
#' @param true.val True values of the taus
#' @param main.bool Return main title. Default is TRUE.
#' @param main.title Main title for the plot. Default is expression for tau.
#' @param x.label X-axis label. Default is "Factor"
#' @param y.label Y-axis label. Default is "Tau"
#' @param factor_labels Character vector of factor names. If NULL, defaults to Factor 1, Factor 2, etc.
#'
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar theme_minimal labs theme element_text
#'
#' @return A ggplot object
#' @export ggplot_tau.CI
#' @export

ggplot_tau.CI <- function(Gibbs,
                        burnin = NA,
                        true.val = NA,
                        main.bool = TRUE,
                        main.title = NULL,
                        x.label = "Factor",
                        y.label = "Tau",
                        factor_labels = NULL){

  if(is.null(Gibbs$tau2)){Gibbs$tau2 <- Gibbs$tau}
  n.iter <- dim(Gibbs$tau2)[1]
  n.factors <- dim(Gibbs$tau2)[2]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  turns <- seq(burnin + 1, n.iter)

  # Set factor labels
  if(is.null(factor_labels)){
    factor_labels <- paste0("Factor ", 1:n.factors)
  } else if(length(factor_labels) != n.factors){
    warning("Length of factor_labels does not match number of factors (k=", n.factors,
            "). Using default labels.")
    factor_labels <- paste0("Factor ", 1:n.factors)
  }

  # Calculate credible intervals
  tau2.CI <- matrix(NA, 3, n.factors)
  tau2.CI[1,] <- apply(Gibbs$tau2[turns,], 2, quantile, 0.025)
  tau2.CI[2,] <- apply(Gibbs$tau2[turns,], 2, mean)
  tau2.CI[3,] <- apply(Gibbs$tau2[turns,], 2, quantile, 0.975)

  # Prepare data
  tau_data <- data.frame(
    Factor = 1:n.factors,
    FactorLabel = factor_labels,
    Mean = tau2.CI[2,],
    Lower = tau2.CI[1,],
    Upper = tau2.CI[3,]
  )

  # Set title
  if(is.null(main.title)){
    main.title <- if(main.bool) expression(paste(tau, " 95% Credible Intervals")) else ""
  }

  # Create plot
  p <- ggplot(tau_data, aes(x = Factor, y = Mean)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
    geom_point(size = 3, color = "black") +
    theme_minimal() +
    labs(title = main.title,
         x = x.label,
         y = y.label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Add true values if available
  if(!is.na(true.val[1])){
    true_data <- data.frame(
      Factor = 1:n.factors,
      TrueValue = true.val
    )
    p <- p + geom_point(data = true_data,
                        aes(x = Factor, y = TrueValue),
                        shape = 17, size = 3.5, color = "blue", alpha = 0.5)
  }

  return(p)
}
