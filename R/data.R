#' Example BCFM Results
#'
#' Pre-computed results from fitting a BCFM model with 4 groups and 4 factors
#' on the simulated dataset. These results are used in the package vignette
#' to demonstrate visualization and interpretation.
#'
#' @format A list with two components:
#' \describe{
#'   \item{Result}{A list containing the MCMC samples and fitted model components:
#'     B (loading matrix), mu (cluster means), X (factor scores), Z (cluster assignments),
#'     probs (cluster probabilities), sigma (error variances), tau (factor variances),
#'     and omega (within-cluster covariances)}
#'   \item{order}{Variable ordering used in the model}
#' }
#'
#' @details
#' This object was created by running BCFM.model.selection with 4 groups and 4 factors
#' on the simulated data for 10,000 MCMC iterations with 5,000 burnin.
#'
#' @source Generated using the BCFM package on simulated data
#'
#' @examples
#' # Load the example results
#' data("SDresult", package = "BCFM")
#'
#' # Check structure
#' str(SDresult$Result, max.level = 1)
#'
#' # Create visualizations
#' ggplot.latent.profiles(SDresult$Result)
#'
"SDresult"
