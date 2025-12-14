#' Simulated BCFM Dataset
#'
#' A simulated dataset for demonstrating the BCFM package functionality.
#'
#' @format A list with components:
#' \describe{
#'   \item{data}{A 3D numeric array (200 x 20 x 5)}
#'   \item{pselected}{Cluster assignments}
#'   \item{X}{Factor scores}
#'   \item{B}{Loading matrix}
#' }
#'
#' @examples
#' data("sim.data", package = "BCFM")
#' dim(sim.data$data)
#'
"sim.data"