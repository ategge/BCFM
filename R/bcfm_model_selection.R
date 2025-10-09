#' @title BCFM Model Selection Over Multiple Groups and Factors
#' @description Performs Bayesian Covariance Factor Model analysis across a grid
#' of group numbers and factor numbers. For each combination, the function fits
#' the BCFM model, calculates BIC, and saves results. This is the primary function
#' for model selection to determine the optimal number of clusters and latent factors.
#'
#' @param data A data frame containing the variables to be analyzed
#' @param cluster.vars A character vector specifying the column names of variables
#'   to be used for clustering
#' @param grouplist A numeric vector specifying the numbers of groups to test
#'   (e.g., c(2, 3, 4, 5))
#' @param factorlist A numeric vector specifying the numbers of latent factors
#'   to test (e.g., c(1, 2, 3))
#' @param n.iter Number of MCMC iterations for each model. Default is 50000.
#' @param vague.mu Logical indicating whether to use vague priors for mu.
#'   Default is FALSE.
#' @param covariance Logical indicating whether to model covariance structure.
#'   Default is TRUE.
#' @param p.exponent The Dirichlet priors exponent for probabilities. Default is 2.
#' @param every Integer specifying the frequency of progress updates during MCMC.
#'   Default is 10.
#' @param burnin Number of initial MCMC iterations to discard when calculating BIC.
#'   If NA, an appropriate burnin is determined automatically.
#'
#' @details The function performs the following steps for each group-factor combination:
#' \enumerate{
#'   \item Preprocesses data using init.data
#'   \item Determines optimal variable ordering using permutation.order
#'   \item Initializes model attributes and hyperparameters
#'   \item Fits BCFM model using BCFM.fit
#'   \item Calculates BIC for model comparison
#'   \item Saves individual results and cumulative BIC matrix
#' }
#'
#' The BIC matrix can be used to identify the optimal model configuration by
#' selecting the combination of groups and factors with the lowest BIC value.
#'
#' @return Invisibly returns NULL. Results are saved to disk:
#' \describe{
#'   \item{results-covarianceF-gX-fY.Rdata}{Individual model results for
#'     each group-factor combination (where X is the number of groups and Y is
#'     the number of factors), containing SDresult and variable order}
#'   \item{BIC.Rdata}{Contains BIC.matrix, md.matrix, timing information, and data.
#'     Load this file to compare models and identify the optimal configuration.}
#' }
#'
#' @note This function can be computationally intensive as it fits multiple models.
#' Consider running on high-performance computing resources for large datasets or
#' extensive model grids. The function includes error handling to continue execution
#' even if individual models fail to converge.
#'
#' @examples
#' \dontrun{
#' # Run model selection over 2-4 groups and 1-3 factors
#' BCFM.model.selection(
#'   data = mydata,
#'   cluster.vars = c("var1", "var2", "var3", "var4"),
#'   grouplist = 2:4,
#'   factorlist = 1:3,
#'   n.iter = 20000,
#'   every = 100,
#'   burnin = 5000
#' )
#'
#' # Load and examine BIC results
#' load("BIC.Rdata")
#' print(BIC.matrix)
#'
#' # Find optimal model
#' optimal <- which(BIC.matrix == min(BIC.matrix, na.rm = TRUE), arr.ind = TRUE)
#' cat("Optimal groups:", optimal[1], "Optimal factors:", optimal[2], "\n")
#' }
#'
#' @seealso \code{\link{BCFM.fit}} for fitting a single model,
#'   \code{\link{init.data}}, \code{\link{initialize.model.attributes}},
#'   \code{\link{initialize.hyp.parm}}
#'
#' @export
BCFM.model.selection <- function(data, cluster.vars, grouplist, factorlist,
                                 n.iter = 50000, vague.mu = FALSE,
                                 covariance = TRUE, p.exponent = 2,
                                 every = 10, burnin = NA) {

  data.pre <- init.data(data, cluster.vars)


  BIC.matrix = matrix(NA, nrow=length(grouplist), ncol=length(factorlist))
  rownames(BIC.matrix) <- paste0("G", grouplist)
  colnames(BIC.matrix) <- paste0("F", factorlist)

  # Compare the run time with loops
  tic2 <- Sys.time()

  for (i in seq_along(grouplist)) {
    ng <- grouplist[i]

    for (j in seq_along(factorlist)) {
      nf <- factorlist[j]
      message(sprintf("Factor: %d, Group: %d", nf, ng))

      # For reproducibility: Setting seed for each dataset/number of groups/number of factors
      set.seed(10 * ng + nf)

      filename.temp.out <- sprintf("results-covarianceF-g%d-f%d.Rdata", ng, nf)

      # Order variables
      order <- permutation.order(data.pre, covariance = covariance, L = nf, fa = TRUE)
      permutation <- permutation.scale(data.pre, permutation = order,
                                       covariance = covariance,
                                       return.array = TRUE, num.layers = 1)
      data <- permutation

      # Initialize model components
      model.attributes <- initialize.model.attributes(
        S = dim(data)[1],
        times = 1,
        R = dim(data)[2],
        L = nf,
        G = ng
      )

      cluster.hyperparms <- initialize.cluster.hyperparms(
        data, model.attributes,
        covariance = covariance,
        diag.Psi = FALSE,
        vague.mu = FALSE
      )

      hyp.parm <- initialize.hyp.parm(
        model.attributes, cluster.hyperparms,
        n.sigma = 2.2,
        n.s2.sigma = 0.1,
        n.tau = 1,
        n.s2.tau = 1,
        p.exponent = 2
      )

      # Fit model with error handling
      tryCatch(
        {
          SDresult <- BCFM.fit(data, model.attributes, hyp.parm,
                               n.iter = n.iter, vague.mu = vague.mu,
                               covariance = covariance, p.exponent = p.exponent,
                               every = every)
          save(SDresult, order, file = filename.temp.out)
            # Compare the run time with loops
          BIC.matrix[i,j] = BIC.like(data, SDresult$Result,
                                     model.attributes, burnin = NA)

        },
        error = function(cond) {
          message("Error message:")
          message(conditionMessage(cond))
          NA
        },
        warning = function(cond) {
          message("Warning message:")
          message(conditionMessage(cond))
          NULL
        },
        finally = {
          message("Done")
        }
      )
    }

    toc2 <- Sys.time()
    run.time2 <- toc2 - tic2
    message(sprintf("Time: %s", format(run.time2)))
    save(BIC.matrix, tic2, toc2, run.time2, data, file = "BIC.Rdata")
  }

  invisible(NULL)
}
