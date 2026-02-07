#' @title Initialize Data Array for BCFM Model
#' @description Prepares the input data by converting it into a 3D array format
#' required by the BCFM model. If data is already in the correct 3D array format,
#' it returns the data as-is. Takes selected clustering variables and creates
#' an array with dimensions (observations, variables, time points).
#'
#' @param data A data frame, matrix, or 3D array containing the data to be used
#'   for clustering. If a 3D array with appropriate dimensions is provided and
#'   cluster.vars is NULL, the function returns the data unchanged.
#' @param cluster.vars A character vector specifying the column names of variables
#'   to be used for clustering (required for data frames). If NULL and data is
#'   already a 3D array, the function returns data as-is. If NULL and data is a
#'   matrix, all columns are used.
#'
#' @return A 3D array with dimensions (n, p, t) where n is the number of
#'   observations, p is the number of clustering variables, and t is the number
#'   of time points (defaults to 1 for cross-sectional data).
#'
#' @examples
#' \donttest{
#' # Example 1: Data frame with variable selection
#' data <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
#' cluster.vars <- c("x1", "x2", "x3")
#' data.pre <- init.data(data, cluster.vars)
#'
#' # Example 2: Matrix (uses all columns)
#' data_matrix <- matrix(rnorm(300), nrow = 100, ncol = 3)
#' data.pre <- init.data(data_matrix)
#'
#' # Example 3: 3D array (returns as-is)
#' data_3d <- array(rnorm(1500), dim = c(100, 3, 5))
#' data.pre <- init.data(data_3d)  # Returns unchanged
#' }
#'
#' @export
init.data <- function(data, cluster.vars = NULL) {

  # Case 1: Data is already a 3D array
  if (is.array(data) && length(dim(data)) == 3) {
    if (is.null(cluster.vars)) {
      # Data is already in correct format, return as-is
      return(data)
    } else {
      # User wants to subset variables from 3D array
      if (is.numeric(cluster.vars)) {
        # Numeric indices
        data.pre <- data[, cluster.vars, , drop = FALSE]
      } else {
        # Character names - not typical for arrays but handle it
        data.pre <- data
      }
      return(data.pre)
    }
  }

  # Case 2: Data is a matrix or 2D array
  if (is.matrix(data) || (is.array(data) && length(dim(data)) == 2)) {
    if (is.null(cluster.vars)) {
      # Use all columns
      n_obs <- nrow(data)
      n_vars <- ncol(data)
      data.pre <- array(dim = c(n_obs, n_vars, 1))
      data.pre[, , 1] <- data
      return(data.pre)
    } else {
      # Subset columns
      if (is.numeric(cluster.vars)) {
        selected_data <- data[, cluster.vars, drop = FALSE]
      } else if (is.character(cluster.vars)) {
        if (is.null(colnames(data))) {
          stop("cluster.vars provided as names but data has no column names")
        }
        selected_data <- data[, cluster.vars, drop = FALSE]
      } else {
        stop("cluster.vars must be numeric indices or character names")
      }

      n_obs <- nrow(selected_data)
      n_vars <- ncol(selected_data)
      data.pre <- array(dim = c(n_obs, n_vars, 1))
      data.pre[, , 1] <- selected_data
      return(data.pre)
    }
  }

  # Case 3: Data is a data frame (original behavior)
  if (is.data.frame(data)) {
    if (is.null(cluster.vars)) {
      # Use all numeric columns
      numeric_cols <- sapply(data, is.numeric)
      if (sum(numeric_cols) == 0) {
        stop("No numeric columns found in data frame")
      }
      cluster.vars <- names(data)[numeric_cols]
      message(sprintf("No cluster.vars specified. Using all %d numeric columns.",
                      length(cluster.vars)))
    }

    # Check if cluster.vars exist in data
    missing_vars <- setdiff(cluster.vars, names(data))
    if (length(missing_vars) > 0) {
      stop(sprintf("Variables not found in data: %s",
                   paste(missing_vars, collapse = ", ")))
    }

    # Prep the data
    n_obs <- nrow(data)
    n_vars <- length(cluster.vars)
    data.pre <- array(dim = c(n_obs, n_vars, 1))
    data.pre[, , 1] <- as.matrix(data[cluster.vars])
    return(data.pre)
  }

  # Case 4: Unsupported data type
  stop(sprintf("Unsupported data type: %s. Expected data.frame, matrix, or 3D array.",
               class(data)[1]))
}
