#' Plot IC Matrix from Model Selection
#'
#' Creates two plots showing IC values across different numbers
#' of groups and factors to help identify the optimal model configuration.
#'
#' @param matrix An IC matrix from BCFM.model.selection, with rows representing
#'   groups and columns representing factors.
#' @param factor_list Numeric vector of factor values corresponding to matrix columns.
#'   Default is 2:6.
#' @param group_list Numeric vector of group values corresponding to matrix rows.
#'   Default is 5:11.
#' @param combine Logical. If TRUE (default), returns a combined plot using ggarrange.
#'   If FALSE, returns a list of two separate ggplot objects.
#'
#' @return If combine = TRUE, a combined ggplot object. If combine = FALSE, 
#'   a list with two elements:
#' \describe{
#'   \item{by_groups}{ggplot showing IC vs. number of groups}
#'   \item{by_factors}{ggplot showing IC vs. number of factors}
#' }
#'
#' @details The function creates two complementary visualizations:
#' \itemize{
#'   \item Plot 1: IC vs. number of groups, with lines for each number of factors
#'   \item Plot 2: IC vs. number of factors, with lines for each number of groups
#' }
#'
#' Lower IC values indicate better model fit.
#'
#' @importFrom dplyr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_bw guides guide_legend
#' @importFrom ggpubr ggarrange
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After running BCFM.model.selection
#' load("IC.Rdata")
#' 
#' # Combined plot (default)
#' plot_IC(IC.matrix, factor_list = 2:4, group_list = 2:4)
#' 
#' # Separate ggplot objects
#' plots <- plot_IC(IC.matrix, factor_list = 2:4, group_list = 2:4, 
#'                   combine = FALSE)
#' plots$by_groups      # First plot
#' plots$by_factors     # Second plot
#' 
#' # Customize individual plots
#' plots$by_groups + ggplot2::ggtitle("My Custom Title")
#' }
ggplot_IC <- function(matrix, factor_list = 2:4, group_list = 2:4, combine = TRUE) {
  
  # Groups on X axis
  rows2 <- data.frame(matrix)
  colnames(rows2) <- factor_list
  
  f <- rows2 %>%
    rownames_to_column('groups') %>%
    pivot_longer(cols = -groups, names_to = 'factors')
  
  f$factors <- factor(f$factors, levels = factor_list)
  f$groups <- factor(f$groups, levels = group_list)
  
  plot1 <- f %>%
    ggplot() + 
    aes(groups, value, color = as.factor(factors), group = as.factor(factors)) + 
    geom_point(size = 2) + 
    geom_line() +
    labs(title = "IC by Number of Groups",
         x = "Number of Groups",
         y = "IC") +
    theme_bw() + 
    guides(color = guide_legend(title = "Factors"))
  
  # Factors on X axis
  plot2 <- f %>%
    ggplot() + 
    aes(factors, value, color = as.factor(groups), group = as.factor(groups)) + 
    geom_point(size = 2) + 
    geom_line() +
    labs(title = "IC by Number of Factors",
         x = "Number of Factors",
         y = "IC") +
    theme_bw() + 
    guides(color = guide_legend(title = "Groups"))
  
  if (combine) {
    return(ggarrange(plot1, plot2, widths = c(1.2, 1.4)))
  } else {
    return(list(by_groups = plot1, by_factors = plot2))
  }
}