# Declare global variables to avoid R CMD check notes
# These are column names created by tidyverse operations
utils::globalVariables(c(
  # ggplot_B.CI
  "Variable", "Mean", "Lower", "Upper", "TrueValue", "xintercept",
  # ggplot_B.trace  
  "Iteration", "Value", "Factor",
  # ggplot_IC (new function)
  "groups", "value", "factors",
  # ggplot_Zit.heatmap
  "group", "assigned", "subject", "x", "y",
  # ggplot_latent.profiles
  "Posterior_Mean", "cluster_name",
  # ggplot_mu.density
  "Cluster",
  # ggplot_omega.density
  "Pair",
  # ggplot_probs
  "Probability",
  # ggplot_variability
  "Proportion"
))
