
# BCFM

<!-- badges: start -->

<!-- badges: end -->

Bayesian Clustering Factor Model (BCFM) for clustering and latent factor
analysis of multivariate cross-sectional data.

## Installation

You can install BCFM from CRAN with:

``` r
install.packages("BCFM")
```

You can also install the development version of BCFM from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ategge/BCFM", build_vignettes = TRUE)
```

## Example

This is a basic example which shows you how to use BCFM:

``` r
library(BCFM)

# Load example data
data("sim.data", package = "BCFM")

# Specify variables to use for clustering
cluster.vars <- paste0("Var", 1:20)

# Create output directory for results
output_dir <- "~/BCFM_results"

# Run model selection
BCFM.model.selection(
  data = sim.data$data,
  cluster.vars = cluster.vars,   # Required parameter
  grouplist = 2:4,               # Try 2, 3, and 4 groups
  factorlist = 2:4,              # Try 2, 3, and 4 factors
  n.iter = 10000,                # Number of MCMC iterations
  burnin = 5000,                 # Burnin for Information Criterion calculations
  every = 1000                   # Progress update frequency
  output_dir = output_dir        # Specify where to save results

)

# Results are saved in output_dir
# Load and visualize BIC results
load(file.path(output_dir, "BIC.Rdata"))
ggplot_BIC(BIC.matrix, factor_list = 2:4, group_list = 2:4)

# Load and visualize model results for 4 groups and 4 factors
load(file.path(output_dir, "results-covarianceF-g4-f4.Rdata"))
ggplot_latent.profiles(SDresult$Result)
```

## Vignette

For a complete workflow tutorial, see the vignette:

``` r
# After installation with build_vignettes = TRUE
vignette("introduction-to-BCFM", package = "BCFM")

# Or browse all vignettes
browseVignettes("BCFM")
```

## Features

- Bayesian clustering with latent factor models
- Model selection across different numbers of groups and factors
- Comprehensive visualization functions
- Example datasets included
- Detailed vignette with complete workflow

## Citation

If you use BCFM in your research, please cite:

    [Add your citation here when you have a publication]

## License

GPL-3
