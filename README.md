
# BCFM

<!-- badges: start -->
<!-- badges: end -->

Bayesian Covariance Factor Model (BCFM) for clustering and latent factor
analysis of multivariate longitudinal data.

## Installation

You can install BCFM from CRAN with:

``` r
install.packages("BCFM")
```

You can also install the development version of BCFM from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ategge/BCFM")
```

## Example

This is a basic example which shows you how to use BCFM:

``` r
library(BCFM)

# Load example data
data("sim.data", package = "BCFM")

# Run model selection
BCFM.model.selection(
  data = sim.data$data,
  grouplist = 2:4,        # Try 2, 3, and 4 groups
  factorlist = 2:4,       # Try 2, 3, and 4 factors
  n.iter = 10000,         # Number of MCMC iterations
  burnin = 5000,          # Burnin for BIC calculation
  every = 1000            # Progress update frequency
)

# Load and visualize results
load("results-covarianceF-g4-f4.Rdata")
ggplot.latent.profiles(SDresult$Result)
```

## Vignette

For a complete workflow tutorial, see the vignette:

``` r
# After installation
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
