## SPEAR - Sparse Supervised Bayesian Factor Model for Multi-Omic Analysis

### Package Requirements

The following packages are required for installing SPEAR and the vignettes. Future versions of SPEAR will install these automatically. All of them can be installed via the `install.packages` function.

**SPEAR Installation:** `Rcpp`, `RcppArmadillo`, `BH`, `glmnet`, `parallel`

**Vignettes:** `jsonlite`, `dplyr`, `reshape2`, `stringr`, `ggplot2`, `cowplot`, `ggthemes`, `maps`

### Installation

To install in R:

`remotes::install_github("jgygi/SPEAR@main", build_vignettes = TRUE)`

**NOTE**: this will install the vignettes described below as well.

### Vignettes

Once installed, use this code to see all vignettes:

`browseVignettes("SPEARcomplete")`

If you want to open vignettes directly from RStudio, run:

`vignette("SPEAR_simulate_data")` - used to simulate multi-omic data for use in the other vignettes

`vignette("SPEAR_running_spear")` - instructions on running SPEAR on multi-omic data

`vignette("SPEAR_downstream_vignette")` - an example of a typical downstream analysis of a SPEAR model
