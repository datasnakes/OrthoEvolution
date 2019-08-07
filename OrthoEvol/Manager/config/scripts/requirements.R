#!/usr/local/apps/R/R-3.4.4/bin/ Rscript
# List Dependencies
required_cran_packages <- c(
  "tidyverse",
  "dplyr",
  "lubridate",
  "ggplot2",
  "tidyr",
  "ggrepel",
  "optparse"
)

# Install needed packages --------------------------------------------------------------------
# Create a list of packages to install.
cran_packages_to_install <- required_cran_packages[!(required_cran_packages %in% installed.packages()[, 1])]

# Install packages from lists.
if (length(cran_packages_to_install) > 0) {
  install.packages(cran_packages_to_install, repos = "https://cran.rstudio.com")
}