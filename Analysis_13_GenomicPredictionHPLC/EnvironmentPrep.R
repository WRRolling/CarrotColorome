prepare_environment_and_data <- function() {
  # Clear the global environment
  rm(list = ls())

  # List required packages
  required_packages <- c("BGLR", "bigmemory", "data.table", "biganalytics",
                         "tidyverse", "parallel", "MASS", "vroom", "tictoc",
			"matrixStats")

  # Install and load packages
  install_and_load_packages <- function(package_name) {
    if (!require(package_name, character.only = TRUE)) {
      install.packages(package_name, dependencies = TRUE)
    }
    library(package_name, character.only = TRUE)
  }

  lapply(required_packages, install_and_load_packages)

  # Set CRAN mirror (optional)
  chooseCRANmirror(ind = 71)  # Adjust mirror index if needed
}


