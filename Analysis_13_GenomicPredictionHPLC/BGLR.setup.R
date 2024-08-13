BGLRSetup <- function(LOOCV = nrow(Pheno)) {
  # Assign values to variables within the function
  nIter <- 20000
  burnIn <- 5000
  numCores <- parallel::detectCores() - 1  # Use parallel:: for clarity
  LOOCV <- LOOCV  # Initialize as NULL to avoid errors
  Models <- c("BayesA", "BayesB", "BayesC", "BL", "BRR", "RKHS")
  Traits <- c("Total", "Lutein", "Alpha", "Beta", "Lycopene", "Phytoene", "Zeta")

  # Return a list containing all the values
  return(list(nIter = nIter, burnIn = burnIn, numCores = numCores,
               LOOCV = LOOCV, Models = Models, Traits = Traits))
}

