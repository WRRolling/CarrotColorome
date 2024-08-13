
run_BGLR_LOOCV <- function(Pheno, GD, Pop.CV, nIter, burnIn, numCores, model) {
  LOOCV <- nrow(Pheno)  # Determine the number of folds for Leave-One-Out CV
   
  # Use parallel::mclapply for parallel processing
  MPL <- parallel::mclapply(1:LOOCV, function(i) {
    PT <- Pheno
    PT$Alpha[i] <- NA  # Set the i-th observation as missing
    PT$Alpha <- as.numeric(unlist(PT$Alpha))  # Ensure numeric type for prediction

    ETA <- list(
      list(X = GD, model = model),
      list(X = Pop.CV$Q1, model = "FIXED"),
      list(X = Pheno$Interior, model = "FIXED")
    )

    fm <- BGLR(y = PT$Alpha, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
    Pred_pheno <- fm$yHat[[i]]
    Pred_pheno[[1]]  # Return the predicted value
  }, mc.cores = numCores)

  unlist(MPL)  # Combine the results from multiple cores
}
