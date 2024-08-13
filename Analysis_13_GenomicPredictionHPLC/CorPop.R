CorePop <- function() {
  Output.dat <- list()
  Output <- lapply(BGLR_settings$Traits, function(trait) {
    tic()
    MPL <- list()
    MPL <- parallel::mclapply(1:BGLR_settings$LOOCV, function(j) {
      # Make an object where you can set one phenotype to NA
      PT <- Pheno
      # Set the jth phenotype to NA for the specified trait
      PT[[trait]][j] <- NA  # Set the jth observation as missing for the chosen trait
      # Ensure the phenotype is numeric
      PT[[trait]] <- as.numeric(unlist(PT[[trait]]))  # Ensure numeric type for prediction
      # Set the model options
      ETA <- list(list(X=Pheno$Interior, model="FIXED"),
                  list(X=Qmat$Q1, model="FIXED"))
      # Run Genomic Prediction
      fm <- BGLR(y = PT[[trait]], ETA = ETA,
                 nIter = BGLR_settings$nIter, burnIn = BGLR_settings$burnIn,
                 verbose = FALSE)
      # Extract Predicted value
      Pred_pheno <- fm$yHat[[j]]
      # Return the predicted value
      MPL[j] <- Pred_pheno[[1]]
    }, mc.cores = BGLR_settings$numCores)
    Commune.Check(MPL, ETA)
    dat <- as.numeric(unlist(MPL))
    Out <- calculate_metrics(dat, Pheno, trait)
    timreg <- toc()
    Out[6] <- timreg$toc - timreg$tic
    return(Out)
    print(Out)
  })
  Output.dat <- do.call(rbind, Output)  # Combine data frames from all iterations
  Output.dat <<- Output.dat
}
