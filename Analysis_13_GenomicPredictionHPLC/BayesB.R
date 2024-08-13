# Each Model Own Script!
# Once Completed Find Best (Accuracy or Time)

BayesB.Mod <- function() {
  for (i in 1:length(Genotype.Files)){
    GD <- vroom(file=paste0(genotype.path, args[2], Genotype.Files[i]), col_names=F)
    GD <- GD[-1,]
    GD <- GD[which(GD$X1 %in% Pheno$Taxa),]
    GD <- GD[order(GD$X1),]
    GD <- GD[,-1]
    GD <- sapply(GD, as.numeric)
    
    Output <- lapply(BGLR_settings$Traits, function(trait) {
      tic()
      MPL <- parallel::mclapply(1:BGLR_settings$LOOCV, function(j) {
        # Make an object where you can set one phenotype to NA
        PT <- Pheno
        # Set the jth phenotype to NA for the specified trait
        PT[[trait]][j] <- NA  # Set the jth observation as missing for the chosen trait
        # Ensure the phenotype is numeric
        PT[[trait]] <- as.numeric(unlist(PT[[trait]]))  # Ensure numeric type for prediction
        # Set the model options
        ETA <- list(list(X = GD, model = "BayesB"))
        # Run Genomic Prediction
        fm <- BGLR(y = PT[[trait]], ETA = ETA,
                   nIter = BGLR_settings$nIter, burnIn = BGLR_settings$burnIn,
                   verbose = FALSE)
        # Extract Predicted value
        Pred_pheno <- fm$yHat[[j]]
        # Return the predicted value
        return(Pred_pheno[[1]])
      }, mc.cores = BGLR_settings$numCores)
      Commune.Check(MPL, ETA)
      dat <- as.numeric(unlist(MPL))
      Out <- calculate_metrics(dat, Pheno, trait)
      timreg <- toc()
      Out[6] <- timreg$toc - timreg$tic
      return(Out)
    })
    Output.2 <- do.call(rbind.data.frame, Output)
    rownames(Output.2) <- BGLR_settings$Traits
    colnames(Output.2) <- c("Cor","Pval", "H2", "Top25", "Bot50", "Time")
    if (i < 2) {
      Output.dat <- Output.2
    } else {
      Output.dat <- rbind(Output.dat, Output.2)
    }
  }
  Output.dat <<- Output.dat
}
