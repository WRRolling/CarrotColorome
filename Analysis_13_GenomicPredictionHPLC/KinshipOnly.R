KinshipOnly <- function() {
  for (i in 1:length(Kinships.Files)){
    Kin <- read.csv(paste0(genotype.path, args[2], Kinships.Files[i]), head=F)
    Kin <- Kin[which(Kin$V1 %in% Pheno$Taxa),
               (1+which(Kin$V1 %in% Pheno$Taxa))]
    Kin <- as.matrix(Kin)
    print(dim(Kin))	    
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
        ETA <- list(list(K = Kin, model = "RKHS"))
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

