Modsummary <- function() {
  
  Mod.dat.list <- c("BayesA.csv",
                    "BayesB.csv",
                    "BayesC.csv",
                    "BayesLasso.csv",
                    "BRR.csv")
  
  
  for (i in 1:length(Mod.dat.list)){ # Iterate through files
    # Read in dataframe
    data <- read.csv(Mod.dat.list[i], row.names = 1)
    summary.outputs <- list()
    # Iterate through traits
    for (j in 1:length(BGLR_settings$Traits)){
      # Work on jth trait
      target_string <- BGLR_settings$Traits[j]
      # Find rows for jth trait
      Average.Ref <- which(str_detect(row.names(data), target_string)==T)
      # subset to just data of interest
      TraitoInt <- data[Average.Ref,]
      # Calculate average and standard deviation for single trait
      summary.output <- c(colMeans(as.matrix(TraitoInt)),
                          colSds(as.matrix(TraitoInt)))
      summary.output <- as.data.frame(summary.output)
      row.names(summary.output)<- c("Avg.Cor","Avg.Pval", "Avg.H2", "Avg.Top25",
                                    "Avg.Bot50", "Avg.Time", "Sd.Cor","Sd.Pval",
                                    "Sd.H2", "Sd.Top25", "Sd.Bot50", "Sd.Time")
      colnames(summary.output) <- paste(Mod.dat.list[i], target_string)
      
      if (j < 2) {
        Output.2 <- summary.output
      } else {
        Output.2 <- cbind(Output.2, summary.output)
      }
    }
    
    if (i < 2) {
      Output.3 <- Output.2
    } else {
      Output.3 <- cbind(Output.2, Output.3)
    }
  }

 write.csv(x=Output.3, file="ModSummary.csv")
  
 ModelsofInterest <- c("BRR",
                        "BayesLasso",
                        "BayesC",
                        "BayesB",
                        "BayesA")
  
 Tmp.Sum <- matrix(nrow = 5, ncol=2)
 Tmp.Sum[,1] <- c("BRR",
                   "BL",
                   "BayesC",
                   "BayesB",
                   "BayesA")
  
  for (l in 1:length(ModelsofInterest)){
    target_string <- ModelsofInterest[l]
    Best.Ref <- which(str_detect(colnames(Output.3), target_string)==T)
    tmp <- mean(as.numeric(Output.3[1,Best.Ref]))
    Tmp.Sum[l,2] <- tmp
  }
  ModelTest <- which(Tmp.Sum[,2] == max(Tmp.Sum[,2]))
  BestMod <- Tmp.Sum[1,ModelTest]
  BestMod <<- BestMod
}
  
