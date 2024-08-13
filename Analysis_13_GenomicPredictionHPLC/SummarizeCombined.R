CombinedModsummary <- function() {
  Mod.dat.list <- c("ModCore.csv",
                    "ModPop.csv",
                    "ModKin.csv",
                    "ModKinPop.csv",
                    "ModKinPopCore.csv")
  
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
  
  write.csv(x=Output.3, file="ModCombinedSummary.csv")
}
