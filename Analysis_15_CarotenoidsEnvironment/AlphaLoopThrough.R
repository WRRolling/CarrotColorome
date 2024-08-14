Summary.Output.Mean <- matrix(nrow=(length(All.Unique.Plots)), ncol=10)
Summary.Output.CI <- matrix(nrow=(length(All.Unique.Plots)), ncol=10)

for (k in 1:length(All.Unique.Plots)){
  Summary.Output.Mean[k,1] <- All.Unique.Plots[k]
  Summary.Output.CI[k,1] <- All.Unique.Plots[k]
  samps <- HPLC.dat[which(HPLC.dat$Plot == All.Unique.Plots[k]),]
  samps <- samps %>%
    drop_na(ug_gramdryAlpha)
  if(dim(samps)[1] < 2){
    print("Forgit abow dit")
  } else {
    Output <- matrix(nrow=4, ncol=10)
    Output[,1] <- c("Number", "Avg Mean", "Av LCI" , "Av UCI")
    Output[1,2:10] <- 2:10
    for (j in 2:dim(samps)[1]) {  
      TmpOut <- matrix(nrow=10,ncol=3)
      for (i in 1:10) {
        samp.1 <- samps[(sample(1:dim(samps)[1], j)),]
        TmpOut[i,1] <- mean(samp.1$ug_gramdryAlpha, na.rm=T)
        CIVAR <- ci_var(samp.1$ug_gramdryAlpha)
        TmpOut[i,2] <- CIVAR$interval[1]
        TmpOut[i,3] <- CIVAR$interval[2]
      }
      Output[2,j] <- mean(TmpOut[,1]) 
      Output[3,j] <- mean(TmpOut[,2])
      Output[4,j] <- mean(TmpOut[,3]) 
    }
    Output[2,10] <- mean(samps$ug_gramdryAlpha)
    CIVAR <- ci_var(samps$ug_gramdryAlpha)
    Output[3,10] <-CIVAR$interval[1]
    Output[4,10] <- CIVAR$interval[2]
    ProduceOutput() 
  }
}

remove <- which(rowSums(is.na(Summary.Output.Mean)) > 5)
Summary.Output.CI <- Summary.Output.CI[-remove,]
Summary.Output.Mean <- Summary.Output.Mean[-remove,]

Big.Diff <- sd(Summary.Output.CI[,10], na.rm=T) * 3
Mean <- mean(Summary.Output.CI[,10], na.rm=T) 
Threshold <- Big.Diff + Mean
remove <- which(Summary.Output.CI[,10] > Threshold)
Summary.Output.CI <- Summary.Output.CI[-remove,]
Summary.Output.Mean <- Summary.Output.Mean[-remove,]

write.csv(x=Summary.Output.CI, file="AlphaCI.csv", row.names=F)
write.csv(x=Summary.Output.Mean, file="AlphaMean.csv", row.names = F)

