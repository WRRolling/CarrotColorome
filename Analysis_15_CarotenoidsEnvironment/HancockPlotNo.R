# Subset to Hancock Data
Hancock.dat <- HPLC.dat[which(HPLC.dat$Location == "Hancock"),]

# Get average of phenotype of plot 
Han.Avg.dat <- Hancock.dat %>%
  group_by(Plot) %>%
  summarise(A450 = mean(avg.ug.fresh450, na.rm=T),
            Alpha = mean(ug_gramfreshAlpha, na.rm=T),
            Beta = mean(ug_gramfreshBeta, na.rm=T),
            Phyto = mean(ug_gramfreshPhyto, na.rm=T),
            Zeta = mean(ug_gramfreshZeta, na.rm=T))
Han.Avg.dat <- Han.Avg.dat[-which(rowSums(is.na(Han.Avg.dat)) > 1),]

# Format Data
for (i in 1:dim(Han.Avg.dat)[1]) {
  ref <- which(Han.Avg.dat$Plot[i]== Hancock.dat$Plot)
  Han.Avg.dat[i,7] <- Hancock.dat$Genotype[ref[1]]
} 
colnames(Han.Avg.dat)[7] <- "Genotype"

# Create Object to Store Output
All.Unique.Geno <- unique(Han.Avg.dat$Genotype)
Summary.Output.Mean <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)
Summary.Output.CI <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)

# Loop through data for Beta 
for (k in 1:length(All.Unique.Geno)){
  Summary.Output.Mean[k,1] <- All.Unique.Geno[k]
  Summary.Output.CI[k,1] <- All.Unique.Geno[k]
  samps <- Han.Avg.dat[which(Han.Avg.dat$Genotype== All.Unique.Geno[k]),]
  samps <- samps %>%
    drop_na(Beta)
    Output <- matrix(nrow=3, ncol=(dim(samps)[1]-1))
    Output[,1] <- c("Avg Mean", "Av LCI" , "Av UCI")
    for (j in 2:(dim(samps)[1]-1)){  
          TmpOut <- matrix(nrow=10,ncol=3)
# Sample x number of plots 10 times
          for (i in 1:10) { 
          samp.1 <- samps[(sample(1:dim(samps)[1],j)),]
          TmpOut[i,1] <- mean(samp.1$Beta, na.rm=T)
          CIVAR <- ci_var(samp.1$Beta)
          TmpOut[i,2] <- CIVAR$interval[1]
          TmpOut[i,3] <- CIVAR$interval[2]
          }
          Output[1,j] <- mean(TmpOut[,1]) 
          Output[2,j] <- mean(TmpOut[,2])
          Output[3,j] <- mean(TmpOut[,3]) 
    }
    for (m in 2:dim(Output)[2]) {
      Summary.Output.Mean[k,m] <- abs(as.numeric(Output[1,m]) - mean(as.numeric(samps$Beta)))
      Summary.Output.CI[k,m] <- as.numeric(Output[3,m]) - as.numeric(Output[2,m])
    }
}



write.csv(x=Summary.Output.CI, file="PlotBetaCI.csv", row.names=F)
write.csv(x=Summary.Output.Mean, file="PlotBetaMean.csv", row.names = F)

# Loop through data for Alpha 
All.Unique.Geno <- unique(Han.Avg.dat$Genotype)
Summary.Output.Mean <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)
Summary.Output.CI <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)


for (k in 1:length(All.Unique.Geno)){
  Summary.Output.Mean[k,1] <- All.Unique.Geno[k]
  Summary.Output.CI[k,1] <- All.Unique.Geno[k]
  samps <- Han.Avg.dat[which(Han.Avg.dat$Genotype== All.Unique.Geno[k]),]
  samps <- samps %>%
    drop_na(Alpha)
  Output <- matrix(nrow=3, ncol=(dim(samps)[1]-1))
  Output[,1] <- c("Avg Mean", "Av LCI" , "Av UCI")
  for (j in 2:(dim(samps)[1]-1)){  
    TmpOut <- matrix(nrow=10,ncol=3)
    # Sample x number of plots 10 times
    for (i in 1:10) { 
      samp.1 <- samps[(sample(1:dim(samps)[1],j)),]
      TmpOut[i,1] <- mean(samp.1$Alpha, na.rm=T)
      CIVAR <- ci_var(samp.1$Alpha)
      TmpOut[i,2] <- CIVAR$interval[1]
      TmpOut[i,3] <- CIVAR$interval[2]
    }
    Output[1,j] <- mean(TmpOut[,1]) 
    Output[2,j] <- mean(TmpOut[,2])
    Output[3,j] <- mean(TmpOut[,3]) 
  }
  for (m in 2:dim(Output)[2]) {
    Summary.Output.Mean[k,m] <- abs(as.numeric(Output[1,m]) - mean(as.numeric(samps$Alpha)))
    Summary.Output.CI[k,m] <- as.numeric(Output[3,m]) - as.numeric(Output[2,m])
  }
}

write.csv(x=Summary.Output.CI, file="PlotAlphaCI.csv", row.names=F)
write.csv(x=Summary.Output.Mean, file="PlotAlphaMean.csv", row.names = F)


# Loop through data for A450 
All.Unique.Geno <- unique(Han.Avg.dat$Genotype)
Summary.Output.Mean <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)
Summary.Output.CI <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)


for (k in 1:length(All.Unique.Geno)){
  Summary.Output.Mean[k,1] <- All.Unique.Geno[k]
  Summary.Output.CI[k,1] <- All.Unique.Geno[k]
  samps <- Han.Avg.dat[which(Han.Avg.dat$Genotype== All.Unique.Geno[k]),]
  samps <- samps %>%
    drop_na(A450)
  Output <- matrix(nrow=3, ncol=(dim(samps)[1]-1))
  Output[,1] <- c("Avg Mean", "Av LCI" , "Av UCI")
  for (j in 2:(dim(samps)[1]-1)){  
    TmpOut <- matrix(nrow=10,ncol=3)
    # Sample x number of plots 10 times
    for (i in 1:10) { 
      samp.1 <- samps[(sample(1:dim(samps)[1],j)),]
      TmpOut[i,1] <- mean(samp.1$A450, na.rm=T)
      CIVAR <- ci_var(samp.1$A450)
      TmpOut[i,2] <- CIVAR$interval[1]
      TmpOut[i,3] <- CIVAR$interval[2]
    }
    Output[1,j] <- mean(TmpOut[,1]) 
    Output[2,j] <- mean(TmpOut[,2])
    Output[3,j] <- mean(TmpOut[,3]) 
  }
  for (m in 2:dim(Output)[2]) {
    Summary.Output.Mean[k,m] <- abs(as.numeric(Output[1,m]) - mean(as.numeric(samps$A450)))
    Summary.Output.CI[k,m] <- as.numeric(Output[3,m]) - as.numeric(Output[2,m])
  }
}

write.csv(x=Summary.Output.CI, file="PlotA450CI.csv", row.names=F)
write.csv(x=Summary.Output.Mean, file="PlotA50Mean.csv", row.names = F)


# Loop through data for Zeta 
All.Unique.Geno <- unique(Han.Avg.dat$Genotype)
Summary.Output.Mean <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)
Summary.Output.CI <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)


for (k in 1:length(All.Unique.Geno)){
  Summary.Output.Mean[k,1] <- All.Unique.Geno[k]
  Summary.Output.CI[k,1] <- All.Unique.Geno[k]
  samps <- Han.Avg.dat[which(Han.Avg.dat$Genotype== All.Unique.Geno[k]),]
  samps <- samps %>%
    drop_na(Zeta)
  Output <- matrix(nrow=3, ncol=(dim(samps)[1]-1))
  Output[,1] <- c("Avg Mean", "Av LCI" , "Av UCI")
  for (j in 2:(dim(samps)[1]-1)){  
    TmpOut <- matrix(nrow=10,ncol=3)
    # Sample x number of plots 10 times
    for (i in 1:10) { 
      samp.1 <- samps[(sample(1:dim(samps)[1],j)),]
      TmpOut[i,1] <- mean(samp.1$Zeta, na.rm=T)
      CIVAR <- ci_var(samp.1$Zeta)
      TmpOut[i,2] <- CIVAR$interval[1]
      TmpOut[i,3] <- CIVAR$interval[2]
    }
    Output[1,j] <- mean(TmpOut[,1]) 
    Output[2,j] <- mean(TmpOut[,2])
    Output[3,j] <- mean(TmpOut[,3]) 
  }
  for (m in 2:dim(Output)[2]) {
    Summary.Output.Mean[k,m] <- abs(as.numeric(Output[1,m]) - mean(as.numeric(samps$Zeta)))
    Summary.Output.CI[k,m] <- as.numeric(Output[3,m]) - as.numeric(Output[2,m])
  }
}

write.csv(x=Summary.Output.CI, file="PlotZetaCI.csv", row.names=F)
write.csv(x=Summary.Output.Mean, file="PlotZetaMean.csv", row.names = F)

# Loop through data for Phyto 

All.Unique.Geno <- unique(Han.Avg.dat$Genotype)
Summary.Output.Mean <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)
Summary.Output.CI <- matrix(nrow=(length(All.Unique.Geno)), ncol=15)

for (k in 1:length(All.Unique.Geno)){
  Summary.Output.Mean[k,1] <- All.Unique.Geno[k]
  Summary.Output.CI[k,1] <- All.Unique.Geno[k]
  samps <- Han.Avg.dat[which(Han.Avg.dat$Genotype== All.Unique.Geno[k]),]
  samps <- samps %>%
    drop_na(Phyto)
  Output <- matrix(nrow=3, ncol=(dim(samps)[1]-1))
  Output[,1] <- c("Avg Mean", "Av LCI" , "Av UCI")
  for (j in 2:(dim(samps)[1]-1)){  
    TmpOut <- matrix(nrow=10,ncol=3)
    # Sample x number of plots 10 times
    for (i in 1:10) { 
      samp.1 <- samps[(sample(1:dim(samps)[1],j)),]
      TmpOut[i,1] <- mean(samp.1$Phyto, na.rm=T)
      CIVAR <- ci_var(samp.1$Phyto)
      TmpOut[i,2] <- CIVAR$interval[1]
      TmpOut[i,3] <- CIVAR$interval[2]
    }
    Output[1,j] <- mean(TmpOut[,1]) 
    Output[2,j] <- mean(TmpOut[,2])
    Output[3,j] <- mean(TmpOut[,3]) 
  }
  for (m in 2:dim(Output)[2]) {
    Summary.Output.Mean[k,m] <- abs(as.numeric(Output[1,m]) - mean(as.numeric(samps$Phyto)))
    Summary.Output.CI[k,m] <- as.numeric(Output[3,m]) - as.numeric(Output[2,m])
  }
}

write.csv(x=Summary.Output.CI, file="PlotPhytoCI.csv", row.names=F)
write.csv(x=Summary.Output.Mean, file="PlotPhytoMean.csv", row.names = F)


