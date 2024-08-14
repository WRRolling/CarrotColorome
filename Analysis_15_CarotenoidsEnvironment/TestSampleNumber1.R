Output <- matrix(nrow=4, ncol=10)
Output[,1] <- c("Number", "Avg Mean", "Av LCI" , "Av UCI")
Output[1,2:10] <- 2:10

for (j in 2:9) {  
TmpOut <- matrix(nrow=10,ncol=3)

for (i in 1:10) {
  samp.1 <- samps[(sample(1:dim(samps)[1], j)),]
  TmpOut[i,1] <- mean(samp.1$ug_gramfreshBeta)
  CIVAR <- ci_var(samp.1$ug_gramfreshBeta)
  TmpOut[i,2] <- CIVAR$interval[1]
  TmpOut[i,3] <- CIVAR$interval[2]
}

  Output[2,j] <- mean(TmpOut[,1]) 
  Output[3,j] <- mean(TmpOut[,2])
  Output[4,j] <- mean(TmpOut[,3]) 
} 

Output[2,10] <- mean(samps$ug_gramfreshBeta)
CIVAR <- ci_var(samps$ug_gramfreshBeta)
Output[3,10] <-CIVAR$interval[1]
Output[4,10] <- CIVAR$interval[2]