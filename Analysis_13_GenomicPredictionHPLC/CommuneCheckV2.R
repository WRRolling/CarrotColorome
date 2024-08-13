Commune.Check <- function (MPL, ETA) {
  if (sum(str_detect(MPL, 'Error')) > 0) {
    Repeats <- which(MPL %like% "Error" == T)
    MPL.1 <- list()
    MPL.1 <- mclapply(Repeats, function (k) {
      PT <- Pheno.WI
      PT[k,i+6] <- NA
      PT[,i+6] <- as.numeric(unlist(PT[,i+6])) 
      fm <- BGLR(y=PT[,i+6], ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
      Pred_pheno <- fm$yHat[[i]]
      MPL.1[[k]]<-Pred_pheno[[1]]
    }, mc.cores = numCores)
    MPL[Repeats] <- MPL.1
  } else {
    print("All Connections Filled")
  }
}
