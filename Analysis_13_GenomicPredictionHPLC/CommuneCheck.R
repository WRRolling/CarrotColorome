Commune.Check <- function (MPL, ETA) {
    if (sum(str_detect(MPL, 'Error')) > 0) {
      Repeats <- which(MPL %like% "Error" == T)
      MPL.1 <- list()
      MPL.1 <- mclapply(Repeats, function (k) {
        PT <- Pheno
        PT[[trait]][k] <- NA
        PT[[trait]] <- as.numeric(unlist(PT[[trait]]))
        fm <- BGLR(y=PT[[trait]], ETA=ETA, nIter=BGLR_settings$nIter, burnIn=BGLR_settings$burnIn, verbose=FALSE)
        Pred_pheno <- fm$yHat[[i]]
        MPL.1[[k]]<-Pred_pheno[[1]]
      }, mc.cores = BGLR_settings$numCores)
      MPL[Repeats] <- MPL.1
    } else {
      print("All Connections Filled")
      }
}

