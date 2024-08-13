 subset_and_check_data <- function(Pheno, Geno.dat, Qmat, Kin) {
  if (!is.data.frame(Geno.dat)) {
    Geno.dat <- as.data.frame(Geno.dat)  # Convert to data frame if needed
  }
  if (!is.data.frame(Kin)) {
    Kin <- as.data.frame(Kin)  # Convert to data frame if needed
  }
  # Ensure consistent column names
  colnames(Geno.dat)[1] <- "Taxa"
  # Convert Taxa to character for consistent comparisons
  Pheno$Taxa <- as.character(Pheno$Taxa)
  # Subset and order data based on common Taxa
  Pheno <- Pheno[Pheno$Taxa %in% Geno.dat$Taxa, ]
  Pheno <- Pheno[order(Pheno$Taxa), ]
  Geno.dat$Taxa <- as.character(Geno.dat$Taxa)
  Geno.dat <- Geno.dat[Geno.dat$Taxa %in% Pheno$Taxa,]
  # Perform checks and provide informative messages
  if (all.equal(Pheno$Taxa, Geno.dat$Taxa)) {
    print("Genotype and Phenotype Match. Gonna let it all hang out!")
  } else {
    stop("Abort: Genotype and Phenotype mismatch.")  # Use stop() for clear error signaling
  }
  Geno.dat <-Geno.dat[,-1]
  Geno.dat <- sapply(Geno.dat, as.numeric)
  Qmat <- Qmat[Qmat$Taxa %in% Pheno$Taxa, ]
  Qmat <- Qmat[order(Qmat$Taxa), ]
  Qmat$Taxa <- as.character(Qmat$Taxa)
  Kin <- Kin[Kin$V1 %in% Pheno$Taxa, ]
  Kin <- Kin[order(Kin$V1), ]
  Kin$V1 <- as.character(Kin$V1)
  
  if (all.equal(Pheno$Taxa, Qmat$Taxa)) {
    print("Qmat and Phenotype Match. Wanna make some noise!")
  } else {
    stop("Abort: Qmat and Phenotype mismatch.")
  }
  
  if (all.equal(Pheno$Taxa, Kin$V1)) {
    print("Kinship and Phenotype Match. Really raise my voice!")
  } else {
    stop("Abort: Kinship and Phenotype mismatch.")
  }
  Kin <- Kin[,-1]
}
