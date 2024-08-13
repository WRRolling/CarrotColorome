subset_and_check_data.Qmat <- function(Pheno, Qmat) {
  # Convert Taxa to character for consistent comparisons
  Pheno$Taxa <- as.character(Pheno$Taxa)
  # Subset and order data based on common Taxa
  Pheno <- Pheno[(Pheno$Taxa %in% Qmat$Taxa),]
  Pheno <- Pheno[order(Pheno$Taxa), ]
  Qmat <- Qmat[(Qmat$Taxa %in% Pheno$Taxa), ]
  Qmat <- Qmat[order(Qmat$Taxa), ]
  Qmat$Taxa <- as.character(Qmat$Taxa)
  
  if (all.equal(Pheno$Taxa, Qmat$Taxa)) {
    print("Qmat and Phenotype Match. Wanna make some noise!")
  } else {
    stop("Abort: Qmat and Phenotype mismatch.")
  }
Pheno <<- Pheno
Qmat <<- Qmat
}
