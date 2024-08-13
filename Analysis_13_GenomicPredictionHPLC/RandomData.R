Random.Accuracy <- function() {
  set.seed(1776)
  Random.Output <- lapply(BGLR_settings$Traits, function(trait) {
    PT <- Pheno
    rows <- sample(nrow(PT))
    PT <- PT[rows,]
    random.dat <- cor.test(Pheno[[trait]],PT[[trait]])
    return(random.dat$estimate) 
  })
  Random.Summary <- cbind(BGLR_settings$Traits,
  as.vector(unlist(Random.Output)))  
  write.csv(x=Random.Summary, 
            file="RandomAccuracy.csv",
            row.names = F)
}
