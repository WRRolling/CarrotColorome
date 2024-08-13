calculate_metrics <- function(Output, Pheno, comparison_column) {
  # Combine predicted values and taxa names
  Predicted.mat <- cbind(Pheno$Taxa, Output)
  
  # Calculate correlation and p-value with chosen column
  Accuracy.dat <- cor.test(Output, Pheno[[comparison_column]])
  posthoc <- Accuracy.dat$estimate^2
  
  # Sort values and data frames based on chosen column
  Predicted.mat.1 <- Predicted.mat[order(Predicted.mat[, 2]), ]
  Pheno.1 <- Pheno[order(Pheno[[comparison_column]]), ]
  
  # Calculate percentile thresholds
  Top25 <- round(length(Pheno.1$Taxa) * 0.75)
  Bottom50 <- round(length(Pheno.1$Taxa) * 0.50)
  
  # Calculate accuracy in top 25% and bottom 50% with chosen column
  Cor25 <- mean(Predicted.mat.1[1:Top25, 1] %in% Pheno.1[1:Top25, "Taxa"])
  Cor50 <- mean(Predicted.mat.1[(Top25 + 1):nrow(Predicted.mat.1), 1] %in% Pheno.1[(Bottom50 + 1):nrow(Pheno.1), "Taxa"])
  
  # Combine results into a vector
  Write.Out <- c(Accuracy.dat$estimate, Accuracy.dat$p.value, posthoc, Cor25, Cor50)
  
  return(Write.Out)
}
