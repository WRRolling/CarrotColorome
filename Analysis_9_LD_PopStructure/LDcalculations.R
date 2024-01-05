######## Part One:Prepare R Environment & Load Phenotypic Data #######
# (1.0) Clear Global R
rm(list = ls())
# 1.1 Load libraries
library('tidyverse')

# 1.3 Load Data
dat <- read_table(file="plink.ld",
                 col_names=T)
######## Part Two:Format for LD vs Physical distance #######
# 2.1: Create a dolumn for Distance in kilobases
colnames(dat)[8] <- "Distance (Kb)"
# 2.2: Calculate Physical Distance Between SNPs 
dat$`Distance (Kb)` <- (dat$BP_B - dat$BP_A )/1000

######## Part Three: LD calculations for SNPs < 1kb apart #######
#3.1 less than 10bp apart
LD[1,1] <- "10bp-25bp"
dat.1 <- dat[which(dat$`Distance (Kb)` < .01),]
LD[1,2] <- mean(dat.1$R2) 
LD[1,3] <- "17964 Comparisons"
#3.2 less than 25bp apart
LD[2,1] <- "10bp-25bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .01 & dat$`Distance (Kb)` < .025),]
LD[2,2] <- mean(dat.1$R2) 
LD[2,3] <- "21157 Comparisons"
#3.3 less than 50bp apart
LD[3,1] <- "25bp-50bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .025 & dat$`Distance (Kb)` < .05),]
LD[3,2] <- mean(dat.1$R2) 
LD[3,3] <- "24457 Comparisons"
#3.4 less than 75bp apart
LD[4,1] <- "50bp-75bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .05 & dat$`Distance (Kb)` < .075),]
LD[4,2] <- mean(dat.1$R2) 
LD[4,3] <- "14288  Comparisons"
#3.5 less than 100bp apart
LD[5,1] <- "75bp-100bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .075 & dat$`Distance (Kb)` < .1),]
LD[5,2] <- mean(dat.1$R2) 
LD[5,3] <- "12341  Comparisons"
#3.6 less than 250bp apart
LD[6,1] <- "100bp-250bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .1 & dat$`Distance (Kb)` < .25),]
LD[6,2] <- mean(dat.1$R2) 
LD[6,3] <- length(dat.1$CHR_A)
#3.7 less than 500bp apart
LD[7,1] <- "250bp-500bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .25 & dat$`Distance (Kb)` < .5),]
LD[7,2] <- mean(dat.1$R2) 
LD[7,3] <- length(dat.1$CHR_A)
#3.8 less than 750bp apart
LD[8,1] <- "500bp-750bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .5 & dat$`Distance (Kb)` < .75),]
LD[8,2] <- mean(dat.1$R2) 
LD[8,3] <- length(dat.1$CHR_A)
#3.8 less than 100bp apart
LD[9,1] <- "750bp-1000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > .75 & dat$`Distance (Kb)` < 1),]
LD[9,2] <- mean(dat.1$R2) 
LD[9,3] <- length(dat.1$CHR_A)
######## Part Four: LD calculations for SNPs < 10kb apart #######
LD[10,1] <- "1000bp-1500bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 1 & dat$`Distance (Kb)` < 1.5),]
LD[10,2] <- mean(dat.1$R2) 
LD[10,3] <- length(dat.1$CHR_A)

LD[11,1] <- "1500bp-2000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 1.5 & dat$`Distance (Kb)` < 2),]
LD[11,2] <- mean(dat.1$R2) 
LD[11,3] <- length(dat.1$CHR_A)

LD[12,1] <- "2000-2500bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 2 & dat$`Distance (Kb)` < 2.5),]
LD[12,2] <- mean(dat.1$R2) 
LD[12,3] <- length(dat.1$CHR_A)

LD[13,1] <- "2500-3000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 2.5 & dat$`Distance (Kb)` < 3),]
LD[13,2] <- mean(dat.1$R2) 
LD[13,3] <- length(dat.1$CHR_A)

LD[14,1] <- "3000-4000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 3 & dat$`Distance (Kb)` < 4),]
LD[14,2] <- mean(dat.1$R2) 
LD[14,3] <- length(dat.1$CHR_A)

LD[15,1] <- "4000-5000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 4 & dat$`Distance (Kb)` < 5),]
LD[15,2] <- mean(dat.1$R2) 
LD[15,3] <- length(dat.1$CHR_A)


LD[16,1] <- "5000-6000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 5 & dat$`Distance (Kb)` < 6),]
LD[16,2] <- mean(dat.1$R2) 
LD[16,3] <- length(dat.1$CHR_A)

LD[17,1] <- "6000 - 7000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 6 & dat$`Distance (Kb)` < 7),]
LD[17,2] <- mean(dat.1$R2) 
LD[17,3] <- length(dat.1$CHR_A)

#common practice threshold for LD at 0.2. 
# Threshold Here: 25Kb
493000000/6000
# This would indicate we need 82166 to cover the genome.


######## Part Five: LD calculations for SNPs < 25kb apart #######

LD[18,1] <- "7000-10000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 7 & dat$`Distance (Kb)` < 10),]
LD[18,2] <- mean(dat.1$R2) 
LD[18,3] <- length(dat.1$CHR_A)


LD[19,1] <- "1000-15000bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 10 & dat$`Distance (Kb)` < 15),]
LD[19,2] <- mean(dat.1$R2) 
LD[19,3] <- length(dat.1$CHR_A)


LD[20,1] <- "15000-17500bp"
dat.1 <- dat[which(dat$`Distance (Kb)` > 15 & dat$`Distance (Kb)` < 17.5),]
LD[20,2] <- mean(dat.1$R2) 
LD[20,3] <- length(dat.1$CHR_A)


LD[21,1] <- "17500-20000"
dat.1 <- dat[which(dat$`Distance (Kb)` > 17.5 & dat$`Distance (Kb)` < 20),]
LD[21,2] <- mean(dat.1$R2) 
LD[21,3] <- length(dat.1$CHR_A)


LD[22,1] <- "20000-25000"
dat.1 <- dat[which(dat$`Distance (Kb)` > 20 & dat$`Distance (Kb)` < 25),]
LD[22,2] <- mean(dat.1$R2) 
LD[22,3] <- length(dat.1$CHR_A)
 
######## Part Seven: Large physical distance LD #######
LD[23,1] <- "25000-30000"
dat.1 <- dat[which(dat$`Distance (Kb)` > 25 & dat$`Distance (Kb)` < 30),]
LD[23,2] <- mean(dat.1$R2) 
LD[23,3] <- length(dat.1$CHR_A)

LD[24,1] <- "30000-40000"
dat.1 <- dat[which(dat$`Distance (Kb)` > 30 & dat$`Distance (Kb)` < 40),]
LD[24,2] <- mean(dat.1$R2) 
LD[24,3] <- length(dat.1$CHR_A)

LD[25,1] <- "40000-50000"
dat.1 <- dat[which(dat$`Distance (Kb)` > 40 & dat$`Distance (Kb)` < 50),]
LD[25,2] <- mean(dat.1$R2) 
LD[25,3] <- length(dat.1$CHR_A)

######## Part Eight: Write Output #######

LD.1 <- as.data.frame(LD)
write.csv(x=LD.1, file="LDtable.csv")


quit() :) 
