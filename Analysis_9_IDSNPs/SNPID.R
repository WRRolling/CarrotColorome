######## Part One:Prepare R Environment & Load Phenotypic Data #######
# (1.0) Clear Global R
rm(list = ls())
# 1.1 Load libraries
library('vroom')
library('tidyverse')
library('bigmemory',verbose=F, quietly = T,warn.conflicts=F)
library('biganalytics',verbose=F, quietly = T,warn.conflicts=F)
library('parallel',verbose=F, quietly = T,warn.conflicts=F)
library('MASS',verbose=F, quietly = T, warn.conflicts=F)


# 1.3 Load Data
RSQ.dat <- vroom('../Analysis_8/SNPs.txt', col_names=F)
GBS.dat <- vroom('GBS.SNPPositions.txt', col_names=F)
dim(GBS.dat)

######## Part Two::Format Data #######
# (2.1)Format Chr Name in WI18 data
RSQ.dat$X1 <- gsub('DCARv3_Chr','',RSQ.dat$X1)
RSQ.dat[1,]
# (2.2) Create Columns with concatenated chr and position
RSQ.dat[,5] <- paste(RSQ.dat$X1,"_", RSQ.dat$X2)
GBS.dat[,5] <- paste(GBS.dat$X1, "_", GBS.dat$X2)
colnames(RSQ.dat)[5] <- "X5"
colnames(GBS.dat)[5] <- "X5"
# (2.3) Find number of matching SNPs
length(which(RSQ.dat$X5 %in% GBS.dat$X5))
length(which(GBS.dat$X5 %in% RSQ.dat$X5))
# 86146

######## Part Three::Find Matching SNPs #######
# (3.1) Subset to matching SNPs
RSQ.dat.1 <- RSQ.dat[which(RSQ.dat$X5 %in% GBS.dat$X5),]
GBS.dat.1 <- GBS.dat[which(GBS.dat$X5 %in% RSQ.dat$X5),]
dim(RSQ.dat.1)
dim(GBS.dat.1)
# (3.2) Write to file
write.csv(x=RSQ.dat.1, file="../Analysis_8/SharedSNPs.csv")
write.csv(x=GBS.dat.1, file="GBSSharedSNPs.csv")

######## Part Four: Find Similiar SNPs #######
# (4.1) Create a "Equivalent" Genotype file
length(which(!GBS.dat$X5 %in% RSQ.dat$X5))
# 74056 SNPs
RSQ.dat.2 <- RSQ.dat[which(!RSQ.dat$X5 %in% GBS.dat$X5),]

# (4.2) Identify which SNPs are not in RSQ data
EqiGBS.dat <-  GBS.dat[which(!GBS.dat$X5 %in% RSQ.dat$X5),]
dim(EqiGBS.dat)
numCores <- (detectCores()-1) 
Mark.No <- dim(EqiGBS.dat)[1]
Equi.Mark <- list()
# (4.3) Get Equivalent Markers
Equi.Mark <- mclapply(1:Mark.No, function (i) {
  Tmp.dat <- RSQ.dat.2[which(RSQ.dat.2$X1 == EqiGBS.dat$X1[i]),]
  OK <- which.min(abs(Tmp.dat$X2 - EqiGBS.dat$X2[i]))
  Equi.Mark[i] <- Tmp.dat$X5[OK]}, mc.cores = numCores)

# (4.4) Identify and Remove duplicated SNP
Remove <- which(duplicated(Equi.Mark)==T)
Equi.Mark.1 <- as.vector(unlist(Equi.Mark))

# (4.5) Write output for Rseq
Equi.Mark.2 <- Equi.Mark.1[-Remove]
write.csv(x=Equi.Mark.1, file="../Analysis_8/Equil.csv")

# (4.6) Write output for GBS
EqiGBS.dat.1 <- EqiGBS.dat[-Remove,]
write.csv(x=EqiGBS.dat.1 , file="GBSEquil.csv")

quit()
#The End :)
