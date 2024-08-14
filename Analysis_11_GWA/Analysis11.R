######## Part One:Prepare R Environment & Format Genotypic Data ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Load GAPIT functions/packages
  # Saved Functions because version updates result in small changes to results
source("InputData/GAPITPackages.R.txt")
source("InputData/GAPITFunctions.R.txt")
library('tidyverse')
library('vroom')
# (1.2) Load Phenotypic Data
myY <- read.csv("InputData/ColorGWA.Pheno.csv", head=T)
myY$Exterior <- as.numeric(myY$Exterior)
myY$Interior <- as.numeric(myY$Interior)
myY$Color <- as.numeric(myY$Color)
# (1.3) Load Marker Data 
myGM <- read.csv("InputData/ContinousGAPITGM.csv", head=T)
myGM$Chromosome <- as.numeric(myGM$Chromosome)
# (1.4) Load Genotype Data
myGD <- vroom("InputData/Probgeno.csv")
colnames(myGD)[1] <- "Taxa"
# (1.5) Load CV
myCV <- read.csv("InputData/Qmat.csv", head=T)

######## Notes ####### 

# I tried Red vs Orange
# I tried white vs yellow
# I tried white + yellow vs Orange
# I did orange vs white but got redundant information as orange v yellow
# The three analyses I am presenting are the best!
# I set the significance on the full genotype file. 
  # This is more stringent then MAF filtered data

######## Part Two: Yellow vs Orange #######
# (2.1) Create director for analyses
dir.create("../YellowvOrange/")
setwd("../YellowvOrange/")
# (2.2) Subset to Just Yellow and Orange Accessions
Colors <- c(1,2)
Keepers <- which(myY$Color %in% Colors)
# (2.3) Subset Phenotype
myY.1 <- myY[Keepers,]
# (2.4) Subset Genotype
myGD.1 <- myGD[which(myGD$Taxa %in% myY.1$Taxa),]
# (2.5) Subset CV
myCV.1 <- myCV[which(myCV$Taxa %in% myY.1$Taxa),]
# (2.6) Run GWAS
myGAPIT=GAPIT(
  Y=myY.1[,c(1,5)], #fist column is ID
  GD=myGD.1,
  GM=myGM,
  PCA.total=2,
  SNP.MAF=0.05, 
  model="MLM")  
######## Part Three: Significance Threshold #######
# (3.1) Format Input data
dat.1 <- myGD
dat.2 <- as.data.frame(t(dat.1))
dat.3 <- mutate_if(dat.2, is.factor, ~ as.numeric(levels(.x))[.x])
dat.4 <- data.matrix(dat.3)
dat.4.5 <- dat.4[-1,]
dat.5 <- scale(dat.4.5)
# (3.2) Load simpleM Functions
PCA_cutoff <- 0.995
Meff_PCA <- function(eigenValues, percentCut){
  totalEigenValues <- sum(eigenValues)
  myCut <- percentCut*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if(myEigenSum <= myCut){
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    }
    else{
      break
    }
  }
  return(index_Eigen)
}
#3.3 infer the cutoff => Meff
inferCutoff <- function(dt_My){
  CLD <- cor(dt_My)
  eigen_My <- eigen(CLD)
  
  # PCA approach
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}
numLoci <- length(dat.5[, 1])
simpleMeff <- NULL
fixLength <- 133
i <- 1
myStart <- 1
myStop <- 1
while(myStop < numLoci){
  myDiff <- numLoci - myStop
  if(myDiff <= fixLength) break
  
  myStop <- myStart + i*fixLength - 1
  snpInBlk <- t(dat.5[myStart:myStop, ])
  MeffBlk <- inferCutoff(snpInBlk)
  simpleMeff <- c(simpleMeff, MeffBlk)
  myStart <- myStop+1
}
# (3.4) Run simpleM
snpInBlk <- t(dat.5[myStart:numLoci, ])
MeffBlk <- inferCutoff(snpInBlk)
simpleMeff <- c(simpleMeff, MeffBlk)
print(simpleMeff)
# (3.5) Record output is here:
Tot.SNP <- cat("Total number of SNPs is: ", numLoci, "\n")
Ind.Test <- cat("Inferred Meff is: ", sum(simpleMeff), "\n")
# 7.968381e-07
# Sig Threshold 6.1
######## Part Four: Manhattan Plot ####### 
# (4.1) Load R Packages
library('qqman')
library('tidyverse')
library('calibrate')  
# (4.2) Format for qqman
dat <- read.csv("WhitevYellow/GAPIT.Association.GWAS_Results.MLM.Color.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

# (4.3) Identify Highest association for y axis 
max(-log10(qq.dat$P)) # 7.773654!
y_axis <- 10

# (4.4)  Highlight "known" loci
Chr03 <- which(qq.dat$CHR == 3)
Chr03 <- qq.dat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5012152 & Chr03$BP < 5133744),]
Chr05 <- which(qq.dat$CHR == 5)
Chr05 <- qq.dat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29982367 & Chr05$BP < 30043214),]   
Chr07 <- which(qq.dat$CHR == 7)
Chr07 <- qq.dat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39189059),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
# (5.5) Combine SNP data
SNPsofFun <- c(OrSNPs, YSNPs, Y2SNPs)
# (5.6) Set colors, shapes and sizes for manhattan plot  
fix(manhattan)

# Transparent background to overlay plots

png(file="WhitevYellow.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Yellow vs White", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("yellow3","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()

# (5.2) Format for qqman
dat <- read.csv("YellowvOrange/GAPIT.Association.GWAS_Results.MLM.Color.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan)

max(-log10(qq.dat$P)) # 7.773654!
y_axis <- 15

# (4.4)  Highlight "known" loci
Chr03 <- which(qq.dat$CHR == 3)
Chr03 <- qq.dat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5012152 & Chr03$BP < 5133744),]
Chr05 <- which(qq.dat$CHR == 5)
Chr05 <- qq.dat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29982367 & Chr05$BP < 30043214),]   
Chr07 <- which(qq.dat$CHR == 7)
Chr07 <- qq.dat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39189059),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
# (4.5) Combine SNP data
SNPsofFun <- c(OrSNPs, YSNPs, Y2SNPs)


png(file="YellowvOrange.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Yellow vs Orange", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()

######## Part Five: Extract Markers for Haploview ####### 
rm(list = ls())
library('vroom')
myGD <- vroom("InputData/Probgeno.csv")
colnames(myGD)[1] <- "Taxa"

# Chr03 PDS
Front <- which(colnames(myGD)=="S3_13947216_T")
End <- which(colnames(myGD)=="S3_14144996_A")
Sub.Front <- Front -50
Sub.End <- End + 50
PDS <- myGD[,Sub.Front:Sub.End]
write.csv(x=PDS, file="../Analysis_12Haploview/PDS.csv", row.names = F)

# Chr03 OR
Or <- which(colnames(myGD)=="S3_5240775_G")
Sub.Front <- Or -50
Sub.End <- Or + 50
ORRegion <- myGD[,Sub.Front:Sub.End]
write.csv(x=ORRegion, file="../Analysis_12Haploview/OR.csv", row.names = F)


# Chr07 
Seven01 <- which(colnames(myGD)=="S7_3677561_G")
Sub.Front <- Seven01 -50
Sub.End <- Seven01 + 50
SevenRegion <- myGD[,Sub.Front:Sub.End]
write.csv(x=SevenRegion, file="../Analysis_12Haploview/Seven.csv", row.names = F)


# Y2
Y2.2 <- which(colnames(myGD)=="S7_39185657_T")
Y2.1 <- which(colnames(myGD)=="S7_39012349_GCCTGCTG")
Sub.Front <- Y2.1 -50
Sub.End <- Y2.2 + 50
Y2 <- myGD[,Sub.Front:Sub.End]
write.csv(x=Y2, file="../Analysis_12Haploview/Y2.csv", row.names = F)

######## Part Six: Core Color #########
# (6.1) Create directory for output
dir.create("CoreColor/")
setwd("CoreColor/")
# (6.2) ID White and Yellows for acurate MAF filtering
Colors <- c(1)
Keepers <- which(myY$Color %in% Colors)
# (6.3) Subset Phenotype
myY.1 <- myY[Keepers,]
myY.1 <- myY.1[-which(myY.1$Interior == 4),]
# (6.4) Subset Genotype
myGD.1 <- myGD[which(myGD$Taxa %in% myY.1$Taxa),]
# (6.5) Subset CV
myCV.1 <- myCV[which(myCV$Taxa %in% myY.1$Taxa),]
# (6.6) Run GWAS
myGAPIT=GAPIT(
  Y=myY.1[,c(1,3)], #fist column is ID
  GD=myGD.1,
  GM=myGM,
  #CV=myCV.1,
  SNP.MAF=0.15, 
  model="MLM")  

######## Part Seven: Wisconsin - Oranges #########
# (7.1) Create directory for output
dir.create("../Wisco18/")
setwd("../Wisco18/")
myY.2 <- myY.1[1:402,]
# (7.2) Run Anaysis
myGAPIT=GAPIT(
  Y=myY.2[,c(1,6:12)], #fist column is ID
  GD=myGD.1,
  GM=myGM,
  CV=myCV.1,
  SNP.MAF=0.15, 
  model="MLM")  
# (7.3) Plot Results!
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Alpha.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
y_axis <- 10
png(file="Alpha.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Beta.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Beta.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Lutein.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Lutein.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Lycopene.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) # Lycopene is diamond - 18
png(file="Lyco.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Phytoene.csv", head=T)

qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Phyto.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Zeta.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Zeta.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Total.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan)

png(file="Total.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()

######## Part Eight: CA - Oranges #########
# (8.1) Create Directory For Output
dir.create("../CA19/")
setwd("../CA19/")
myY.2 <- myY.1[403:658,]
# (8.2) Run Analysis
myGAPIT=GAPIT(
  Y=myY.2[,c(1,6:12)], #fist column is ID
  GD=myGD.1,
  GM=myGM,
  CV=myCV.1,
  SNP.MAF=0.15, 
  model="MLM")  
# (8.3) Plot Results!
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Alpha.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
y_axis <- 10
png(file="Alpha.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Beta.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Beta.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Lutein.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Lutein.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Lycopene.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) # Lycopene is diamond - 18
png(file="Lyco.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Phytoene.csv", head=T)

qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Phyto.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Zeta.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
fix(manhattan) 
png(file="Zeta.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Total.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

png(file="Total.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()


######## Part Nine: Combined Results #########
# (8.1) Create Directory For Output
dir.create("../CombinedHPLC/")
setwd("../CombinedHPLC/")
# (8.2) Run Analysis
myGAPIT=GAPIT(
  Y=myY.1[,c(1,6:12)], #fist column is ID
  GD=myGD.1,
  GM=myGM,
  CV=myCV.1,
  SNP.MAF=0.15, 
  model="MLM")  
# (8.3) Plot Results!
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Alpha.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
y_axis <- 10
png(file="Alpha.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
fix(manhattan)
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Total.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")
png(file="Total.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)