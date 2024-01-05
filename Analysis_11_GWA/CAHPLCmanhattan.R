######## Part One: Load Packages ####### 
# (1.1) Load R Packages
library('qqman')
library('tidyverse')
library('calibrate')  

# Part 2 - Create Manhattan Plots for the WI_2018_Trial 
# (2.0) Notes to Create a manhattan plot for alpha
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Alpha.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

# (2.1) Set limit for y axis
y_axis <- 10

# (2.2)  Highlight "known" loci
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
# (2.2) Combine SNP data
SNPsofFun <- c(OrSNPs, YSNPs, Y2SNPs)
# (2.3) Set colors, shapes and sizes for manhattan plot  
fix(manhattan) # alpha = 3

# (2.4) Create manhattan plot
png(file="Alpha.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()


#### 3.0 Beta-carotene concentration 
# (3.1) Format for qqman
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Beta.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

# (3.2) Set colors, shapes and sizes for manhattan plot  
fix(manhattan) # 19 - circle

# (3.3) Create manhattan plot
png(file="Beta.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()


###############
#### 4.0 Lutein concentration 
# (4.1) Format for qqman
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Lutein.csv", head=T)
qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

# (4.2) Set colors, shapes and sizes for manhattan plot  
fix(manhattan) # lutein - 4 - 'x'

# (4.3) Create manhattan plot 
png(file="Lutein.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()

#### 4.0 Lycopene concentration 
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

###############
#### 5.0 Phytoene concentration 
dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Phytoene.csv", head=T)

qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

fix(manhattan) # # Phytoene is Triangle -17

png(file="Phyto.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()

###############
### 6.0 Zeta - 15 - square 
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


###############
### 7.0 Ratio - upside down trialngle - 6 

dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Ratio.csv", head=T)

qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

# (5.6) Set colors, shapes and sizes for manhattan plot  
fix(manhattan)

png(file="Ratio.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()


###############
### 8.0 Total - star - 8 

dat <- read.csv("GAPIT.Association.GWAS_Results.MLM.Total.csv", head=T)

qq.dat <- dat[,1:5]
colnames(qq.dat) <- c("SNP",
                      "CHR",
                      "BP",
                      "P",
                      "zscore")

# (5.6) Set colors, shapes and sizes for manhattan plot  
fix(manhattan)

png(file="Total.png", bg="transparent", width = 1000)
manhattan(qq.dat, main="Carotenoids", ylim=c(0,y_axis), cex = 1.2, cex.lab=2, 
          cex.axis = 2, cex.main=2, col=c("orange","blue"), 
          suggestiveline=F, genomewideline=-log10(7.968381e-07),
          chrlabs=c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6",
                    "Chr7","Chr8","Chr9"),
          annotatePval  = F, highlight = SNPsofFun)
dev.off()



