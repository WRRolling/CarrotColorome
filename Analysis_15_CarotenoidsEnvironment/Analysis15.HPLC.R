######## Part One:Prepare R Environment & Load Phenotypic Data ####### 
# (1.0) Clear Global R 
  rm(list = ls())
# (1.1) Load Libraries
  library('agricolae')
  library('lme4')
  library('lmerTest')
  library('tidyverse')
  library('confintr')
  source("ProduceOutput.R")
# (1.2) Load Data 
  HPLC.dat <- read.csv("Carotenoids.csv", head=T)
# (1.3) 
  
  HPLC.dat <- HPLC.dat[-which(HPLC.dat$ug_gramdryBeta > 5000),] 
  
  
  HPLC.dat <- HPLC.dat[-which(HPLC.dat$Genotype == "Brava"),] 
  HPLC.dat <- HPLC.dat[-which(HPLC.dat$Genotype == "NB3999"),] 
  
# (1.4) Check Data is Laoded OK
  All.Unique.Plots <- unique(HPLC.dat$Plot)
  (dim(HPLC.dat)[1])/length(All.Unique.Plots)
# average 9.34 samples per plot!

######## Part Two: How much variation within plot?####### 
# Question: How many plants we should phenotype per plot?
  source("TestVariance1.R.txt")
  source("TestSampleNumber1.R")
  view(Output)
  source("BetaLoopThrough.R") # 8 Samples
  source("AlphaLoopThrough.R") # 8 Samples
  source("ZetaLoopThrough.R") # 8 Samples
  source("PhytoLoopThrough.R") # 8 Samples
  source("A450LoopThrough.R") # The More the Better. 
######## Part Three: How much within environment?#######
# Question: How many plots do you need? 
  source("HancockPlotNo.R")
  source("ElCentro22PlotNo.R")
  source("ElCentro23PlotNo.R")
 


######## Anova for Trial on a per genotype basis -Beta ######## 
HPLC.dat1408 <- HPLC.dat[which(HPLC.dat$Genotype == "L1408"),]
model.L1408 <- aov(ug_gramfreshBeta ~ Trial, data = HPLC.dat1408)
summary(model.L1408)  
L1408.res <- HSD.test(model.L1408, "Trial")

HPLC.datNb4001<- HPLC.dat[which(HPLC.dat$Genotype == "Nb4001"),]
model.Nb4001<- aov(ug_gramfreshBeta ~ Trial, data = HPLC.datNb4001)
summary(model.Nb4001)  
Nb4001.res <- HSD.test(model.Nb4001, "Trial")

HPLC.datNb8483 <- HPLC.dat[which(HPLC.dat$Genotype == "Nb8483"),]
model.Nb8483<- aov(ug_gramfreshBeta ~ Trial, data = HPLC.datNb8483)
summary(model.Nb8483)  
Nb8483.res <- HSD.test(model.Nb8483, "Trial")


HPLC.datNs5154 <- HPLC.dat[which(HPLC.dat$Genotype == "Ns5154"),]
model.Ns5154 <- aov(ug_gramfreshBeta ~ Trial, data = HPLC.datNs5154)
summary(model.Ns5154)  
Ns5154.res <- HSD.test(model.Ns5154, "Trial")

HPLC.dat.mod <- aov(ug_gramfreshBeta~Genotype, data=HPLC.dat)
summary(HPLC.dat.mod)
HPLC.dat.res <- HSD.test(HPLC.dat.mod, "Genotype")

######## Anova for Trial on a per genotype basis -Alpha ######## 
HPLC.dat1408 <- HPLC.dat[which(HPLC.dat$Genotype == "L1408"),]
model.L1408 <- aov(ug_gramfreshAlpha ~ Trial, data = HPLC.dat1408)
summary(model.L1408)  
L1408.res <- HSD.test(model.L1408, "Trial")

HPLC.datNb4001<- HPLC.dat[which(HPLC.dat$Genotype == "Nb4001"),]
model.Nb4001<- aov(ug_gramfreshAlpha ~ Trial, data = HPLC.datNb4001)
summary(model.Nb4001)  
Nb4001.res <- HSD.test(model.Nb4001, "Trial")

HPLC.datNb8483 <- HPLC.dat[which(HPLC.dat$Genotype == "Nb8483"),]
model.Nb8483<- aov(ug_gramfreshAlpha ~ Trial, data = HPLC.datNb8483)
summary(model.Nb8483)  
Nb8483.res <- HSD.test(model.Nb8483, "Trial")


HPLC.datNs5154 <- HPLC.dat[which(HPLC.dat$Genotype == "Ns5154"),]
model.Ns5154 <- aov(ug_gramfreshAlpha ~ Trial, data = HPLC.datNs5154)
summary(model.Ns5154)  
Ns5154.res <- HSD.test(model.Ns5154, "Trial")

HPLC.dat.mod <- aov(ug_gramfreshAlpha~Genotype, data=HPLC.dat)
summary(HPLC.dat.mod)
HPLC.dat.res <- HSD.test(HPLC.dat.mod, "Genotype")



######## Figure  ######## 

All.dat.Beta <-  ggplot(HPLC.dat, aes(x=Trial, y=ug_gramfreshBeta)) +
  geom_violin() +
  geom_boxplot(width=0.25) +
  ylab("β-carotene") +
  facet_grid(~Genotype)+
  theme_classic()

All.dat.Alpha <-  ggplot(HPLC.dat, aes(x=Trial, y=ug_gramfreshAlpha)) +
  geom_violin() +
  geom_boxplot(width=0.25) +
  ylab("α-carotene") +
  facet_grid(~Genotype)+
  theme_classic()

######## Figure  ######## 
HPLC.dat.CA_22 <- HPLC.dat[which(HPLC.dat$Trial == "CA_22"),]
HPLC.dat.CA_22.Ns5154 <- HPLC.dat.CA_22[which(HPLC.dat.CA_22$Genotype == "Ns5154"),]

Ns5154.Plot <-  ggplot(HPLC.dat.CA_22.Ns5154, aes(x=as.factor(Plot), y=ug_gramfreshBeta)) +
  geom_dotplot(binaxis='y', stackdir='center')+
  ylab("β-carotene") +
  xlab("Plot") +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 22,
    fill = "grey")+
  theme_classic()

# average 
samps <- which(HPLC.dat.CA_22.Ns5154$Plot == "20265")
mean(HPLC.dat.CA_22.Ns5154$ug_gramfreshBeta[samps]) # 46.7 ug/g
samps <- which(HPLC.dat.CA_22.Ns5154$Plot == "20253")
mean(HPLC.dat.CA_22.Ns5154$ug_gramfreshBeta[samps]) # 29.9 

Ns5154.Plot <-  ggplot(HPLC.dat.CA_22.Ns5154, aes(x=as.factor(Plot), y=ug_gramfreshAlpha)) +
  geom_dotplot(binaxis='y', stackdir='center')+
  ylab("α-carotene") +
  xlab("Plot") +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 22,
    fill = "grey")+
  theme_classic()

# average 18
samps <- which(HPLC.dat.CA_22.Ns5154$Plot == "20265")
mean(HPLC.dat.CA_22.Ns5154$ug_gramfreshAlpha[samps]) # 22.3 ug/g
samps <- which(HPLC.dat.CA_22.Ns5154$Plot == "20253")
mean(HPLC.dat.CA_22.Ns5154$ug_gramfreshAlpha[samps]) # 15.3 (68%)

###### average data! ######  
# Trial
Avg.dat <- HPLC.dat %>%
  group_by(Trial) %>%
  summarise(A450 = mean(avg.ug.fresh450, na.rm=T),
            Alpha = mean(ug_gramfreshAlpha, na.rm=T),
            Beta = mean(ug_gramfreshBeta, na.rm=T),
            Phyto = mean(ug_gramfreshPhyto, na.rm=T),
            Zeta = mean(ug_gramfreshZeta, na.rm=T))
# Genotype
Avg.dat <- HPLC.dat %>%
  group_by(Genotype) %>%
  summarise(A450 = mean(avg.ug.fresh450, na.rm=T),
            Alpha = mean(ug_gramfreshAlpha, na.rm=T),
            Beta = mean(ug_gramfreshBeta, na.rm=T),
            Phyto = mean(ug_gramfreshPhyto, na.rm=T),
            Zeta = mean(ug_gramfreshZeta, na.rm=T))



###### ###### 


Avg.dat <- Avg.dat[-which(rowSums(is.na(Avg.dat)) > 1),]

# Format Data
for (i in 1:dim(Avg.dat)[1]) {
  ref <- which(Avg.dat$Plot[i]== HPLC.dat$Plot)
  Avg.dat[i,7] <- HPLC.dat$Genotype[ref[1]]
} 
colnames(Avg.dat)[7] <- "Genotype"

for (i in 1:dim(Avg.dat)[1]) {
  ref <- which(Avg.dat$Plot[i]== HPLC.dat$Plot)
  Avg.dat[i,8] <- HPLC.dat$Location[ref[1]]
} 
colnames(Avg.dat)[8] <- "Location"

for (i in 1:dim(Avg.dat)[1]) {
  ref <- which(Avg.dat$Plot[i]== HPLC.dat$Plot)
  Avg.dat[i,9] <- HPLC.dat$Year[ref[1]]
} 
colnames(Avg.dat)[9] <- "Year"

Avg.dat.1 <- Avg.dat[-which(Avg.dat$Genotype == 'NB3999'),]  

mod1 <- aov(Beta ~ Genotype + Location + Year, data=Avg.dat.1)  

genotype <- HSD.test(mod1, "Genotype")
genotype

Location <- HSD.test(mod1, "Location")
Location  

Avg.dat.2 <- Avg.dat.1[-which(Avg.dat.1$Location == "Hancock"),]  
mod2 <- aov(ug_gramfreshBeta ~ Genotype + Year, data = Avg.dat.2)
summary(mod2)  
year <- HSD.test(mod2, "Year")
year



