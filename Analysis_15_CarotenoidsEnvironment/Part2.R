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
library(ggpubr)
# (1.2) Load Data 
HPLC.dat <- read.csv("Carotenoids.csv", head=T)
# (1.3) 
hist(HPLC.dat$ug_gramdryBeta)
HPLC.dat <- HPLC.dat[-which(HPLC.dat$ug_gramdryBeta > 5000),] 
hist(HPLC.dat$ug_gramdryBeta)
HPLC.dat <- HPLC.dat[-which(HPLC.dat$Genotype == "Brava"),] 
HPLC.dat <- HPLC.dat[-which(HPLC.dat$Genotype == "NB3999"),] 
# (1.4) Check Data is Laoded OK
All.Unique.Plots <- unique(HPLC.dat$Plot)
(dim(HPLC.dat)[1])/length(All.Unique.Plots)
# average 9.45 samples per plot!

######## Part Two:Significant differences based on Genotype? ####### 
HPLC.dat.mod <- aov(ug_gramdryBeta ~ Genotype + Genotype*Trial, data=HPLC.dat)
summary(HPLC.dat.mod)
HPLC.dat.res <- HSD.test(HPLC.dat.mod, "Genotype")
#groups
#ug_gramdryBeta groups
#L1408        757.6782      a
#Nb8483       530.8291      b
#Nb4001       514.1526      b
#Ns5154       422.2372      c

HPLC.dat.mod.alpha <- aov(ug_gramdryAlpha~Genotype + Genotype * Trial, data=HPLC.dat)
summary(HPLC.dat.mod.alpha)
HPLC.dat.res.alpha <- HSD.test(HPLC.dat.mod, "Genotype")

#ug_gramdryBeta groups
#L1408        757.6782      a
#Nb8483       530.8291      b
#Nb4001       514.1526      b
#Ns5154       422.2372      c


######## Part Three: Significant differences between trials ####### 
HPLC.dat1408 <- HPLC.dat[which(HPLC.dat$Genotype == "L1408"),]
model.L1408 <- aov(ug_gramdryBeta ~ Trial, data = HPLC.dat1408)
summary(model.L1408)  
L1408.res <- HSD.test(model.L1408, "Trial")

HPLC.datNb4001<- HPLC.dat[which(HPLC.dat$Genotype == "Nb4001"),]
model.Nb4001<- aov(ug_gramdryBeta ~ Trial, data = HPLC.datNb4001)
summary(model.Nb4001)  
Nb4001.res <- HSD.test(model.Nb4001, "Trial")

HPLC.datNb8483 <- HPLC.dat[which(HPLC.dat$Genotype == "Nb8483"),]
model.Nb8483<- aov(ug_gramdryBeta ~ Trial, data = HPLC.datNb8483)
summary(model.Nb8483)  
Nb8483.res <- HSD.test(model.Nb8483, "Trial")

HPLC.datNs5154 <- HPLC.dat[which(HPLC.dat$Genotype == "Ns5154"),]
model.Ns5154 <- aov(ug_gramdryBeta ~ Trial, data = HPLC.datNs5154)
summary(model.Ns5154)  
Ns5154.res <- HSD.test(model.Ns5154, "Trial")

######## Anova for Trial on a per genotype basis -Alpha ######## 
HPLC.dat1408 <- HPLC.dat[which(HPLC.dat$Genotype == "L1408"),]
model.L1408 <- aov(ug_gramdryAlpha ~ Trial, data = HPLC.dat1408)
summary(model.L1408)  
L1408.res <- HSD.test(model.L1408, "Trial")

HPLC.datNb4001<- HPLC.dat[which(HPLC.dat$Genotype == "Nb4001"),]
model.Nb4001<- aov(ug_gramdryAlpha ~ Trial, data = HPLC.datNb4001)
summary(model.Nb4001)  
Nb4001.res <- HSD.test(model.Nb4001, "Trial")

HPLC.datNb8483 <- HPLC.dat[which(HPLC.dat$Genotype == "Nb8483"),]
model.Nb8483<- aov(ug_gramdryAlpha ~ Trial, data = HPLC.datNb8483)
summary(model.Nb8483)  
Nb8483.res <- HSD.test(model.Nb8483, "Trial")

HPLC.datNs5154 <- HPLC.dat[which(HPLC.dat$Genotype == "Ns5154"),]
model.Ns5154 <- aov(ug_gramdryAlpha ~ Trial, data = HPLC.datNs5154)
summary(model.Ns5154)  
Ns5154.res <- HSD.test(model.Ns5154, "Trial")

########  Create Figure - Trial Specific Results  ######## 

All.dat.Beta <-  ggplot(HPLC.dat, aes(x=Trial, y=ug_gramdryBeta, fill=Group.Beta)) +
  geom_violin() +
  geom_boxplot(width=0.25) +
  ylab("β-carotene") +
  facet_grid(~Genotype)+
  scale_fill_manual(values=c("grey", "orange", "yellow")) +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") 



All.dat.Alpha <-  ggplot(HPLC.dat, aes(x=Trial, y=ug_gramdryAlpha,fill=Group.Alpha)) +
  geom_violin() +
  geom_boxplot(width=0.25) +
  ylab("α-carotene") +
  facet_grid(~Genotype)+
  scale_fill_manual(values=c("grey", "orange", "yellow")) +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") 



source("ProduceOutput.R")
source("AlphaLoopThrough.R")
Alpha.IntraPlotVar <- as.data.frame(Output)
Alpha.IntraPlotVar2 <- as.data.frame(t(Alpha.IntraPlotVar))
colnames(Alpha.IntraPlotVar2) <- Alpha.IntraPlotVar2[1,]
Alpha.IntraPlotVar2 <- Alpha.IntraPlotVar2[-1,]
Alpha.IntraPlotVar2$Number <- as.numeric(unlist(Alpha.IntraPlotVar2$Number))
Alpha.IntraPlotVar2$`Avg Mean` <- as.numeric(Alpha.IntraPlotVar2$`Avg Mean`)
Alpha.IntraPlotVar2$`Av LCI` <- as.numeric(Alpha.IntraPlotVar2$`Av LCI`)
Alpha.IntraPlotVar2$`Av UCI` <- as.numeric(Alpha.IntraPlotVar2$`Av UCI`)
Alpha.IntraPlotVar2$Difference <- Alpha.IntraPlotVar2$`Av UCI` -Alpha.IntraPlotVar2$`Av LCI`
Alpha.IntraPlotVar2$Percent <- Alpha.IntraPlotVar2$Difference/Alpha.IntraPlotVar2$Difference[1] * 100
Alpha.IntraPlotVar2$Group <- "Alpha"

Alpha.CI <- ggplot(data=Alpha.IntraPlotVar2, aes(x=Number, y=(`Av UCI`-`Av LCI`))) +
  geom_point() +
  xlab(label = "Sample Number") +
  ylab(label = "Confidence Interval") +
  theme_classic()

source("BetaLoopThrough.R") # 8 Samples
Beta.IntraPlotVar <-  as.data.frame(Output)
Beta.IntraPlotVar <- as.data.frame(t(Beta.IntraPlotVar))
colnames(Beta.IntraPlotVar) <- Beta.IntraPlotVar[1,]
Beta.IntraPlotVar2 <- Beta.IntraPlotVar[-c(1,2),]
Beta.IntraPlotVar2$Number <- as.numeric(unlist(Beta.IntraPlotVar2$Number))
Beta.IntraPlotVar2$`Avg Mean` <- as.numeric(Beta.IntraPlotVar2$`Avg Mean`)
Beta.IntraPlotVar2$`Av LCI` <- as.numeric(Beta.IntraPlotVar2$`Av LCI`)
Beta.IntraPlotVar2$`Av UCI` <- as.numeric(Beta.IntraPlotVar2$`Av UCI`) 
Beta.IntraPlotVar2$Difference <- Beta.IntraPlotVar2$`Av UCI` -Beta.IntraPlotVar2$`Av LCI`
Beta.IntraPlotVar2$Percent <- (Beta.IntraPlotVar2$Difference/Beta.IntraPlotVar2$Difference[1]) * 100
Beta.IntraPlotVar2$Group <- "Beta"

Plot.dat <- rbind(Beta.IntraPlotVar2, Alpha.IntraPlotVar2)


ConfIn <- ggplot(data=Plot.dat, aes(x=Number, y=Percent, color=Group)) +
  geom_point(aes(size = 2), show.legend = F) +
  xlab(label = "Sample Number") +
  ylab(label = "Percent") +
  scale_color_manual(values=c("Orange","Grey")) +
  theme(text = element_text(size = 30)) +
  theme(text = element_text(size = 30),
        legend.position = c(.92, .85))

ConfIn + guides(color = guide_legend(override.aes = list(size=5)))

######## Part Two: How much variation within plot?####### 
# Question: How many plants we should phenotype per plot?
source("TestVariance1.R.txt")
source("TestSampleNumber1.R")
view(Output)
source("ZetaLoopThrough.R") # 8 Samples
source("PhytoLoopThrough.R") # 8 Samples
source("A450LoopThrough.R") # The More the Better. 
######## Part Three: How much within environment?#######
# Question: How many plots do you need? 
source("HancockPlotNo.R")
source("ElCentro22PlotNo.R")
source("ElCentro23PlotNo.R")

######## Figure  ######## 
HPLC.dat.CA_22 <- HPLC.dat[which(HPLC.dat$Trial == "CA_22"),]
HPLC.dat.CA_22.L1408 <- HPLC.dat.CA_22[which(HPLC.dat.CA_22$Genotype == "L1408"),]
# average 
samps <- which(HPLC.dat.CA_22.Ns5154$Plot == "20265")
mean(HPLC.dat.CA_22.Ns5154$ug_gramdryBeta[samps]) # 46.7 ug/g
samps <- which(HPLC.dat.CA_22.Ns5154$Plot == "20253")
mean(HPLC.dat.CA_22.Ns5154$ug_gramdryBeta[samps]) # 29.9 

Ns5154.Plot <-  ggplot(HPLC.dat.CA_22.Ns5154, aes(x=as.factor(Plot), y=ug_gramdryBeta)) +
  geom_dotplot(binaxis='y', stackdir='center')+
  ylab("β-carotene") +
  xlab("Plot") +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 22,
    fill = "grey")


  ggtitle(label="Ns5154")+
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 90, hjust = 1))


Ns5154.Plot.Alpha <-  ggplot(HPLC.dat.CA_22.Ns5154, aes(x=as.factor(Plot), y=ug_gramdryAlpha)) +
  geom_dotplot(binaxis='y', stackdir='center')+
  ylab("Alpha") +
  xlab("Plot") +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 22,
    fill = "grey")+
  ggtitle(label="Ns5154")+
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1))



png(file="Figure1.png", width = 800, height=1200)
ggarrange(All.dat.Beta,ConfIn,Ns5154.Plot, 
          nrow=3, 
          ncol = 1, 
          labels = c("A","B","C"))
dev.off()


Avg.dat <- HPLC.dat %>%
  group_by(Genotype) %>%
  summarise(A450 = mean(avg.ug.dry450, na.rm=T),
            Alpha = mean(ug_gramdryAlpha, na.rm=T),
            Beta = mean(ug_gramdryBeta, na.rm=T),
            Phyto = mean(ug_gramdryPhyto, na.rm=T),
            Zeta = mean(ug_gramdryZeta, na.rm=T))


Avg.dat2 <- HPLC.dat %>%
  group_by(Trial) %>%
  summarise(A450 = mean(avg.ug.dry450, na.rm=T),
            Alpha = mean(ug_gramdryAlpha, na.rm=T),
            Beta = mean(ug_gramdryBeta, na.rm=T),
            Phyto = mean(ug_gramdryPhyto, na.rm=T),
            Zeta = mean(ug_gramdryZeta, na.rm=T))

Avg.dat3 <- HPLC.datNs5154 %>%
  group_by(Plot) %>%
  summarise(A450 = mean(avg.ug.dry450, na.rm=T),
            Alpha = mean(ug_gramdryAlpha, na.rm=T),
            Beta = mean(ug_gramdryBeta, na.rm=T),
            Phyto = mean(ug_gramdryPhyto, na.rm=T),
            Zeta = mean(ug_gramdryZeta, na.rm=T))

dat2 <- read.csv("PlotVar.csv", head=F)
colnames(dat2) <- c("inbred", "Concentration", "Plots", "Trait")

Plot.dat <- dat2 %>%
  group_by(Plots,Trait) %>%
  summarise(Average = mean(Concentration, na.rm=T),
            Stdev = sd(Concentration, na.rm=T))

Plot.dat2 <- Plot.dat[-(1:2),]

ConfIn.Plot <- ggplot(data=Plot.dat2, aes(x=Plots, y=Average, color=Trait)) +
  geom_point(aes(size = 2), show.legend = F) +
  geom_errorbar(aes(ymin=Average-Stdev, ymax=Average+Stdev), width=.2)+
  scale_color_manual(values=c("Orange","Grey")) +
  xlab(label = "Plot Number") +
  ylab(label = "Confidence Interval") +
  theme(text = element_text(size = 30),
        legend.position = c(.92, .85))

ConfIn.Plot + guides(color = guide_legend(override.aes = list(size=5)))
