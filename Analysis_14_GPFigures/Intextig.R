######## Prepare R environment: Load Data and Library #######
library(tidyverse)
gp.dat <- read.csv("FigureData.csv", head=T)
######## Prepare R environment: CA Alpha data #######

CA.dat <-gp.dat[which(gp.dat$Dataset == "CA"),]
CA.dat.alpha <- CA.dat[which(CA.dat$Trait == "Alpha"),]

Random <- -0.0659
Pop <- -0.024
Core <- 0.165
PopCore <- 0.153
CA.dat.alpha.2 <- CA.dat.alpha[5:9,]

alpha.fig <- ggplot(dat=CA.dat.alpha.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="grey")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.38, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="α-carotene") +
  theme_classic()


alpha.fig <- alpha.fig + theme(text=element_text(size=36),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.position = "none") 

######## Prepare R environment: CA Beta data #######

CA.dat.beta <- CA.dat[which(CA.dat$Trait == "Beta"),]

Random <- CA.dat.beta$Correlation[which(CA.dat.beta$Data == "Random")] 
Pop <- CA.dat.beta$Correlation[which(CA.dat.beta$Data == "Pop")] 
Core <- CA.dat.beta$Correlation[which(CA.dat.beta$Data == "Core")] 
PopCore <- CA.dat.beta$Correlation[which(CA.dat.beta$Data == "Pop|Core")] 

CA.dat.beta.2 <- CA.dat.beta[5:9,]

beta.fig <- ggplot(dat=CA.dat.beta.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="grey")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.28, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="β-carotene") +
  theme_classic()

beta.fig <- beta.fig + theme(text=element_text(size=36),
                               axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.title.y=element_blank(),
                               legend.position = "none") 


######## Prepare R environment: CA lutein data #######

CA.dat.trait.1 <- CA.dat[which(CA.dat$Trait == "Lutein"),]

Random <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Random")] 
Pop <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop")] 
Core <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Core")] 
PopCore <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop|Core")] 

CA.dat.trait.2 <- CA.dat.trait.1[5:9,]

lutein.fig <- ggplot(dat=CA.dat.trait.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="yellow")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.24, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="Lutein") +
  theme_classic()


lutein.fig <- lutein.fig + theme(text=element_text(size=36),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.title.y=element_blank(),
                             legend.position = "none") 

######## Prepare R environment: CA Lycopene  data #######
CA.dat.trait.1 <- CA.dat[which(CA.dat$Trait == "Lycopene"),]

Random <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Random")] 
Pop <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop")] 
Core <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Core")] 
PopCore <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop|Core")] 

CA.dat.trait.2 <- CA.dat.trait.1[5:9,]

Lyco.fig <- ggplot(dat=CA.dat.trait.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="yellow")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.21, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="Lycopene") +
  theme_classic()


Lyco.fig <- Lyco.fig + theme(text=element_text(size=36),
                                 axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.title.y=element_blank(),
                                 legend.position = "none") 

######## Prepare R environment: CA phytoene data #######

CA.dat.trait.1 <- CA.dat[which(CA.dat$Trait == "Phytoene"),]

Random <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Random")] 
Pop <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop")] 
Core <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Core")] 
PopCore <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop|Core")] 

CA.dat.trait.2 <- CA.dat.trait.1[5:9,]

Phyto.fig <- ggplot(dat=CA.dat.trait.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="yellow")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.29, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="Lycopene") +
  theme_classic()


Phyto.fig <- Phyto.fig + theme(text=element_text(size=36),
                             axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.title.y=element_blank(),
                             legend.position = "none") 
######## Prepare R environment: CA Total data #######
CA.dat.trait.1 <- CA.dat[which(CA.dat$Trait == "Total"),]

Random <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Random")] 
Pop <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop")] 
Core <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Core")] 
PopCore <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop|Core")] 

CA.dat.trait.2 <- CA.dat.trait.1[5:9,]

Total.fig <- ggplot(dat=CA.dat.trait.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="yellow")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.33, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="Total") +
  theme_classic()


Total.fig <- Total.fig + theme(text=element_text(size=36),
                               axis.title.y=element_blank(),
                               axis.text.x = element_text(angle = 90),
                               legend.position = "none") 

######## Prepare R environment: CA phytoene data #######
CA.dat.trait.1 <- CA.dat[which(CA.dat$Trait == "Zeta"),]

Random <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Random")] 
Pop <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop")] 
Core <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Core")] 
PopCore <- CA.dat.trait.1$Correlation[which(CA.dat.trait.1$Data == "Pop|Core")] 

CA.dat.trait.2 <- CA.dat.trait.1[5:9,]

Zeta.fig <- ggplot(dat=CA.dat.trait.2, aes(x=Data, y=Correlation, size=2, fatten=2))+
  geom_pointrange(aes(ymin=Correlation-StDev, ymax=Correlation+StDev)) +
  geom_hline(yintercept = Random, linetype="dashed", size=2, color="red") +
  geom_hline(yintercept = Pop, linetype="dashed",size=2, color="orange") +
  geom_hline(yintercept = Core, linetype="dashed", size=2,color="yellow")+ 
  geom_hline(yintercept = PopCore, linetype="dashed",size=2, color="white")+
  geom_hline(yintercept = 0.24, linetype="dashed",size=2, color="black") +
  labs(x="Additions to Model", title="ζ-carotene") +
  theme_classic()


Zeta.fig <- Zeta.fig + theme(text=element_text(size=36),
                               axis.title.y=element_blank(),
                               axis.text.x = element_text(angle = 90),
                               legend.position = "none") 
library(ggpubr)
png(filename="CAFig.png", height = 1500, width=2000)
ggarrange(alpha.fig,beta.fig,lutein.fig,
          Lyco.fig, Phyto.fig, Total.fig, Zeta.fig,
          labels = c("A", "B", "C","D","E","F","G"),
          ncol = 2, nrow = 4)
dev.off()




#################################





Fig.CA <- ggplot(data=CA.dat, aes(x=Trait, y=Correlation, group=Data))+
  geom_point(size=4,aes(shape=Data, color=Trait)) +
  scale_shape_manual(values=c(15,16,17,18,19,13,3,4,10))+
  scale_color_manual(values=c('goldenrod1','darkorange', 'yellow',
                              "tomato1","white","grey","gold1"))+
  ylim(-0.25,0.5)+
  guides(color="none")+
  labs(title="CA_2019_730PI", x="Trait")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))

WI.dat <-gp.dat[which(gp.dat$Dataset == "WI"),]


Fig.WI <- ggplot(data=WI.dat, aes(x=Trait, y=Correlation, group=Data))+
  geom_point(size=4,aes(shape=Data, color=Trait)) +
  scale_shape_manual(values=c(15,16,17,18,19,13,3,4,10))+
  scale_color_manual(values=c('goldenrod1','darkorange', 'yellow',
                              "tomato1","white","grey","gold1"))+
  ylim(-0.25,0.5)+
  guides(color="none")+
  labs(title="WI_2018_605PI", x="Trait")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))



Fig.WI

  library(ggpubr)
png(filename="Intextfig.png", height = 600, width=1800)
ggarrange(Fig.CA,Fig.WI, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()
