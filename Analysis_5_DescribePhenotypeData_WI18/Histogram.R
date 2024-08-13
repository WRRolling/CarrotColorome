######## Part One:Prepare R Environment & Load Phenotypic Data ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Load libraries
library('tidyverse')
library(rstatix)
library(grid)
library(gridExtra)
library(ggpubr)
# Load phenotype Data

HPLC.dat <- read.csv("../Analysis_4_FigsHPLCdat_WI18/HPLC_dataWI18MU2.csv",
                     head = T)
HPLC.dat[,8] <-(HPLC.dat$Alpha..ug.g. + HPLC.dat$Beta..ug.g.+ HPLC.dat$Phyto..ug.g.
                  + HPLC.dat$Zeta..ug.g.)
colnames(HPLC.dat) <- c("Sample.Name","Alpha (ug/g)", "Beta (ug/g)",
                        "Lut (ug/g)", "Lycopene (ug/g)", "Phyto (ug/g)",
                        "Zeta (ug/g)", "Total (ug/g)")

Alpha <-   ggplot(HPLC.dat, aes(x=`Alpha (ug/g)`)) +
  geom_histogram(binwidth = 50, position="dodge") +
  theme_classic()
# make font much larger!
Alpha <- Alpha + theme(text = element_text(size=36))

# Repeat for other traits!
Beta <- ggplot(HPLC.dat, aes(x=`Beta (ug/g)`)) +
  geom_histogram(binwidth = 100, position="dodge")+
  theme_classic()

Beta <- Beta + theme(text = element_text(size=36))


Lut <-  ggplot(HPLC.dat, aes(x=`Lut (ug/g)`)) +
  geom_histogram(binwidth = 25, position="dodge")+
  theme_classic()

Lut <- Lut + theme(text = element_text(size=36))


Lyco <- ggplot(HPLC.dat, aes(x=`Lycopene (ug/g)`)) +
  geom_histogram(binwidth = 25, position="dodge")+
  theme_classic()

Lyco <- Lyco + theme(text = element_text(size=36))


Phyto <- ggplot(HPLC.dat, aes(x=`Phyto (ug/g)`)) +
  geom_histogram(binwidth = 50, position="dodge")+
  theme_classic()

Phyto <- Phyto + theme(text = element_text(size=36))


Zeta <- ggplot(HPLC.dat, aes(x=`Zeta (ug/g)`)) +
  geom_histogram(binwidth = 25, position="dodge")+
  theme_classic()

Zeta <- Zeta + theme(text = element_text(size=36))


Total <- ggplot(HPLC.dat, aes(x=`Total (ug/g)`)) +
  geom_histogram(binwidth = 500, position="dodge")+
  theme_classic()

Total <- Total + theme(text = element_text(size=36))

# Load packages so I can annotate the figure
library(grid)
library(gridExtra)

# Add (a) to annotate figure
Alpha <- arrangeGrob(Alpha, top = textGrob("(a)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))
Beta <- arrangeGrob(Beta, top = textGrob("(b)", x = unit(0, "npc")
                                         , y   = unit(1, "npc"), just=c("left","top"),
                                         gp=gpar(col="black", fontsize=36, 
                                                 fontfamily="Times Roman")))
Lut <- arrangeGrob(Lut, top = textGrob("(c)", x = unit(0, "npc")
                                       , y   = unit(1, "npc"), just=c("left","top"),
                                       gp=gpar(col="black", fontsize=36, 
                                               fontfamily="Times Roman")))
Lyco <- arrangeGrob(Lyco, top = textGrob("(d)", x = unit(0, "npc")
                                         , y   = unit(1, "npc"), just=c("left","top"),
                                         gp=gpar(col="black", fontsize=36, 
                                                 fontfamily="Times Roman")))
Phyto <- arrangeGrob(Phyto, top = textGrob("(e)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))

Zeta <- arrangeGrob(Zeta, top = textGrob("(f)", x = unit(0, "npc")
                                         , y   = unit(1, "npc"), just=c("left","top"),
                                         gp=gpar(col="black", fontsize=36, 
                                                 fontfamily="Times Roman")))

Total <- arrangeGrob(Total, top = textGrob("(g)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))

# Create Multipane Figure
SupFigHist <- ggarrange(Alpha, 
                        Beta, 
                        Lut, 
                        Lyco, 
                        Phyto, 
                        Zeta, 
                        Total, 
                        ncol=2,
                        nrow=4)

png(filename = "SupHistFig.png", width=1200, height = 1800)
SupFigHist
dev.off()
