######## Part One:Prepare R Environment & Load Phenotypic Data ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Load packages required for analysis
Pckg.Lst <-c("tidyverse", "ggbeeswarm", "RColorBrewer",
             "ggpubr", "corrplot", "agricolae")
package.check <- lapply(
  Pckg.Lst,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})

# (1.2) Load HPLC Data
HPLC.dat <- read_csv("../Analysis_3_HPLCDat_WI18/HPLC_dataWI18.csv", col_names = T)

# (1.3) Load Color Data
Color.dat <- read_csv("ColorScores.csv", col_names=T)
colnames(Color.dat)[1] <- "Sample Name"

# Clean R enviornment
rm(Pckg.Lst, package.check)

######## Part Two: Count Number with each color score ######
unique(Color.dat$Outer.Color)
length(which(Color.dat$Outer.Color == "Orange")) # 423
length(which(Color.dat$Outer.Color == "Red")) # 4
length(which(is.na(Color.dat$Outer.Color))==T) # 1

unique(Color.dat$Inner.Color)
length(which(Color.dat$Inner.Color == "Orange")) # 364
length(which(Color.dat$Inner.Color == "Yellow")) # 63
length(which(Color.dat$Inner.Color == "White")) # 1
length(which(is.na(Color.dat$Inner.Color))==T) # 0

unique(Color.dat$Color)
length(which(Color.dat$`Color` == "Orange")) # 428

# Manual entered into Supplemental Table


########  Identify Problem samples where HPLC & Figures do not match #######
# Any orange carrots without alpha or Beta?
min(HPLC.dat$`Beta (ug/g)`, na.rm=T)
min(HPLC.dat$`Alpha (ug/g)`, na.rm=T)
# Yes there is a problem there. 
# Manually Filter in Excel/Check pictures for accruate color score





######## Reload data & Recalculate Counts #######
HPLC.dat <- read_csv("HPLC_dataWI18MU.csv", col_names = T)
colnames(Color.dat)[1] <- "Sample Name"

######## HPLC Histograms ######
# Seven Panels of histograms for HPLC data! 

# Create Histograms for Supplemental 

# Parts of ggplto
  # Tell data: Pheno.dat
  # State the Trait (aes(x="Alpha (ug/g)"))
  # Suggest histograms & change bin width 
  # make a simple plot with theme_classic
  
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
# Write Figure to file 
png(filename = "SupHistFig.png", width=1200, height = 1800)
SupFigHist
dev.off()

# Clean R environment
rm(Alpha,
   Beta, 
   Lut, 
   Lyco, 
   Phyto,
   Zeta,
   Total)

# The end :) 

