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

# (1.4) Clean R enviornment
rm(Pckg.Lst, package.check)

######## Part Two: Count Number with each color score ######
# (2.1) Count the Color Scores for the exterior
unique(Color.dat$Exterior) # Check which colors are present
length(which(Color.dat$Exterior == "Orange")) # 466
length(which(Color.dat$Exterior == "Yellow")) # 99
length(which(Color.dat$Exterior == "Red")) # 5
length(which(Color.dat$Exterior == "White")) #38
length(which(is.na(Color.dat$Exterior))==T) # 0

# (2.2) Check interior color scores
length(which(Color.dat$Interior == "Orange")) # 405
length(which(Color.dat$Interior== "Yellow")) # 160
length(which(Color.dat$Interior == "Red")) # 0
length(which(Color.dat$Interior == "White")) #43
length(which(is.na(Color.dat$Interior))==T) # 0

# (2.3) Count the categorical color score
length(which(Color.dat$`Color` == "Orange")) # 473
length(which(Color.dat$`Color` == "Yellow")) # 98
length(which(Color.dat$`Color` == "Red")) # 1
length(which(Color.dat$`Color` == "White")) # 36
length(which(is.na(Color.dat$`Color`))==T) # 0


# (2.4) Store count data in an R object
Count.dat <- matrix(nrow=6, ncol=4)
Count.dat[1,] <- c("Score", "Outer.Color", "Inner.Color", "Color")
Count.dat[2,] <- c("Orange", "466", "405", "473")
Count.dat[3,] <- c("Yellow", "99", "160", "98")
Count.dat[4,] <- c("Red", "5", "0", "1")
Count.dat[5,] <- c("White", "38", "43", "36")

write.csv(x=Count.dat, "ColorCounts.csv", row.names = F)




######## Part Three: Identify where HPLC & Color do not match #######
# (3.1) Identify samples with both a color score and HPLC data
Color.dat.sub <- Color.dat[which(Color.dat$`Sample Name` %in% HPLC.dat$Sample.Name),]
HPLC.dat.sub <- HPLC.dat[which(HPLC.dat$Sample.Name %in% Color.dat$`Sample Name`),]

# (3.2) Ensure that orders at the same
Color.dat.sub <- Color.dat.sub[order(Color.dat.sub$`Sample Name`),]
HPLC.dat.sub <- HPLC.dat.sub [order(HPLC.dat.sub $Sample.Name),]

# (3.3) Combine HPLC data and color scores
Pheno.dat <- cbind(Color.dat.sub, HPLC.dat.sub)
Pheno.dat <- Pheno.dat[,-6]

# (3.4) Create Histograms to check color vs HPLC data
ggplot(Pheno.dat, aes(x=Pheno.dat$Color, y=`Lycopene (ug/g)`, fill=Color)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Orange", "Red", "White", "Yellow")) 
  ## Not super informative other than White and Yellow are OK!

# (3.5) Ensure this is correct. 
which(Pheno.dat$Color == "White" & Pheno.dat$`Lycopene (ug/g)` > 0)
which(Pheno.dat$Color == "Yellow" & Pheno.dat$`Lycopene (ug/g)` > 0)

# (3.6) Make these updates
Pheno.dat$ID[11] # "Ames 25040" - Change Exterior Color to Orange
Pheno.dat$ID[206] # "PI 271473" - Change Exterior Color to Orange

# (3.7) Make sure Beta Carotene looks correct. 
ggplot(Pheno.dat, aes(x=Color, y=Pheno.dat$`Beta (ug/g)`, fill=Color)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Orange", "Red", "White", "Yellow")) 

# (3.8) Check which sample is not quite right
which(Pheno.dat$Color == "White" & Pheno.dat$`Beta (ug/g)` > 0)
# 133
Pheno.dat$ID[133] # Leave as is... 

# (3.9) Manual Updates to correct "Yellow" Scores
YellowChk <- which(Pheno.dat$Color == "Yellow" & Pheno.dat$`Beta (ug/g)` > 0)
Pheno.dat.Yellow.Check <- Pheno.dat[YellowChk,]
# Changed ~ 50 scores for yellow to orange,yellow, orange

# (3.10) Reload Color Scores Color Scores. 
Color.dat <- read_csv("ColorScores.csv", col_names=T)

# (3.11) Recreate Phenotype Object
Color.dat.sub <- Color.dat[which(Color.dat$`Sample Name` %in% HPLC.dat$Sample.Name),]
HPLC.dat.sub <- HPLC.dat[which(HPLC.dat$Sample.Name %in% Color.dat$`Sample Name`),]
Color.dat.sub <- Color.dat.sub[order(Color.dat.sub$`Sample Name`),]
HPLC.dat.sub <- HPLC.dat.sub [order(HPLC.dat.sub $Sample.Name),]
Pheno.dat <- cbind(Color.dat.sub, HPLC.dat.sub)
Pheno.dat <- Pheno.dat[,-6]

######## Part Four: Check Orange Exterior and Yellow Interior Samples  #######
# (4.1) Subset to Orange Exterior / Yellow Interior
OrangeChk <- which(Pheno.dat$Exterior == "Orange" & Pheno.dat$Interior == "Yellow")
PhenodatOrange <- Pheno.dat[OrangeChk,]

# (4.2) Check that these samples of alpha/beta
mean(PhenodatOrange$`Beta (ug/g)`, na.rm=T)
mean(PhenodatOrange$`Alpha (ug/g)`, na.rm=T)

# (4.3) Find those extremely low samples
OrangeChk <- which(PhenodatOrange$`Beta (ug/g)` < 10) #47 samples
PhenodatOrange <- PhenodatOrange[OrangeChk,]

# (4.4) Write to file so I can recheck the colors of these samples
write.csv(x=PhenodatOrange, file="OrangeCheckLow.csv") # Remove this HPLC data

# (4.5) Clean R environment
rm(OrangeChk)
rm(PhenodatOrange)

# Updates were needed to the color scores!

# (4.6) Check Color score correctness up to 25 ug/g
OrangeChk <- which(PhenodatOrange$`Beta (ug/g)` > 10 & PhenodatOrange$`Beta (ug/g)` < 25) #47 samples
PhenodatOrange <- PhenodatOrange[OrangeChk,]
write.csv(x=PhenodatOrange, file="OrangeCheckLow2.csv") # Remove this HPLC data
rm(OrangeChk)
rm(PhenodatOrange)
# These Looked OK!

######## Part Five: Check Orange Exterior and Orange Interior Samples  #######
# (5.1) Subset to all orange samples
OrangeChk <- which(Pheno.dat$Exterior == "Orange" & Pheno.dat$Interior == "Orange")
PhenodatOrange <- Pheno.dat[OrangeChk,]
mean(PhenodatOrange$`Beta (ug/g)`, na.rm=T)
hist(PhenodatOrange$`Beta (ug/g)`)

# (5.2) Identify samples with low concentrations of B-carotene
OrangeChk <- which(PhenodatOrange$`Beta (ug/g)` < 25)
PhenodatOrange <- PhenodatOrange[OrangeChk,]
write.csv(x=PhenodatOrange, file="All.OrangeLow.csv") 
# These samples look like the HPLC data was wrong (i.e., color correct)

# (5.3) Check next level 
OrangeChk <- which(PhenodatOrange$`Beta (ug/g)` > 25 & PhenodatOrange$`Beta (ug/g)` < 50)
PhenodatOrange <- PhenodatOrange[OrangeChk,]
write.csv(x=PhenodatOrange, file="All.OrangeLow2.csv") # Remove this HPLC data
# These samples look like the HPLC data was wrong (i.e., color correct)

# (5.4) Check next level 
OrangeChk <- which(PhenodatOrange$`Beta (ug/g)` > 50 & PhenodatOrange$`Beta (ug/g)` < 60)
PhenodatOrange <- PhenodatOrange[OrangeChk,]
write.csv(x=PhenodatOrange, file="All.OrangeLow3.csv") 
# These samples look like the HPLC data was wrong (i.e., color correct)

# (5.5) Check next level 
OrangeChk <- which(PhenodatOrange$`Beta (ug/g)` > 60 & PhenodatOrange$`Beta (ug/g)` < 70)
PhenodatOrange <- PhenodatOrange[OrangeChk,]
write.csv(x=PhenodatOrange, file="All.OrangeLow4.csv") # Remove this HPLC data
# These samples look correct!


######## Part Six: Notes!  #######

# Track Changes
# Remove Color Score and HPLC data for 80186 / PI 223360
# PI 652401 - very orange exterior, but definitely yellow interior. Remove HPLC Data
# NSL 26503 - Remove HPLC Data 
# PI 256065 - Remove HPLC Data
# PI 652205 - HPLC is wrong!
# PI 652255 - HPLC is wrong! 

# Change PI 267091 to Yellow

# Orange Changes in File
# PI 226043 - Orange Throughout
# PI 234619 - Orange Throughout



######## Part Seven: Reload data & Recalculate Counts #######
HPLC.dat <- read_csv("HPLC_dataWI18MU2.csv", col_names = T)
colnames(HPLC.dat)[1] <- "Sample Name"

Color.dat <- read_csv("ColorScoresMU.csv", col_names = T)
colnames(Color.dat)[1] <- "Sample Name"

# 6.3 Calculate a Total of alpha, beta, phytoene, and zeta
HPLC.dat[,8] <- HPLC.dat$`Alpha (ug/g)` + HPLC.dat$`Beta (ug/g)`+
  HPLC.dat$`Phyto (ug/g)`+ HPLC.dat$`Zeta (ug/g)`
colnames(HPLC.dat)[8] <- "Total (ug/g)"

# 6.4 Calculate a ratio of alpha & Beta to Lutein
HPLC.dat[,9] <- (HPLC.dat$`Alpha (ug/g)` + HPLC.dat$`Beta (ug/g)`)/
  HPLC.dat$`Lut (ug/g)`
colnames(HPLC.dat)[9] <- "AB.Ratio"


# Remake matrix to store Color Counts!
Count.dat <- matrix(nrow=6, ncol=4)
Count.dat[1,] <- c("Score", "Outer.Color", "Inner.Color", "Color")
Count.dat[2:6,1] <- c("Orange", "Yellow", "Red", "White", "NaN")

# Count Oranges
Count.dat[2,2] <- length(which(Color.dat$Exterior == "Orange")) 
Count.dat[2,3] <- length(which(Color.dat$Interior == "Orange"))
Count.dat[2,4] <- length(which(Color.dat$`Color` == "Orange")) 

# Count Yellows
Count.dat[3,2] <- length(which(Color.dat$Exterior == "Yellow")) 
Count.dat[3,3] <- length(which(Color.dat$Interior == "Yellow"))
Count.dat[3,4] <- length(which(Color.dat$`Color` == "Yellow")) 

# Count Reds
Count.dat[4,2] <- length(which(Color.dat$Exterior == "Red")) 
Count.dat[4,3] <- length(which(Color.dat$Interior== "Red"))
Count.dat[4,4] <- length(which(Color.dat$`Color` == "Red")) 

# Count White
Count.dat[5,2] <- length(which(Color.dat$Exterior == "White")) 
Count.dat[5,3] <- length(which(Color.dat$Interior == "White"))
Count.dat[5,4] <- length(which(Color.dat$`Color` == "White")) 

# Count Unscored
Count.dat[6,2] <- length(which(is.na(Color.dat$Exterior) == T)) 
Count.dat[6,3] <- length(which(is.na(Color.dat$Interior) == T))
Count.dat[6,4] <- length(which(is.na(Color.dat$`Color`) == T)) 

# Rewrite to file
write.csv(x=Count.dat, "ColorCounts2.csv", row.names = F)
write.csv(x=HPLC.dat, file="UpdatedHPLCdat.csv")


######## Part Eight: HPLC Histograms ######
# Seven Panels of histograms for HPLC data for Supplemental 
# Parts of ggplot
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

