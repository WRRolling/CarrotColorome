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
HPLC.dat <- read_csv("../Analysis_1_HPLCDat_CA19/HPLC_data.csv", col_names = T)

# (1.3) Load Color Data
Color.dat <- read_csv("Colors.csv", col_names=T)
colnames(Color.dat)[1] <- "Sample Name"

# Clean R enviornment
rm(Pckg.Lst, package.check)

######## Part Two: Count Number with each color score ######
unique(Color.dat$`Outer Color`)
length(which(Color.dat$`Outer Color` == "Orange")) # 505
length(which(Color.dat$`Outer Color` == "Yellow")) # 147
length(which(Color.dat$`Outer Color` == "Red")) # 31
length(which(Color.dat$`Outer Color` == "White")) #25
length(which(is.na(Color.dat$`Outer Color`))==T) # 22


length(which(Color.dat$`Inner Color` == "Orange")) # 386
length(which(Color.dat$`Inner Color` == "Yellow")) # 291
length(which(Color.dat$`Inner Color` == "Red")) # 19
length(which(Color.dat$`Inner Color` == "White")) #34
length(which(is.na(Color.dat$`Inner Color`))==T) # 0

length(which(Color.dat$`Color` == "Orange")) # 508
length(which(Color.dat$`Color` == "Yellow")) # 162
length(which(Color.dat$`Color` == "Red")) # 29
length(which(Color.dat$`Color` == "White")) # 31
length(which(is.na(Color.dat$`Color`))==T) # 0


# Updates Here! 
# Store count data
Count.dat <- matrix(nrow=6, ncol=4)
Count.dat[1,] <- c("Score", "Outer.Color", "Inner.Color", "Color")
Count.dat[2,] <- c("Orange", "505", "386", "508")
Count.dat[3,] <- c("Yellow", "147", "292", "162")
Count.dat[4,] <- c("Red", "31", "19", "29")
Count.dat[5,] <- c("White", "25", "34", "31")
Count.dat[6,] <- c("NA", "22", "0", "0")

write.csv(x=Count.dat, "ColorCounts.csv", row.names = F)

########  Identify Problem samples where HPLC & Figures do not match #######
Color.dat.sub <- Color.dat[which(Color.dat$`Sample Name` %in% HPLC.dat$`Sample Name`),]
HPLC.dat.sub <- HPLC.dat[which(HPLC.dat$`Sample Name` %in% Color.dat$`Sample Name`),]
Pheno.dat <- cbind(Color.dat.sub, HPLC.dat.sub)
Pheno.dat <- Pheno.dat[,-5]

# How should we present the data? 
# Create histograms

ggplot(Pheno.dat, aes(x=Color, y=`Lycopene (ug/g)`, fill=Color)) + 
geom_boxplot() +
scale_fill_manual(values = c("Orange", "Red", "White", "Yellow")) 


ggplot(Pheno.dat, aes(x=Color, y=Pheno.dat$`Beta (ug/g)`, fill=Color)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Orange", "Red", "White", "Yellow")) 

# Samples 90001, 90361, 90512, 90657, 90659 
  # HPLC data does not match picture ... Set All HPLC to NA
  # B-carotene to low for Orange Carrots
# Samples 90674, 90766
  # HPLC data does not match picture ... Set All HPLC to NA
  # B-carotene to high for non orange carrrot
# 90118, 90119, 90165 Color score wrong
  # Outer = NA (Purple), Inner =Yellow, Color =NA
# 90700, Color score wrong
  # NA, Orange, Orange
# 90734
  # Orange, Yellow, Orange
# 90677 
  # HPLC data does not match picture ... Set All HPLC to NA
  # Carrot is all yellow but with high lycopene HPLC. Set to NA
# 90267 and 90672 Color Score Wrong
  # Red, Yellow, Red
# 90373 Color Score Wrong
  # NA, Yellow, NA
# 90010 - Red x3
# 90127 Red Yellow Red
# 90149 Red Orange Orange
# 90675 Red Orange Orange
# 90759 NA, Orange, Yellow

# 90055 90062, 90442, 90481, 90544 - Bad HPLC
  # Red Color low Lycopene. 

# Update Excel Files - did this in excel :p 

######## Reload data & Recalculate Counts #######
HPLC.dat <- read_csv("HPLC_dataMU.csv", col_names = T)
Color.dat <- read_csv("ColorsMU.csv", col_names = T)
colnames(Color.dat)[1] <- "Sample Name"
# Updates Here! 
# Remake matrix to store Color Counts!
Count.dat <- matrix(nrow=6, ncol=4)
Count.dat[1,] <- c("Score", "Outer.Color", "Inner.Color", "Color")
Count.dat[2:6,1] <- c("Orange", "Yellow", "Red", "White", "NaN")

# Count Oranges
Count.dat[2,2] <- length(which(Color.dat$`Outer Color` == "Orange")) 
Count.dat[2,3] <- length(which(Color.dat$`Inner Color` == "Orange"))
Count.dat[2,4] <- length(which(Color.dat$`Color` == "Orange")) 

# Count Yellows
Count.dat[3,2] <- length(which(Color.dat$`Outer Color` == "Yellow")) 
Count.dat[3,3] <- length(which(Color.dat$`Inner Color` == "Yellow"))
Count.dat[3,4] <- length(which(Color.dat$`Color` == "Yellow")) 

# Count Reds
Count.dat[4,2] <- length(which(Color.dat$`Outer Color` == "Red")) 
Count.dat[4,3] <- length(which(Color.dat$`Inner Color` == "Red"))
Count.dat[4,4] <- length(which(Color.dat$`Color` == "Red")) 

# Count White
Count.dat[5,2] <- length(which(Color.dat$`Outer Color` == "White")) 
Count.dat[5,3] <- length(which(Color.dat$`Inner Color` == "White"))
Count.dat[5,4] <- length(which(Color.dat$`Color` == "White")) 

# Count Unscored
Count.dat[6,2] <- length(which(is.na(Color.dat$`Outer Color`) == T)) 
Count.dat[6,3] <- length(which(is.na(Color.dat$`Inner Color`) == T))
Count.dat[6,4] <- length(which(is.na(Color.dat$`Color`) == T)) 

# Rewrite to file
write.csv(x=Count.dat, "ColorCounts.csv", row.names = F)

#Clean R environment #
rm(Count.dat)

######## Recreate Phenotype File ######
Color.dat.sub <- Color.dat[which(Color.dat$`Sample Name` %in% HPLC.dat$`Sample Name`),]
HPLC.dat.sub <- HPLC.dat[which(HPLC.dat$`Sample Name` %in% Color.dat$`Sample Name`),]
Color.dat.sub <- Color.dat.sub[order(Color.dat.sub$`Sample Name`),]
HPLC.dat.sub <-  HPLC.dat.sub[order(HPLC.dat.sub$`Sample Name`),]
Pheno.dat <- cbind(Color.dat.sub, HPLC.dat.sub)
Pheno.dat <- Pheno.dat[,-5]

# Clean R environment
rm(Color.dat,
   Color.dat.sub,
   HPLC.dat,
   HPLC.dat.sub)
######## HPLC Histograms ######
# Seven Panels of histograms for HPLC data! 

# Create Histograms for Supplemental 

# Parts of ggplto
  # Tell data: Pheno.dat
  # State the Trait (aes(x="Alpha (ug/g)"))
  # Suggest histograms & change bin width 
  # make a simple plot with theme_classic
  
Alpha <-   ggplot(Pheno.dat, aes(x=`Alpha (ug/g)`)) +
  geom_histogram(binwidth = 50, position="dodge") +
  theme_classic()
# make font much larger!
Alpha <- Alpha + theme(text = element_text(size=36))

# Repeat for other traits!
Beta <- ggplot(Pheno.dat, aes(x=`Beta (ug/g)`)) +
  geom_histogram(binwidth = 100, position="dodge")+
  theme_classic()

Beta <- Beta + theme(text = element_text(size=36))


Lut <-  ggplot(Pheno.dat, aes(x=`Lut (ug/g)`)) +
  geom_histogram(binwidth = 25, position="dodge")+
  theme_classic()

Lut <- Lut + theme(text = element_text(size=36))


Lyco <- ggplot(Pheno.dat, aes(x=`Lycopene (ug/g)`)) +
  geom_histogram(binwidth = 25, position="dodge")+
  theme_classic()

Lyco <- Lyco + theme(text = element_text(size=36))


Phyto <- ggplot(Pheno.dat, aes(x=`Phyto (ug/g)`)) +
  geom_histogram(binwidth = 50, position="dodge")+
  theme_classic()

Phyto <- Phyto + theme(text = element_text(size=36))


Zeta <- ggplot(Pheno.dat, aes(x=`Zeta (ug/g)`)) +
  geom_histogram(binwidth = 25, position="dodge")+
  theme_classic()

Zeta <- Zeta + theme(text = element_text(size=36))


 Total <- ggplot(Pheno.dat, aes(x=`Total (ug/g)`)) +
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

######## Remove any potential outliers ######

library(rstatix)

# Identify & Remove any outliers

Outliers <- Pheno.dat %>%
  identify_outliers(`Total (ug/g)`)
Pheno.dat.aov <- Pheno.dat[-which(Pheno.dat$`Sample Name` 
                                  %in% Outliers$`Sample Name`),]
Outliers.Alpha <- Pheno.dat %>%
  identify_outliers(`Alpha (ug/g)`)

Pheno.dat.aov <- Pheno.dat.aov[-which(Pheno.dat.aov$`Sample Name` 
                                  %in% Outliers.Alpha$`Sample Name`),]

Outliers.Beta <- Pheno.dat %>%
  identify_outliers(`Beta (ug/g)`)


Pheno.dat.aov <- Pheno.dat.aov[-which(Pheno.dat.aov$`Sample Name` 
                                      %in% Outliers.Beta$`Sample Name`),]

Outliers.Lut <- Pheno.dat %>%
  identify_outliers(`Lut (ug/g)`)

Pheno.dat.aov <- Pheno.dat.aov[-which(Pheno.dat.aov$`Sample Name` 
                                      %in% Outliers.Lut$`Sample Name`),]

Outliers.Lyco <- Pheno.dat %>%
  identify_outliers(`Lycopene (ug/g)`) # Do not remove these because Reds!

Outliers.Phyto <- Pheno.dat %>%
  identify_outliers(`Phyto (ug/g)`) # Do not remove these because Reds!

Pheno.dat.aov <- Pheno.dat.aov[-which(Pheno.dat.aov$`Sample Name` 
                                      %in% Outliers.Phyto$`Sample Name`),]

# Save an HPLC file w/o outliers so I can run a second GWA w/o outliers
write.csv(x=Pheno.dat.aov, file="PhenotypeNoOutliers.csv",
          row.names = F)
# Clean R environment:
rm(Outliers,
   Outliers.Alpha, 
   Outliers.Beta, 
   Outliers.Lut, 
   Outliers.Lyco,
   Outliers.Phyto)

######## Check Normality #######
# Final Filter Remove Those Samples I did not assign a color score too. 
Pheno.dat.aov <- Pheno.dat.aov[-which(is.na(Pheno.dat.aov$Color) == T ),]

# Check Total for normality. 
Mod.Total <- lm(`Total (ug/g)` ~ Color, data=Pheno.dat.aov)
ggqqplot(residuals(Mod.Total)) # Bad 
shapiro_test(residuals(Mod.Total)) # failed :( )
ggqqplot(Pheno.dat.aov, "Total (ug/g)", facet.by="Color")
# No dice!
### Orange & Red are mostly OK. 
# Fits with histogram though, lots of scores near zero for white / yellow


# try a transformation
Pheno.dat.aov[,13] <- log((Pheno.dat.aov$`Total (ug/g)`+1)) 
hist(Pheno.dat.aov$V13)
Mod.Total.Trans <- lm(V13 ~ Color, data=Pheno.dat.aov)
ggqqplot(residuals(Mod.Total.Trans)) # Bad 
# No Dice

# Test Lutein has nearest to normal distribution 
Mod.Lut <- lm(`Lut (ug/g)` ~ Color, data=Pheno.dat.aov)
ggqqplot(residuals(Mod.Lut)) # Decent! 
shapiro_test(residuals(Mod.Lut)) # failed, but less than Total
ggqqplot(Pheno.dat.aov, "Lut (ug/g)", facet.by="Color") # OK
# Fits with Histogram, closest to Gaussian of all traits. 

# OK this data is not normally distributed so it needs non-parametric analysis
# kruskal.test
# pairwise.wilcox. 

# Clean R environment
rm(Mod.Lut, Mod.Total.Trans, Mod.Total)
#remove extra column I created when transforming data
Pheno.dat.aov <- Pheno.dat.aov[,-13]
######## Non-parametric mean separation & Boxplots - alpha #######
# alpha
alpha.np <- kruskal.test(`Alpha (ug/g)`~ Color, data = Pheno.dat.aov)
alpha.np$p.value # pvalue = 2.952569e-78
pairwise.wilcox.test(Pheno.dat.aov$`Alpha (ug/g)`, 
                     Pheno.dat.aov$Color, 
                     p.adjust.method = "BH") 
# Orange significant different from others
# Red, White, and Yellow not significantly different

alpha.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Alpha (ug/g)`, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
  theme_classic()
  alpha.box.1 <- alpha.box + 
                           annotate(geom="text",
                                      x="Orange",
                                      y=125,
                                      label="a") +
                           annotate(geom="text", 
                                       x="Red",
                                       y=125,
                                       label="b") +
                           annotate(geom="text",
                                     x="White",
                                     y=125,
                                     label="b") +
                           annotate(geom="text",
                                       x="Yellow",
                                       y=125,
                                       label="b")


  alpha.box.2 <- alpha.box.1 + theme(text = element_text(size=36)) +
                                theme(legend.position = "none")
  
####### Beta-carotene distribution & mean separation ######  
  beta.np <- kruskal.test(`Beta (ug/g)`~ Color, data = Pheno.dat.aov)
  beta.np$p.value # pvalue = 4.399654e-78
  pairwise.wilcox.test(Pheno.dat.aov$`Beta (ug/g)`, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  beta.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Beta (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  beta.box.1 <- beta.box + 
    annotate(geom="text",
             x="Orange",
             y=1500,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=1500,
             label="b") +
    annotate(geom="text",
             x="White",
             y=1500,
             label="c") +
    annotate(geom="text",
             x="Yellow",
             y=1500,
             label="d")
  
  
  beta.box.2 <- beta.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### Lutein distribution & mean separation #######  
  
  Lut.np <- kruskal.test(`Lut (ug/g)`~ Color, data = Pheno.dat.aov)
  Lut.np$p.value # pvalue = 3.764503e-09
  pairwise.wilcox.test(Pheno.dat.aov$`Lut (ug/g)`, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  lut.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Lut (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  lut.box.1 <- lut.box + 
    annotate(geom="text",
             x="Orange",
             y=150,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=150,
             label="ab") +
    annotate(geom="text",
             x="White",
             y=150,
             label="b") +
    annotate(geom="text",
             x="Yellow",
             y=150,
             label="a")
  
  
  lut.box.2 <- lut.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
  
  
####### Lycopene distribution & mean separation ######  
  Lyco.np <- kruskal.test(`Lycopene (ug/g)`~ Color, data = Pheno.dat.aov)
  Lyco.np$p.value # pvalue = 4.398363e-17
  pairwise.wilcox.test(Pheno.dat.aov$`Lycopene (ug/g)`, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  lyco.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Lycopene (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  lyco.box.1 <- lyco.box + 
    annotate(geom="text",
             x="Orange",
             y=750,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=750,
             label="b") +
    annotate(geom="text",
             x="White",
             y=750,
             label="c") +
    annotate(geom="text",
             x="Yellow",
             y=750,
             label="c")
  
  
  lyco.box.2 <- lyco.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### Phyotoene distribution & mean separation #######  
  
  
  Phyto.np <- kruskal.test(`Phyto (ug/g)`~ Color, data = Pheno.dat.aov)
  Phyto.np$p.value # pvalue = 1.484789e-62
  pairwise.wilcox.test(Pheno.dat.aov$`Phyto (ug/g)`, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  phyto.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Phyto (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  phyto.box.1 <- phyto.box + 
    annotate(geom="text",
             x="Orange",
             y=200,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=200,
             label="a") +
    annotate(geom="text",
             x="White",
             y=200,
             label="b") +
    annotate(geom="text",
             x="Yellow",
             y=200,
             label="b")
  
  
  phyto.box.2 <- phyto.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### Zeta distribution & mean separation #######   
  
  
  zeta.np <- kruskal.test(`Zeta (ug/g)`~ Color, data = Pheno.dat.aov)
  zeta.np$p.value # pvalue = 2.410986e-06
  pairwise.wilcox.test(Pheno.dat.aov$`Zeta (ug/g)`, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  zeta.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Zeta (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  zeta.box.1 <- zeta.box + 
    annotate(geom="text",
             x="Orange",
             y=100,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=100,
             label="ab") +
    annotate(geom="text",
             x="White",
             y=100,
             label="ab") +
    annotate(geom="text",
             x="Yellow",
             y=100,
             label="b")
  
  
  zeta.box.2 <- zeta.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### Total distribution & mean separation ######

  Total.np <- kruskal.test(`Total (ug/g)`~ Color, data = Pheno.dat.aov)
  Total.np$p.value # pvalue = 7.461452e-77
  pairwise.wilcox.test(Pheno.dat.aov$`Total (ug/g)`, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all < 0.05
  # all groups significantly different 
  Total.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=`Total (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  Total.box.1 <- Total.box + 
    annotate(geom="text",
             x="Orange",
             y=2500,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=2500,
             label="b") +
    annotate(geom="text",
             x="White",
             y=2500,
             label="c") +
    annotate(geom="text",
             x="Yellow",
             y=2500,
             label="d")
  
  
  Total.box.2 <- Total.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
####### Ab.ratio distributions & separation #######  
  
  AB.np <- kruskal.test(AB.Ratio ~ Color, data = Pheno.dat.aov)
  AB.np$p.value # pvalue = 4.399654e-78
  pairwise.wilcox.test(Pheno.dat.aov$AB.Ratio, 
                       Pheno.dat.aov$Color, 
                       p.adjust.method = "BH") 
  # all < 0.05
  # all groups significantly different 
  AB.box <- ggplot(data= Pheno.dat.aov, aes(x=Color, y=AB.Ratio, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  AB.box.1 <- AB.box + 
    annotate(geom="text",
             x="Orange",
             y=25,
             label="a") +
    annotate(geom="text", 
             x="Red",
             y=25,
             label="b") +
    annotate(geom="text",
             x="White",
             y=25,
             label="c") +
    annotate(geom="text",
             x="Yellow",
             y=25,
             label="d")
  
  
  AB.box.2 <- AB.box.1 + theme(text = element_text(size=36))   +
    theme(legend.position = "none")
  
####### Combine Figures for manuscript #######  
# Add (a) to annotate figure
  alpha.box.3 <- arrangeGrob(alpha.box.2, top = textGrob("(a)", x = unit(0, "npc")
                                             , y   = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=36, 
                                                     fontfamily="Times Roman")))
  beta.box.3 <- arrangeGrob(beta.box.2 , top = textGrob("(b)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))
  lut.box.3 <- arrangeGrob(lut.box.2, top = textGrob("(c)", x = unit(0, "npc")
                                         , y   = unit(1, "npc"), just=c("left","top"),
                                         gp=gpar(col="black", fontsize=36, 
                                                 fontfamily="Times Roman")))
  lyco.box.3 <- arrangeGrob(lyco.box.2, top = textGrob("(d)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))
  phyto.box.3 <- arrangeGrob(phyto.box.2, top = textGrob("(e)", x = unit(0, "npc")
                                               , y   = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=36, 
                                                     fontfamily="Times Roman")))
  
  zeta.box.3 <- arrangeGrob(zeta.box.2, top = textGrob("(f)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))
  
  Total.box.3 <- arrangeGrob(Total.box.2, top = textGrob("(g)", x = unit(0, "npc")
                                             , y   = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=36, 
                                                     fontfamily="Times Roman")))
  
  AB.box.3 <- arrangeGrob(AB.box.2, top = textGrob("(h)", x = unit(0, "npc")
                                                         , y   = unit(1, "npc"), just=c("left","top"),
                                                         gp=gpar(col="black", fontsize=36, 
                                                                 fontfamily="Times Roman")))
  
  # Create Multipane Figure
  Boxfigs <- ggarrange(alpha.box.3, 
                          beta.box.3, 
                          lut.box.3, 
                          lyco.box.3, 
                          phyto.box.3, 
                          zeta.box.3, 
                          Total.box.3,
                          AB.box.3,
                          ncol=2,
                          nrow=4)
  # Write Figure to file 
  png(filename = "BoxFigsManu.png", width=1200, height = 1800)
  Boxfigs
  dev.off()

####### Correlation Between Traits ##### 
  