######## Part One: Prepare R Environment & Load Phenotypic Data ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Load packages required for analysis
library(grid)
library(gridExtra)
library(tidyverse)
library(rstatix)
library(ggpubr)

# (1.2) Load HPLC Data
HPLC.dat <- read_csv("../Analysis_2_FigsHPLCdat_CA19/HPLC_dataMU.csv", col_names = T)

# (1.3) Load Color Data
Color.dat <- read_csv("../Analysis_2_FigsHPLCdat_CA19/ColorsMU.csv", col_names=T)
colnames(Color.dat)[1]<- "Sample Name"

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

######## Part Two:Calculate Ratio of alpha + Beta relative to lutein
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Orange")], na.rm=T) # 15.3
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Yellow")], na.rm=T) # 0.1
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Red")], na.rm=T) # 2.2
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "White")], na.rm=T) # 0
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Orange")], na.rm=T) # 7.9
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Yellow")], na.rm=T) # 0.4
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Red")], na.rm=T) # 2.0
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "White")], na.rm=T) # 0

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

######## Part Two: Identify Orange Outliers #######
# (2.1) Subset to orange accessions
Orange.dat <- Pheno.dat[Pheno.dat$Color == "Orange",]
# (2.2) Identify Alpha Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Alpha (ug/g)`)
Outliers.Alpha <- Outliers[which(Outliers$is.extreme == T),]
# (2.3) Identify Beta Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Beta (ug/g)`)
Outliers.Beta <- Outliers[which(Outliers$is.extreme == T),]
# (2.4) Identify Lutein Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Lut (ug/g)`)
Outliers.Lut <- Outliers[which(Outliers$is.extreme == T),] # None
# (2.5) Identify Lycopene Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Lycopene (ug/g)`)
Outliers.Lyco <- Outliers[which(Outliers$is.extreme == T),]
# (2.6) Identify Phytoene Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Phyto (ug/g)`)
Outliers.Phyto <- Outliers[which(Outliers$is.extreme == T),]
# (2.7) Identify Phytoene Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Zeta (ug/g)`)
Outliers.Zeta <- Outliers[which(Outliers$is.extreme == T),] 
  # Don't remove "84 samples (20%) were "outliers"
    # this doesn't seem like outliers then... 
# (2.8) Identify Phytoene Outliers
Outliers <- Orange.dat %>%
  identify_outliers(`Total (ug/g)`)
Outliers.Tot <- Outliers[which(Outliers$is.extreme == T),]

# (2.9) Create Object to fitler samples with!
FilteredSamples <- rbind(Outliers.Alpha, 
                         Outliers.Beta,
                         Outliers.Lyco,
                         Outliers.Phyto,
                         Outliers.Tot)
FilteredSamples <- FilteredSamples[-c(5,6,8,10,15:23),]



# (2.10) Clean R environment:
rm(Outliers,
   Outliers.Alpha, 
   Outliers.Beta, 
   Outliers.Lut, 
   Outliers.Lyco,
   Outliers.Phyto)

# (2.11) Filtered Orange Data
Orange.dat.2 <- Orange.dat[-which(Orange.dat$`Sample Name` %in% 
                                    FilteredSamples$`Sample Name`),]

# Write these samples to a file for subsequent analyses
write.csv(x=Orange.dat.2,
            file="FilteredOrangeDat.csv",
          row.names=F)

# (2.11) Clean R environment:
  rm(Outliers.Zeta,
     Outliers.Tot,
     Orange.dat,
     Orange.dat.2,
     HPLC.dat)

  
  
  
  
######## Part Three: Test for Normality #######
# (3.1) Create object with filtered data
Pheno.dat.Filt <- Pheno.dat[-which(Pheno.dat$`Sample Name`
                                   %in% FilteredSamples$`Sample Name`),]  
Pheno.dat.Filt <- Pheno.dat.Filt[-which(is.na(Pheno.dat.Filt$Color)==T),]
# (3.2) Check Total for normality. 
Mod.Total <- lm(`Total (ug/g)` ~ Color, data=Pheno.dat.Filt)
ggqqplot(residuals(Mod.Total)) # Bad 
shapiro_test(residuals(Mod.Total)) # failed :( )
ggqqplot(Pheno.dat.Filt, "Total (ug/g)", facet.by="Color")
# No dice!
### Orange & Red are mostly OK. 
# Fits with histogram though, lots of scores near zero for white / yellow
# (3.3) try a transformation
Pheno.dat.Filt[,13] <- log((Pheno.dat.Filt$`Total (ug/g)`+1)) 
hist(Pheno.dat.Filt$V13)
Mod.Total.Trans <- lm(V13 ~ Color, data=Pheno.dat.Filt)
ggqqplot(residuals(Mod.Total.Trans)) # Bad 
# No Dice

# (3.4) Test Lutein has nearest to normal distribution 
Mod.Lut <- lm(`Lut (ug/g)` ~ Color, data=Pheno.dat.Filt)
ggqqplot(residuals(Mod.Lut)) # Decent! 
shapiro_test(residuals(Mod.Lut)) # failed, but less than Total
ggqqplot(Pheno.dat.Filt, "Lut (ug/g)", facet.by="Color") # OK
# Fits with Histogram, closest to Gaussian of all traits. 

# (3.5) Clean R environment
rm(Mod.Lut, Mod.Total.Trans, Mod.Total)
#remove extra column I created when transforming data
Pheno.dat.Filt <- Pheno.dat.Filt[,-13]


######## Part Four: Nonparametric comparison of means #######
# OK this data is not normally distributed so it needs non-parametric analysis
# kruskal.test
# pairwise.wilcox. 

######## (4.1) Non-parametric mean separation & Boxplots - alpha #######
# alpha
alpha.np <- kruskal.test(`Alpha (ug/g)`~ Color, data = Pheno.dat.Filt)
alpha.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$`Alpha (ug/g)`, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
# Orange significant different from others
# Red, White, and Yellow not significantly different

alpha.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=`Alpha (ug/g)`, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
  theme_classic()
  alpha.box.1 <- alpha.box + 
                           annotate(geom="text",
                                      x="Orange",
                                      y=125,
                                      label="a",
                                    size=16) +
                           annotate(geom="text", 
                                       x="Red",
                                       y=125,
                                       label="b",
                                    size=16) +
                           annotate(geom="text",
                                     x="White",
                                     y=125,
                                     label="b",
                                    size=16) +
                           annotate(geom="text",
                                       x="Yellow",
                                       y=125,
                                       label="b",
                                    size=16)


  alpha.box.2 <- alpha.box.1 + theme(text = element_text(size=36)) +
                                theme(legend.position = "none")
  
####### (4.2) Beta-carotene distribution & mean separation ######  
  beta.np <- kruskal.test(`Beta (ug/g)`~ Color, data = Pheno.dat.Filt)
  beta.np$p.value # pvalue = 2.903925e-81
  pairwise.wilcox.test(Pheno.dat.Filt$`Beta (ug/g)`, 
                       Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  beta.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=`Beta (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  beta.box.1 <- beta.box + 
    annotate(geom="text",
             x="Orange",
             y=1500,
             label="a",
             size=16) +
    annotate(geom="text", 
             x="Red",
             y=1500,
             label="b",
             size=16) +
    annotate(geom="text",
             x="White",
             y=1500,
             label="c",
             size=16)+
    annotate(geom="text",
             x="Yellow",
             y=1500,
             label="d",
             size=16)
  
  
  beta.box.2 <- beta.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### (4.3) Lutein distribution & mean separation #######  
  
  Lut.np <- kruskal.test(`Lut (ug/g)`~ Color, data = Pheno.dat.Filt)
  Lut.np$p.value # pvalue = 3.218703e-09
  pairwise.wilcox.test(Pheno.dat.Filt$`Lut (ug/g)`, 
                       Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  lut.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=`Lut (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  lut.box.1 <- lut.box + 
    annotate(geom="text",
             x="Orange",
             y=150,
             label="a",
             size=16) +
    annotate(geom="text", 
             x="Red",
             y=150,
             label="a",
             size=16) +
    annotate(geom="text",
             x="White",
             y=150,
             label="b",
             size=16) +
    annotate(geom="text",
             x="Yellow",
             y=150,
             label="a",
             size=16)
  
  
  lut.box.2 <- lut.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
  
  
####### (4.4) Lycopene distribution & mean separation ######  
  Lyco.np <- kruskal.test(`Lycopene (ug/g)`~ Color, data = Pheno.dat.Filt)
  Lyco.np$p.value # pvalue = 5.804908e-19
  pairwise.wilcox.test(Pheno.dat.Filt$`Lycopene (ug/g)`, 
                       Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  lyco.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=`Lycopene (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  lyco.box.1 <- lyco.box + 
    annotate(geom="text",
             x="Orange",
             y=750,
             label="a",
             size=16) +
    annotate(geom="text", 
             x="Red",
             y=750,
             label="b",
             size=16) +
    annotate(geom="text",
             x="White",
             y=750,
             label="c",
             size=16) +
    annotate(geom="text",
             x="Yellow",
             y=750,
             label="c",
             size=16)
  
  
  lyco.box.2 <- lyco.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### (4.5) Phyotoene distribution & mean separation #######  
  Phyto.np <- kruskal.test(`Phyto (ug/g)`~ Color, data = Pheno.dat.Filt)
  Phyto.np$p.value # pvalue = 6.631561e-66
  pairwise.wilcox.test(Pheno.dat.Filt$`Phyto (ug/g)`, 
                       Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  phyto.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=`Phyto (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  phyto.box.1 <- phyto.box + 
    annotate(geom="text",
             x="Orange",
             y=200,
             label="a",
             size = 16) +
    annotate(geom="text", 
             x="Red",
             y=200,
             label="a",
             size = 16) +
    annotate(geom="text",
             x="White",
             y=200,
             label="b",
             size = 16) +
    annotate(geom="text",
             x="Yellow",
             y=200,
             label="b",
             size = 16)
  
  
  phyto.box.2 <- phyto.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### (4.6) Zeta distribution & mean separation #######   
  zeta.np <- kruskal.test(`Zeta (ug/g)`~ Color, data = Pheno.dat.Filt)
  zeta.np$p.value # pvalue = 2.107543e-08
  pairwise.wilcox.test(Pheno.dat.Filt$`Zeta (ug/g)`, 
                       Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all groups significantly different 
  zeta.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=`Zeta (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  zeta.box.1 <- zeta.box + 
    annotate(geom="text",
             x="Orange",
             y=100,
             label="a",
             size=16) +
    annotate(geom="text", 
             x="Red",
             y=100,
             label="b",
             size=16) +
    annotate(geom="text",
             x="White",
             y=100,
             label="b",
             size=16) +
    annotate(geom="text",
             x="Yellow",
             y=100,
             label="b",
             size=16)
  
  
  zeta.box.2 <- zeta.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
  
####### (4.7) Total distribution & mean separation ######

  Total.np <- kruskal.test(`Total (ug/g)`~ Color, data = Pheno.dat.Filt)
  Total.np$p.value # pvalue = 8.43771e-80
  pairwise.wilcox.test( Pheno.dat.Filt$`Total (ug/g)`, 
                        Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all < 0.05
  # all groups significantly different 
  Total.box <- ggplot(data=  Pheno.dat.Filt, aes(x=Color, y=`Total (ug/g)`, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  Total.box.1 <- Total.box + 
    annotate(geom="text",
             x="Orange",
             y=2500,
             label="a",
             size=16) +
    annotate(geom="text", 
             x="Red",
             y=2500,
             label="b",
             size=16) +
    annotate(geom="text",
             x="White",
             y=2500,
             label="c",
             size=16) +
    annotate(geom="text",
             x="Yellow",
             y=2500,
             label="d",
             size=16)
  
  
  Total.box.2 <- Total.box.1 + theme(text = element_text(size=36)) +
    theme(legend.position = "none")
####### (4.8) Ab.ratio distributions & separation #######  
  AB.np <- kruskal.test(AB.Ratio ~ Color, data = Pheno.dat.Filt)
  AB.np$p.value # pvalue = 5.638406e-72
  pairwise.wilcox.test(Pheno.dat.Filt$AB.Ratio, 
                       Pheno.dat.Filt$Color, 
                       p.adjust.method = "BH") 
  # all < 0.05
  # all groups significantly different 
  AB.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=AB.Ratio, fill=Color))+
    geom_boxplot()+
    scale_fill_manual(values=c("Orange", "Red", "Grey", "Yellow")) + 
    theme_classic()
  AB.box.1 <- AB.box + 
    annotate(geom="text",
             x="Orange",
             y=25,
             label="a",
             size=16) +
    annotate(geom="text", 
             x="Red",
             y=25,
             label="b",
             size=16) +
    annotate(geom="text",
             x="White",
             y=25,
             label="c",
             size=16) +
    annotate(geom="text",
             x="Yellow",
             y=25,
             label="d",
             size=16)
  
  
  AB.box.2 <- AB.box.1 + theme(text = element_text(size=36))   +
    theme(legend.position = "none")
  
####### Part Five: Write Output #######  
# (5.1) write the phenotype file 
  write.csv(x=Pheno.dat.Filt, file="AllColPhenoFilt.csv",
            row.names=F)
  
# (5.2) Add (a) to annotate figure
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

  
  
####### Part Six: Focus on orange accessions only. ####### 
rm(list = ls())
dat <- read_csv("FilteredOrangeDat.csv", col_names=T)
  
mean(dat$`Alpha (ug/g)`, na.rm=T) # 329.3875
sd(dat$`Alpha (ug/g)`, na.rm=T)  # 249.7071
mean(dat$`Beta (ug/g)`, na.rm=T) # 895.2632
sd(dat$`Beta (ug/g)`, na.rm=T) # 576.8396
mean(dat$`Lut (ug/g)`, na.rm=T) # 89.36006
sd(dat$`Lut (ug/g)`, na.rm=T) #  48.85191
mean(dat$`Lycopene (ug/g)`, na.rm=T)  # 52.04231
sd(dat$`Lycopene (ug/g)`, na.rm=T) #32.48107
mean(dat$`Phyto (ug/g)`, na.rm=T) # 123.3972
sd(dat$`Phyto (ug/g)`, na.rm=T) # 97.30621
mean(dat$`Zeta (ug/g)`, na.rm=T) # 12.28356
sd(dat$`Zeta (ug/g)`, na.rm=T) # 30.32293
mean(dat$`Total (ug/g)`, na.rm=T) # 1363.661
sd(dat$`Total (ug/g)`, na.rm=T) # 879.505

####### Part Seven Correlation Between Traits ##### 
cor.test(dat$`Alpha (ug/g)`, dat$`Beta (ug/g)`)
cor.test(dat$`Alpha (ug/g)`, dat$`Lut (ug/g)`)
cor.test(dat$`Alpha (ug/g)`, dat$`Lycopene (ug/g)`)
cor.test(dat$`Alpha (ug/g)`, dat$`Phyto (ug/g)`)
cor.test(dat$`Alpha (ug/g)`, dat$`Zeta (ug/g)`)
cor.test(dat$`Alpha (ug/g)`, dat$`Total (ug/g)`)


cor.test(dat$`Beta (ug/g)`, dat$`Lut (ug/g)`)
cor.test(dat$`Beta (ug/g)`, dat$`Lycopene (ug/g)`)
cor.test(dat$`Beta (ug/g)`, dat$`Phyto (ug/g)`)
cor.test(dat$`Beta (ug/g)`, dat$`Zeta (ug/g)`)
cor.test(dat$`Beta (ug/g)`, dat$`Total (ug/g)`)


cor.test(dat$`Lut (ug/g)`, dat$`Lycopene (ug/g)`)
cor.test(dat$`Lut (ug/g)`, dat$`Phyto (ug/g)`)
cor.test(dat$`Lut (ug/g)`, dat$`Zeta (ug/g)`)
cor.test(dat$`Lut (ug/g)`, dat$`Total (ug/g)`)

cor.test(dat$`Lycopene (ug/g)`, dat$`Phyto (ug/g)`)
cor.test(dat$`Lycopene (ug/g)`, dat$`Zeta (ug/g)`)
cor.test(dat$`Lycopene (ug/g)`, dat$`Total (ug/g)`)


cor.test(dat$`Phyto (ug/g)`, dat$`Zeta (ug/g)`)
cor.test(dat$`Phyto (ug/g)`, dat$`Total (ug/g)`)

cor.test(dat$`Zeta (ug/g)`, dat$`Total (ug/g)`)

####### Part Eight Compare Core Color amounts  ###### 

# (4.2) Non-parametric test - alpha
# Determine if core color represents carotenoid accumulation
alpha.np <- kruskal.test(`Alpha (ug/g)` ~ `Inner Color`, data = dat)
alpha.np$p.value #5.666231e-06
pairwise.wilcox.test(dat$`Alpha (ug/g)`, 
                     dat$`Inner Color`, 
                     p.adjust.method = "BH") 
# (4.2) Non-parametric test - beta
beta.np <- kruskal.test(`Beta (ug/g)`~ `Inner Color`, data = dat)
beta.np$p.value # 3.924782e-06
 
# (4.3) Non-parametric test - lutein
lut.np <- kruskal.test(`Lut (ug/g)` ~ `Inner Color`, data = dat)
lut.np$p.value # 0.02
# (4.4) Non-parametric test - lycopene
lyco.np <- kruskal.test(`Lycopene (ug/g)` ~ `Inner Color`, data = dat)
lyco.np$p.value # 0.0001737949
# (4.4) Non-parametric test - phyotene
phyto.np <- kruskal.test(`Phyto (ug/g)` ~ `Inner Color`, data = dat)
phyto.np$p.value # 0.22
# (4.4) Non-parametric test - zeta
zeta.np <- kruskal.test(`Zeta (ug/g)` ~ `Inner Color`, data = dat)
zeta.np$p.value # 0.72
# (4.4) Non-parametric test - total
Tot.np <- kruskal.test(`Total (ug/g)` ~ `Inner Color`, data = dat)
Tot.np$p.value # 1.456005e-05


####### Part Nine: Create Box plots to present in supplemental ###### 
dat.1 <- dat[-which(is.na(dat$`Inner Color`)==T),]



# (4.5) . 

beta.box <- ggplot(data=dat.1, 
                   aes(x=`Inner Color`, y=`Beta (ug/g)`, fill=`Inner Color`))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

beta.1 <- beta.box + labs(x="Inner Color", y="Beta(ug/g)")


alpha.box <- ggplot(data=dat.1, 
                    aes(x=`Inner Color`, y=`Alpha (ug/g)`, fill=`Inner Color`))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

alpha.1 <-  alpha.box + labs(x="Inner Color", y="Alpha(ug/g)")


# (4.5) Create Box plots to present in supplemental. 

lut.box <- ggplot(data=dat.1, 
                   aes(x=`Inner Color`, y=`Alpha (ug/g)`, fill=`Inner Color`))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

lut.box.1 <- beta.box + labs(x="Inner Color", y="Lutein(ug/g)")


lyco.box <- ggplot(data=dat, 
                    aes(x=`Inner Color`, y=`Lycopene (ug/g)`, fill=`Inner Color`))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

lyco.box.1 <-  alpha.box + labs(x="Inner Color", y="Lycopene(ug/g)")




# (4.6) Create Multipane Figure
Alpha <- arrangeGrob(alpha.1, top = textGrob("c", x = unit(0, "npc")
                                             , y   = unit(1, "npc"), just=c("right","top"),
                                             gp=gpar(col="black", fontsize=36, 
                                                     fontfamily="Times Roman")))
Beta <- arrangeGrob(beta.1, top = textGrob("d", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("right","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                     fontfamily="Times Roman")))

Lut <- arrangeGrob(lut.box.1, top = textGrob("e", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("right","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))
Lyco <- arrangeGrob(lyco.box.1, top = textGrob("f", x = unit(0, "npc")
                                             , y   = unit(1, "npc"), just=c("right","top"),
                                             gp=gpar(col="black", fontsize=36, 
                                                     fontfamily="Times Roman")))

SupFigBoxColor <- ggarrange(Alpha, Beta, Lut, Lyco)

# (4.7)  Write to file!
png(filename = "CA19SupFigBoxCoreColor.png", width=1950, height = 1300)
SupFigBoxColor
dev.off()

#