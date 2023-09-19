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
Dat.1 <- read.csv("../Analysis_4_FigsHPLCdat_WI18/HPLC_dataWI18MU.csv",
                  head = T)

CLR <- read_csv("WI19ColorScores.csv", col_names = T)

######## Part Two: Identify Outliers ####### 
# (2.1) identify outliers for alpha
Outliers <- Dat.1 %>%
  identify_outliers(Alpha..ug.g.)
# (2.2) identify extreme outliers
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1 
# # (2.3 )Remove extreme outliers
Dat.2 <- Dat.1[-which(Dat.1$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (2.4) Repeat for Beta
Outliers <- Dat.1 %>%
  identify_outliers(Beta..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
                     %in% Outliers.1$`Sample.Name`))
#(2.5)Repeat for Lutein
Outliers <- Dat.1 %>%
  identify_outliers(Lut..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
             %in% Outliers.1$`Sample.Name`))
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (2.6) Repeat for Lycopene
Outliers <- Dat.1 %>%
  identify_outliers(Lycopene..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
             %in% Outliers.1$`Sample.Name`))
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (2.7) Repeat for Phytoene
Outliers <- Dat.1 %>%
  identify_outliers(Phyto..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (2.8) Repeat for Zeta
Outliers <- Dat.1 %>%
  identify_outliers(Zeta..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]

# (2.9) Repeat for Total
Outliers <- Dat.1 %>%
  identify_outliers(Total..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
             %in% Outliers.1$`Sample.Name`))

# (2.10) write output to file for subsequent analyses
write.csv(x=Dat.2, 
          file="WI18PhenoNoOutLie.csv", 
          row.names=F)

######## Part Three: Calculate Mean and Sd ####### 
mean(Dat.2$Alpha..ug.g., na.rm=T) # 177.0747
sd(Dat.2$Alpha..ug.g., na.rm=T)  # 188.5045
mean(Dat.2$Beta..ug.g., na.rm=T) # 290.6945
sd(Dat.2$Beta..ug.g., na.rm=T) # 264.5246
mean(Dat.2$Lut..ug.g., na.rm=T) # 54.78286
sd(Dat.2$Lut..ug.g., na.rm=T) # 48.58097
mean(Dat.2$Lycopene..ug.g., na.rm=T) # 21.79776
sd(Dat.2$Lycopene..ug.g., na.rm=T) # 15.68609
mean(Dat.2$Phyto..ug.g., na.rm=T)  # 10.8824
sd(Dat.2$Phyto..ug.g., na.rm=T) # 15.40891
mean(Dat.2$Zeta..ug.g., na.rm=T) # 12.30215
sd(Dat.2$Zeta..ug.g., na.rm=T) # 20.41982
mean(Dat.2$Total..ug.g., na.rm=T) # 459.4159
sd(Dat.2$Total..ug.g., na.rm=T) # 430.0208
 
######## Part Four: Nonparametric - Core Color ####### 
# (4.1) Add Color Scores to Phenotype
CLR.1 <- CLR[which(CLR$`Sample Name` %in% Dat.2$Sample.Name),]
colnames(CLR.1)[1] <- "Sample.Name"
Pheno.dat <- merge(Dat.2, CLR.1, by="Sample.Name")

# (4.2) Non-parametric test - alpha
  # Determine if core color represents carotenoid accumulation
alpha.np <- kruskal.test(Alpha..ug.g. ~ Inner.Color, data = Pheno.dat)
alpha.np$p.value # 0.01345261
pairwise.wilcox.test(Pheno.dat$Alpha..ug.g., 
                     Pheno.dat$Inner.Color, 
                     p.adjust.method = "BH") 
# (4.2) Non-parametric test - beta
beta.np <- kruskal.test(Beta..ug.g. ~ Inner.Color, data = Pheno.dat)
beta.np$p.value # 0.010
pairwise.wilcox.test(Pheno.dat$Beta..ug.g., 
                     Pheno.dat$Inner.Color, 
                     p.adjust.method = "BH") 
# (4.3) Non-parametric test - lutein
lut.np <- kruskal.test(Lut..ug.g. ~ Inner.Color, data = Pheno.dat)
lut.np$p.value # 0.22
# (4.4) Non-parametric test - lycopene
lyco.np <- kruskal.test(Lycopene..ug.g. ~ Inner.Color, data = Pheno.dat)
lyco.np$p.value # 0.42
# (4.4) Non-parametric test - phyotene
phyto.np <- kruskal.test(Phyto..ug.g. ~ Inner.Color, data = Pheno.dat)
phyto.np$p.value # 0.27
# (4.4) Non-parametric test - zeta
zeta.np <- kruskal.test(Zeta..ug.g. ~ Inner.Color, data = Pheno.dat)
zeta.np$p.value # 0.57
# (4.4) Non-parametric test - total
Tot.np <- kruskal.test(Total..ug.g. ~ Inner.Color, data = Pheno.dat)
Tot.np$p.value # 0.058

# (4.5) Create Box plots to present in supplemental. 

beta.box <- ggplot(data=Pheno.dat, 
                   aes(x=Inner.Color, y=Beta..ug.g., fill=Inner.Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

beta.1 <- beta.box + labs(x="Inner Color", y="Beta(ug/g)")


alpha.box <- ggplot(data=Pheno.dat, 
                   aes(x=Inner.Color, y=Alpha..ug.g., fill=Inner.Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

alpha.1 <-  alpha.box + labs(x="Inner Color", y="Alpha(ug/g)")

# (4.6) Create Multipane Figure
Alpha <- arrangeGrob(alpha.1, top = textGrob("(a)", x = unit(0, "npc")
                                           , y   = unit(1, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=36, 
                                                   fontfamily="Times Roman")))
Beta <- arrangeGrob(beta.1, top = textGrob("(b)", x = unit(0, "npc")
                                         , y   = unit(1, "npc"), just=c("left","top"),
                                         gp=gpar(col="black", fontsize=36, 
                                                 fontfamily="Times Roman")))
SupFigBoxColor <- ggarrange(Alpha, Beta)

# (4.7)  Write to file!
png(filename = "WI18SupFigBoxColor.png", width=1950, height = 650)
SupFigBoxColor
dev.off()


######## Part Five: Create a file for GWA ####### 
# (5.1) Remove unneeded columns
Pheno.dat <- Pheno.dat[,-10]
Pheno.dat <- Pheno.dat[,-13]
#(5.2) Change column names to avoid format issues
colnames(Pheno.dat) <- c("Sample.Name",
                         "Alpha", 
                         "Beta", 
                         "Lutein",
                         "Lycopene",
                         "Phytoene",
                         "Zeta", 
                         "Total",
                         "AB.ratio",
                         "Outer.Color",
                         "Inner.Color",
                         "Color")

# (5.3) Change Outer Colors to Numeric
Pheno.dat$Outer.Color[which(Pheno.dat$Outer.Color == "Orange")] <- 1
Pheno.dat$Outer.Color[which(Pheno.dat$Outer.Color == "Yellow")] <- 2
Pheno.dat$Outer.Color[which(Pheno.dat$Outer.Color == "Red")] <- 3
Pheno.dat$Outer.Color[which(Pheno.dat$Outer.Color == "Purple")] <- NA

# (5.4) Change Inner Colors to Numeric
Pheno.dat$Inner.Color[which(Pheno.dat$Inner.Color == "Orange")] <- 1
Pheno.dat$Inner.Color[which(Pheno.dat$Inner.Color == "Yellow")] <- 2
Pheno.dat$Inner.Color[which(Pheno.dat$Inner.Color == "Red")] <- 3
Pheno.dat$Inner.Color[which(Pheno.dat$Inner.Color == "Purple")] <- NA
# (5.5) Write to file! 
write.csv(x=Pheno.dat,
          file="WI18.GWA.pheno.csv",
          row.names=F)



######## Part Six: Correlations of traits #######
# (6.1) Alpha Correlations
cor.test(Pheno.dat$Alpha, Pheno.dat$Beta)
cor.test(Pheno.dat$Alpha, Pheno.dat$Lutein)
cor.test(Pheno.dat$Alpha, Pheno.dat$Lycopene)
cor.test(Pheno.dat$Alpha, Pheno.dat$Phytoene)
cor.test(Pheno.dat$Alpha, Pheno.dat$Zeta)
cor.test(Pheno.dat$Alpha, Pheno.dat$Total)
# (6.2) Beta Correlations
cor.test(Pheno.dat$Beta, Pheno.dat$Lutein)
cor.test(Pheno.dat$Beta, Pheno.dat$Lycopene)
cor.test(Pheno.dat$Beta, Pheno.dat$Phytoene)
cor.test(Pheno.dat$Beta, Pheno.dat$Zeta)
cor.test(Pheno.dat$Beta, Pheno.dat$Total)
# (6.3) Lutein Correlations
cor.test(Pheno.dat$Lutein, Pheno.dat$Lycopene)
cor.test(Pheno.dat$Lutein, Pheno.dat$Phytoene)
cor.test(Pheno.dat$Lutein, Pheno.dat$Zeta)
cor.test(Pheno.dat$Lutein, Pheno.dat$Total)
# (6.4) Lycopene Correlations
cor.test(Pheno.dat$Lycopene, Pheno.dat$Phytoene)
cor.test(Pheno.dat$Lycopene, Pheno.dat$Zeta)
cor.test(Pheno.dat$Lycopene, Pheno.dat$Total)
# (6.5) Phytoene Correlations
cor.test(Pheno.dat$Phytoene, Pheno.dat$Zeta)
cor.test(Pheno.dat$Phytoene, Pheno.dat$Total)
# (6.5) Zeta Correlations
cor.test(Pheno.dat$Zeta, Pheno.dat$Total)
# Manually Entered into First Table
# The end :) 
