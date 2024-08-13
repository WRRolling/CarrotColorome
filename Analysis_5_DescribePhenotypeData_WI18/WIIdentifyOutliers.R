######## Part One:Prepare R Environment & Load Phenotypic Data ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Load libraries
library('tidyverse')
library(rstatix)
library(grid)
library(gridExtra)
library(ggpubr)
# (1.2) Load phenotype Data
HPLC.dat <- read.csv("../Analysis_4_FigsHPLCdat_WI18/HPLC_dataWI18MU2.csv",
                  head = T)
# (1.3) Load Color Data
CLR <- read_csv("WI18ColorScores2.csv", col_names = T)

# (1.4) Subset data to only include those samples with color and HPLC
Color.dat.sub <- CLR[which(CLR$Number %in% HPLC.dat$Sample.Name),]
HPLC.dat.sub <- HPLC.dat[which(HPLC.dat$Sample.Name %in% Color.dat.sub$Number),]
Color.dat.sub <- Color.dat.sub[order(Color.dat.sub$Number),]
HPLC.dat.sub <-  HPLC.dat.sub[order(HPLC.dat.sub$Sample.Name),]
Pheno.dat <- cbind(Color.dat.sub, HPLC.dat.sub)
Pheno.dat <- Pheno.dat[,-6]
colnames(Pheno.dat)[1] <- "Sample.Name"

# (1.5) Clean R environment
rm(CLR,
   Color.dat.sub,
   HPLC.dat,
   HPLC.dat.sub)

# (1.6) Create AB ratio and total carotenoid phenotypes
Pheno.dat[,12] <- (Pheno.dat$Alpha..ug.g.+Pheno.dat$Beta..ug.g.)/Pheno.dat$Lut..ug.g.
colnames(Pheno.dat)[12] <- "AB.Ratio"
Pheno.dat[,13] <- Pheno.dat$Alpha..ug.g.+
                  Pheno.dat$Beta..ug.g. +
                  Pheno.dat$Phyto..ug.g. +
                  Pheno.dat$Zeta..ug.g.
colnames(Pheno.dat)[13] <- "Total"



######## Part Two:Calculate Ratio of alpha + Beta relative to lutein #######
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Orange")], na.rm=T) # 9.9
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Yellow")], na.rm=T) # 0.0
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Red")], na.rm=T) # .13
mean(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "White")], na.rm=T) # 0.16
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Orange")], na.rm=T) # 8.5
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Yellow")], na.rm=T) # 0.05
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "Red")], na.rm=T) # .18
sd(Pheno.dat$AB.Ratio[which(Pheno.dat$Color == "White")], na.rm=T) # 0.09

######## Part Three: Identify Outliers in Orange Samples ####### 
### (3.1) Subset to Orange Data
Dat.1 <- Pheno.dat[Pheno.dat$Color == "Orange",]
# (3.2) identify outliers for alpha
Outliers <- Dat.1 %>%
  identify_outliers(Alpha..ug.g.)
# (3.3) identify extreme outliers
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1 
# (3.4) Remove extreme outliers
Dat.2 <- Dat.1[-which(Dat.1$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (3.5) Repeat for Beta
Outliers <- Dat.1 %>%
  identify_outliers(Beta..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
                     %in% Outliers.1$`Sample.Name`))
#(3.6)Repeat for Lutein
Outliers <- Dat.1 %>%
  identify_outliers(Lut..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
             %in% Outliers.1$`Sample.Name`))
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (3.7) Repeat for Lycopene
Outliers <- Dat.1 %>%
  identify_outliers(Lycopene..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
             %in% Outliers.1$`Sample.Name`))
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (3.8) Repeat for Phytoene
Outliers <- Dat.1 %>%
  identify_outliers(Phyto..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]
# (3.9) Repeat for Zeta
Outliers <- Dat.1 %>%
  identify_outliers(Zeta..ug.g.)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
Dat.2 <- Dat.2[-which(Dat.2$`Sample.Name` 
                      %in% Outliers.1$`Sample.Name`),]

# (3.10) Repeat for Total
Outliers <- Dat.1 %>%
  identify_outliers(Total)
Outliers.1 <- Outliers[which(Outliers$is.extreme == T),]
Outliers.1
length(which(Dat.2$`Sample.Name`
             %in% Outliers.1$`Sample.Name`))

# (3.11) write output to file for subsequent analyses
write.csv(x=Dat.2, 
          file="WI18PhenoNoOutLie.csv", 
          row.names=F)






######## Part Four: Calculate Mean and Sd ####### 
mean(Dat.2$Alpha..ug.g., na.rm=T) # 167.76
sd(Dat.2$Alpha..ug.g., na.rm=T)  # 186.43
mean(Dat.2$Beta..ug.g., na.rm=T) # 281.1
sd(Dat.2$Beta..ug.g., na.rm=T) # 265.72
mean(Dat.2$Lut..ug.g., na.rm=T) # 56.23
sd(Dat.2$Lut..ug.g., na.rm=T) # 46.9
mean(Dat.2$Lycopene..ug.g., na.rm=T) # 20.75
sd(Dat.2$Lycopene..ug.g., na.rm=T) # 16.91519
mean(Dat.2$Phyto..ug.g., na.rm=T)  # 10.34585
sd(Dat.2$Phyto..ug.g., na.rm=T) # 15.13154
mean(Dat.2$Zeta..ug.g., na.rm=T) # 11.71479
sd(Dat.2$Zeta..ug.g., na.rm=T) # 20.08579
mean(Dat.2$Total, na.rm=T) # 435.0628
sd(Dat.2$Total, na.rm=T) # 423.7077
 
######## Part Five: Nonparametric - Core Color ####### 
# (5.1) Non-parametric test - alpha
  # Determine if core color represents carotenoid accumulation
alpha.np <- kruskal.test(Alpha..ug.g. ~ Interior, data = Dat.2)
alpha.np$p.value # 8.261558e-19

# (5.2) Non-parametric test - beta
beta.np <- kruskal.test(Beta..ug.g. ~ Interior, data =  Dat.2)
beta.np$p.value # 1.663683e-19

# (5.3) Non-parametric test - lutein
lut.np <- kruskal.test(Lut..ug.g. ~ Interior, data =  Dat.2)
lut.np$p.value # 0.001620314
# (5.4) Non-parametric test - lycopene
lyco.np <- kruskal.test(Lycopene..ug.g. ~ Interior, data =  Dat.2)
lyco.np$p.value # 0.06719715
# (5.5) Non-parametric test - phyotene
phyto.np <- kruskal.test(Phyto..ug.g. ~ Interior, data =  Dat.2)
phyto.np$p.value #  3.857733e-07
# (5.6) Non-parametric test - zeta
zeta.np <- kruskal.test(Zeta..ug.g. ~ Interior, data =  Dat.2)
zeta.np$p.value # 0.0003160275
# (5.7) Non-parametric test - total
Tot.np <- kruskal.test(Total ~ Interior, data =  Dat.2)
Tot.np$p.value # 4.636416e-18



# (4.5) Create Box plots to present in supplemental. 

beta.box <- ggplot(data=Orange.dat, 
                   aes(x=Interior, y=Beta..ug.g., fill=Interior))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Yellow")) + 
  theme(text = element_text(size = 20))

beta.1 <- beta.box + labs(x="Inner Color", y="Beta(ug/g)")


alpha.box <- ggplot(data=Orange.dat, 
                   aes(x=Interior, y=Alpha..ug.g., fill=Interior))+
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


######## Part Six: Create a file for GWA ####### 
# (6.1) Remove unneeded columns
Orange.dat <- Dat.2[,-2]
#(6.2) Change column names to avoid format issues
colnames(Orange.dat) <- c("Sample.Name",
                         "Outer.Color",
                         "Inner.Color",
                         "Color",
                         "Alpha", 
                         "Beta", 
                         "Lutein",
                         "Lycopene",
                         "Phytoene",
                         "Zeta", 
                         "AB.ratio",
                         "Total")

# (6.3) Change Outer Colors to Numeric
Orange.dat$Outer.Color[which(Orange.dat$Outer.Color == "Orange")] <- 1
Orange.dat$Outer.Color[which(Orange.dat$Outer.Color == "Yellow")] <- 2
Orange.dat$Outer.Color[which(Orange.dat$Outer.Color == "Red")] <- NA
Orange.dat$Outer.Color[which(Orange.dat$Outer.Color == "White")] <- NA

# (6.4) Change Inner Colors to Numeric
Orange.dat$Inner.Color[which(Orange.dat$Inner.Color == "Orange")] <- 1
Orange.dat$Inner.Color[which(Orange.dat$Inner.Color == "Yellow")] <- 2
Orange.dat$Inner.Color[which(Orange.dat$Inner.Color == "Red")] <- NA
Orange.dat$Inner.Color[which(Orange.dat$Inner.Color == "White")] <- NA
# (6.5) Write to file! 
write.csv(x=Pheno.dat,
          file="WI18.HPLC.GWA.pheno.csv",
          row.names=F)



######## Part Seven: Correlations of traits #######
# (7.1) Alpha Correlations
cor.test(Orange.dat$Alpha, Orange.dat$Beta)
cor.test(Orange.dat$Alpha, Orange.dat$Lutein)
cor.test(Orange.dat$Alpha, Orange.dat$Lycopene)
cor.test(Orange.dat$Alpha, Orange.dat$Phytoene)
cor.test(Orange.dat$Alpha, Orange.dat$Zeta)
cor.test(Orange.dat$Alpha, Orange.dat$Total)
# (7.2) Beta Correlations
cor.test(Orange.dat$Beta, Orange.dat$Lutein)
cor.test(Orange.dat$Beta, Orange.dat$Lycopene)
cor.test(Orange.dat$Beta, Orange.dat$Phytoene)
cor.test(Orange.dat$Beta, Orange.dat$Zeta)
cor.test(Orange.dat$Beta, Orange.dat$Total)
# (7.3) Lutein Correlations
cor.test(Orange.dat$Lutein, Orange.dat$Lycopene)
cor.test(Orange.dat$Lutein, Orange.dat$Phytoene)
cor.test(Orange.dat$Lutein, Orange.dat$Zeta)
cor.test(Orange.dat$Lutein, Orange.dat$Total)
# (7.4) Lycopene Correlations
cor.test(Orange.dat$Lycopene, Orange.dat$Phytoene)
cor.test(Orange.dat$Lycopene, Orange.dat$Zeta)
cor.test(Orange.dat$Lycopene, Orange.dat$Total)
# (7.5) Phytoene Correlations
cor.test(Orange.dat$Phytoene, Orange.dat$Zeta)
cor.test(Orange.dat$Phytoene, Orange.dat$Total)
# (7.6) Zeta Correlations
cor.test(Orange.dat$Zeta, Orange.dat$Total)
# Manually Entered into First Table
# The end :) 
















######## Part Eight: clean R environment and take Notes #######
rm(phyto.np)
rm(zeta.np)
rm(lyco.np)
rm(lut.np)
rm(beta.np)
rm(alpha.np)
rm(Outliers)
rm(Outliers.1)
rm(Tot.np)

# (8.2) Create List of Filtered Samples
Filtered.Samples <- c(80209, 80320,80325,80326,80334,80376,80424,
                      80426, 80451, 80566,80579,80583,80589,80601,
                      80624,80648,80675,80696,80697,80732,80735,
                      80744,80748,80749,80753,80760,80762,80765,80782)

######## Part Nine: Mean Separation and histograms all colors #######
# (9.1) Remove Outliers
Pheno.dat.Filt <-  Pheno.dat[-which(Pheno.dat$Sample.Name %in% 
                                    Filtered.Samples),]
colnames(Pheno.dat.Filt) <- c( "Sample.Name", "ID", "Exterior.Color", "Interior.Color",
                               "Color", "Alpha", "Beta", "Lutein", "Lycopene",
                               "Phytoene", "Zeta", "AB.Ratio", "Total")
# (9.2) Remove Any sample without a color score, or red score
Pheno.dat.Filt <- Pheno.dat.Filt[-which(is.na(Pheno.dat.Filt$Color)==T),]
Pheno.dat.Filt <- Pheno.dat.Filt[-which(Pheno.dat.Filt$Color == "Red"),]

# (9.3) Complete Non-parametric mean separation for alpha
alpha.np <- kruskal.test(Alpha ~ Color, data = Pheno.dat.Filt)
alpha.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Alpha, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
# (9.4) Create Boxplot of alpha data
alpha.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Alpha, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
alpha.box.1 <- alpha.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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

######## Part Ten,: Color Mean Separation - Beta #######
beta.np <- kruskal.test(Beta ~ Color, data = Pheno.dat.Filt)
beta.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Beta, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
beta.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Beta, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
beta.box.1 <- beta.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
beta.box.2 <- beta.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")

######## Part Eleven: Color Mean Separation - Lutein #######

Lutein.np <- kruskal.test(Lutein ~ Color, data = Pheno.dat.Filt)
Lutein.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Lutein, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
Lutein.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Lutein, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
Lutein.box.1 <- Lutein.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
Lutein.box.2 <- Lutein.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")

######## Part Twelve: Color Mean Separation - Lycopene #######

Lycopene.np <- kruskal.test(Lycopene ~ Color, data = Pheno.dat.Filt)
Lycopene.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Lycopene, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
Lycopene.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Lycopene, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
Lycopene.box.1 <- Lycopene.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
Lycopene.box.2 <-Lycopene.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")

######## Part Thirteen: Color Mean Separation - Phyto #######

Phytoene.np <- kruskal.test(Phytoene ~ Color, data = Pheno.dat.Filt)
Phytoene.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Phytoene, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
# (7.5) Create Boxplot of alpha data
Phytoene.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Phytoene, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
Phytoene.box.1 <- Phytoene.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
Phytoene.box.2 <-Phytoene.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")

######## Part Fourteen: Color Mean Separation - Zeta #######

Zeta.np <- kruskal.test(Zeta ~ Color, data = Pheno.dat.Filt)
Zeta.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Zeta, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
Zeta.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Zeta, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
Zeta.box.1 <-Zeta.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
Zeta.box.2 <-Zeta.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")


######## Part Fifteen: Color Mean Separation - Total #######
Total.np <- kruskal.test(Total ~ Color, data = Pheno.dat.Filt)
Total.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$Total, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
Total.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=Total, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
Total.box.1 <-Total.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
Total.box.2 <-Total.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")

######## Part Sixteen: Color Mean Separation - AB.Ratio #######
AB.ratio.np <- kruskal.test(AB.Ratio ~ Color, data = Pheno.dat.Filt)
AB.ratio.np$p.value # pvalue = 2.274252e-82
pairwise.wilcox.test(Pheno.dat.Filt$AB.Ratio, 
                     Pheno.dat.Filt$Color, 
                     p.adjust.method = "BH") 
AB.ratio.box <- ggplot(data= Pheno.dat.Filt, aes(x=Color, y=AB.Ratio, fill=Color))+
  geom_boxplot()+
  scale_fill_manual(values=c("Orange", "Grey", "Yellow")) + 
  theme_classic()
AB.ratio.box.1 <-AB.ratio.box + 
  annotate(geom="text",
           x="Orange",
           y=125,
           label="a",
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
AB.ratio.box.2 <-AB.ratio.box.1 + theme(text = element_text(size=36)) +
  theme(legend.position = "none")



######## Part Seventeen: Create Plot for Manuscript #######
alpha.box.3 <- arrangeGrob(alpha.box.2, top = textGrob("(a)", x = unit(0, "npc")
                                                       , y   = unit(1, "npc"), just=c("left","top"),
                                                       gp=gpar(col="black", fontsize=36, 
                                                               fontfamily="Times Roman")))
beta.box.3 <- arrangeGrob(beta.box.2 , top = textGrob("(b)", x = unit(0, "npc")
                                                      , y   = unit(1, "npc"), just=c("left","top"),
                                                      gp=gpar(col="black", fontsize=36, 
                                                              fontfamily="Times Roman")))
lut.box.3 <- arrangeGrob(Lutein.box.2, top = textGrob("(c)", x = unit(0, "npc")
                                                   , y   = unit(1, "npc"), just=c("left","top"),
                                                   gp=gpar(col="black", fontsize=36, 
                                                           fontfamily="Times Roman")))
lyco.box.3 <- arrangeGrob(Lycopene.box.2, top = textGrob("(d)", x = unit(0, "npc")
                                                     , y   = unit(1, "npc"), just=c("left","top"),
                                                     gp=gpar(col="black", fontsize=36, 
                                                             fontfamily="Times Roman")))
phyto.box.3 <- arrangeGrob(Phytoene.box.2, top = textGrob("(e)", x = unit(0, "npc")
                                                       , y   = unit(1, "npc"), just=c("left","top"),
                                                       gp=gpar(col="black", fontsize=36, 
                                                               fontfamily="Times Roman")))

zeta.box.3 <- arrangeGrob(Zeta.box.2, top = textGrob("(f)", x = unit(0, "npc")
                                                     , y   = unit(1, "npc"), just=c("left","top"),
                                                     gp=gpar(col="black", fontsize=36, 
                                                             fontfamily="Times Roman")))

Total.box.3 <- arrangeGrob(Total.box.2, top = textGrob("(g)", x = unit(0, "npc")
                                                       , y   = unit(1, "npc"), just=c("left","top"),
                                                       gp=gpar(col="black", fontsize=36, 
                                                               fontfamily="Times Roman")))

AB.box.3 <- arrangeGrob(AB.ratio.box.2, top = textGrob("(h)", x = unit(0, "npc")
                                                 , y   = unit(1, "npc"), just=c("left","top"),
                                                 gp=gpar(col="black", fontsize=36, 
                                                         fontfamily="Times Roman")))

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

# The end :) 
