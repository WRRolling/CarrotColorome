###### Part One: Prepare R environment ######
# 1.1 Clean R environment
rm(list = ls())
library('tidyverse')
library(readxl)
# 1.2 Load HPLC data and PI anems
WI18 <- read_csv("../Analysis_5_DescribePhenotypeData_WI18/WI18.HPLC.GWA.pheno.csv",
                 col_names=T)
CA19 <- read_csv("FilteredOrangeDat.csv", col_names = T)
which(WI18$Sample.Name %in% CA19$`Sample Name`)
WIPIs <- read_xlsx("../../Text/Supplemental Tablesv3.xlsx", sheet = 1)
CAPIs <- read_xlsx("../../Text/Supplemental Tablesv3.xlsx", sheet = 2)

# 1.3 Subset to only include PIs/Accessions found in both trials
WIPIs <- WIPIs[which(WIPIs$`Plot Number` %in% WI18$Sample.Name),]
CAPIs <- CAPIs[which(CAPIs$`Plot Number` %in% CA19$`Sample Name`),]
CA19[,13] <- CAPIs$Accession
colnames(CA19)[13] <- "PI"
CA19.1<- CA19[which(CA19$PI %in% WI18$ID),]
which(duplicated(CA19.1$PI) == T)
CA19.1 <- CA19.1[-234,]
CA19.1 <- CA19.1[-189,]
CA19.1 <- CA19.1[order(CA19.1$PI),] # same order
WI18.1 <- WI18[which(WI18$ID %in% CA19$PI),]             
WI18.1 <- WI18.1[order(WI18.1$ID),] # same order
CA19.1 <- CA19.1[which(CA19.1$PI %in% WI18.1$ID),]             
WI18.1 <- WI18.1[which(WI18.1$ID %in% CA19.1$PI),]             

# 1.4 Cor test between trials
cor.test(WI18.1$Alpha..ug.g., CA19.1$`Alpha (ug/g)`,  method = 'spearman')
cor.test(WI18.1$Beta..ug.g., CA19.1$`Beta (ug/g)`,  method = 'spearman')
cor.test(WI18.1$Lut..ug.g., CA19.1$`Lut (ug/g)`,  method = 'spearman')
cor.test(WI18.1$Lycopene..ug.g., CA19.1$`Lycopene (ug/g)`,  method = 'spearman')
cor.test(WI18.1$Phyto..ug.g., CA19.1$`Phyto (ug/g)`,  method = 'spearman')
cor.test(WI18.1$Zeta..ug.g., CA19.1$`Zeta (ug/g)`,  method = 'spearman')
cor.test(WI18.1$Total, CA19.1$`Total (ug/g)`,  method = 'spearman')
cor.test(WI18.1$AB.Ratio, CA19.1$AB.Ratio,  method = 'spearman')

# 1.5 Means and Standard Deviations
mean(WI18.1$Alpha..ug.g., na.rm=T)
sd(WI18.1$Alpha..ug.g., na.rm=T)
mean(WI18.1$Beta..ug.g., na.rm=T)
sd(WI18.1$Beta..ug.g., na.rm=T)
mean(WI18.1$Lut..ug.g., na.rm=T)
sd(WI18.1$Lut..ug.g., na.rm=T)
mean(WI18.1$Lycopene..ug.g., na.rm=T)
sd(WI18.1$Lycopene..ug.g., na.rm=T)
mean(WI18.1$Phyto..ug.g., na.rm=T)
sd(WI18.1$Phyto..ug.g., na.rm=T)
mean(WI18.1$Zeta..ug.g., na.rm=T)
sd(WI18.1$Zeta..ug.g., na.rm=T)
mean(WI18.1$Total, na.rm=T)
sd(WI18.1$Total, na.rm=T)
mean(WI18.1$AB.Ratio, na.rm=T)
sd(WI18.1$AB.Ratio, na.rm=T)


mean(CA19.1$`Alpha (ug/g)`, na.rm=T)
sd(CA19.1$`Alpha (ug/g)`)
mean(CA19.1$`Beta (ug/g)`, na.rm=T)
sd(CA19.1$`Beta (ug/g)`, na.rm=T)
mean(CA19.1$`Lut (ug/g)`, na.rm=T)
sd(CA19.1$`Lut (ug/g)`, na.rm=T)
mean(CA19.1$`Lycopene (ug/g)`, na.rm=T)
sd(CA19.1$`Lycopene (ug/g)`, na.rm=T)
mean(CA19.1$`Phyto (ug/g)`, na.rm=T)
sd(CA19.1$`Phyto (ug/g)`, na.rm=T)
mean(CA19.1$`Zeta (ug/g)`, na.rm=T)
sd(CA19.1$`Zeta (ug/g)`, na.rm=T)
mean(CA19.1$`Total (ug/g)`, na.rm=T)
sd(CA19.1$`Total (ug/g)`, na.rm=T)
mean(CA19.1$AB.Ratio, na.rm=T)
sd(CA19.1$AB.Ratio, na.rm=T)

#Write Output to file
write.csv(x=WI18.1, file="WI18Match.csv")
write.csv(x=CA19.1, file="CA19Match.csv")

####### Part Eleven: Mean Separation Between Years ######
# 11.1 Clean R environment
rm(list = ls())
# 11.2 Load a File For mean separation
dat <- read_csv(file="YearMeanSeparation.csv", col_names = T)
alpha.np <- kruskal.test(Alpha ~ Year, data = dat)
alpha.np$p.value 

beta.np <- kruskal.test(Beta ~ Year, data = dat)
beta.np$p.value 

lut.np <- kruskal.test(Lutein ~ Year, data = dat)
lut.np$p.value 

lyco.np <- kruskal.test(Lycopene ~ Year, data = dat)
lyco.np$p.value 

phyto.np <- kruskal.test(Phytoene ~ Year, data = dat)
phyto.np$p.value 

zeta.np <- kruskal.test(Zeta ~ Year, data = dat)
zeta.np$p.value

total.np <- kruskal.test(Total ~ Year, data = dat)
total.np$p.value

ratio.np <- kruskal.test(AB.Ratio ~ Year, data = dat)
ratio.np$p.value 


CA19.2 <- CA19.1[which(WI18.1$Interior == CA19.1$`Inner Color`),]
WI18.2 <- WI18.1[which(WI18.1$Interior == CA19.1$`Inner Color`),]
cor.test(CA19.2$`Alpha (ug/g)`, WI18.2$Alpha..ug.g.) #22.4
cor.test(CA19.2$`Beta (ug/g)`, WI18.2$Beta..ug.g.) # .202
cor.test(CA19.2$`Lut (ug/g)`, WI18.2$Lut..ug.g.) # 0.0
cor.test(CA19.2$`Lycopene (ug/g)`, WI18.2$Lycopene..ug.g.) # 0.14
cor.test(CA19.2$`Phyto (ug/g)`, WI18.2$Phyto..ug.g.) # 0.31
cor.test(CA19.2$`Zeta (ug/g)`, WI18.2$Zeta..ug.g.) # 0.007
cor.test(CA19.2$`Total (ug/g)`, WI18.2$Total) # 0.22
                  
                 