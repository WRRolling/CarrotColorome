####### Part One:Prepare R Environment & Load Phenotypic Data ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Load packages required for analysis
Pckg.Lst <-c("vroom","readxl","tidyverse", "PerformanceAnalytics",
             "lme4","lmerTest", "ggpubr", "corrplot")
package.check <- lapply(
  Pckg.Lst,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})
# (1.2) Load DataSheet. Data collected by Kevin Coe & Shelby Ellison 
Raw.data <- read_csv("RawDataWI18.csv", col_names=T)
View(Raw.data) # loaded correctly
# (1.3) Subset to only include those samples with "Orange" coloring
Color.dat <- read_csv("ColorScores.csv", col_names = T)
SubRaw.data <- Raw.data[which(Raw.data$Sample.Name %in% Color.dat$`Number`),]
# (1.4) get a list of columns in the dataframe
Raw.data %>% 
  colnames()
# (1.4) Subset to only columns needed for analyses
Raw.data.1 <- SubRaw.data %>%
  select(c("Sample.Name","Standard", "ug/g dry Lutein", "ug/g dry Lycopene",
           "ug/g dry Zeta","ug/g dry Alpha","ug/g dry Beta", 
           "ug/g dry Phytoene"))
# (1.5) Identify samples w/o tech replicate
Tech.Replciated <- 
  Raw.data.1$`Sample.Name`[duplicated(Raw.data.1$`Sample.Name`, )]
No.Tech.Rep <- which(Raw.data.1$`Sample.Name` %in% Tech.Replciated == FALSE)
No.Tech.Rep.1 <- Raw.data.1[No.Tech.Rep,] 
write.csv(No.Tech.Rep.1, file="NoTechRep.csv", row.names = F)
# Not applicable. 
rm(No.Tech.Rep)
rm(No.Tech.Rep.1)
####### Part Two: Check for tech rep issues from outliers  ####### 
# (2.1) Store ggplot in Raw.Lutein, specify data frame and carotenoid of interest
Raw.Lutein <- ggplot(data=Raw.data.1, aes(Raw.data.1$`ug/g dry Lutein`)) + 
  geom_histogram() + # Specify a ggplot - histogram
  xlab("Lutein (μg/g)") + # Update x axis label
  ylab("Observations") + # Change y axis label
  theme_classic() # Make the plot simple 
# (2.2) Create a histogram of the estimated concentration of each carotenoid
Raw.Lycopene <- ggplot(data=Raw.data.1, aes(Raw.data.1$`ug/g dry Lycopene`)) + 
  geom_histogram() +
  xlab("Lycopene (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.3) Create a histogram of the estimated concentration of each carotenoid
Raw.Zeta <- ggplot(data=Raw.data.1, aes(Raw.data.1$`ug/g dry Zeta`)) + 
  geom_histogram() +
  xlab("Zeta (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.4) Create a histogram of the estimated concentration of each carotenoid
Raw.Alpha <- ggplot(data=Raw.data.1, aes(Raw.data.1$`ug/g dry Alpha`)) + 
  geom_histogram() +
  xlab("Alpha (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.5) Create a histogram of the estimated concentration of each carotenoid
Raw.Beta <- ggplot(data=Raw.data.1, aes(Raw.data.1$`ug/g dry Beta`)) + 
  geom_histogram() +
  xlab("Beta (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.6) Create a histogram of the estimated concentration of each carotenoid
Raw.Phyto <- ggplot(data=Raw.data.1, aes(Raw.data.1$`ug/g dry Phytoene`)) + 
  geom_histogram() +
  xlab("Phytoene (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.7) Smoosh Plots together for quick viewing
Raw.data.histograms <- ggarrange(Raw.Alpha,Raw.Beta,Raw.Lutein, 
                                 Raw.Lycopene, Raw.Phyto, Raw.Zeta,
                                 labels = c("A", "B", "C", "D", "E", "F"),
                                 ncol=3, nrow=2)
# (2.7) Clean Up the environment
rm(Raw.Alpha, Raw.Beta, Raw.Lutein, Raw.Lycopene, Raw.Phyto, Raw.Zeta) 
# (2.8) View Histogram  
png(filename="RawHPLCHistogram.png",
    height=720, width = 720)
Raw.data.histograms
dev.off()

# (2.9) Identify samples with high Beta
Questionable.Sample <- which(x=Raw.data.1$`ug/g dry Beta` > 3000)  
Check.These.Closely <-   Raw.data.1$Sample.Name[Questionable.Sample]
#80675
# (2.10) Identify samples with high Alpha
Questionable.Sample <- which(x=Raw.data.1$`ug/g dry Alpha` > 2000)  
Check.These.Closely[2:3] <- Raw.data.1$Sample.Name[Questionable.Sample]
# (2.11) Identify samples with high Phyotene
Questionable.Sample <- which(x=Raw.data.1$`ug/g dry Phytoene` > 200)  
Raw.data.1$Sample.Name[Questionable.Sample]
Check.These.Closely[4:5] <- Raw.data.1$Sample.Name[Questionable.Sample]
# (2.12) Identify samples with high Zeta
Questionable.Sample <- which(x=Raw.data.1$`ug/g dry Zeta` > 200)  
Raw.data.1$Sample.Name[Questionable.Sample]
Check.These.Closely[6:8] <- Raw.data.1$Sample.Name[Questionable.Sample]
# (2.13) Identify samples with high Lycopene
Questionable.Sample <- which(x=Raw.data.1$`ug/g dry Lycopene` > 100)  
Raw.data.1$Sample.Name[Questionable.Sample]
Check.These.Closely[9:18] <- Raw.data.1$Sample.Name[Questionable.Sample]
# (2.14) Check List
Check.These.Closely
  #80675, 80585 only one tech rep listed as outlier
# (2.15) Check 80675 for issue
Checking <- which(x=Raw.data.1$Sample.Name == 80675)  
Tst <- Raw.data.1[Checking,]
# No problema yet
# (2.16) Check 80675 for issue
Checking <- which(x=Raw.data.1$Sample.Name == 80585)  
Tst <- Raw.data.1[Checking,]
# Huge difference between tech reps... remove sample
Raw.data.2 <- Raw.data.1[-Checking,]

####### Part Three: Check for a batch issue ####### 
# Not all HPLC data was collected at the same time. 
# A standard was included to correct for "Batch" effects, 
# Technical replicates were completed on the same HPLC run 
# All the first technical replicates going first
# HPLC-machine issues are identifiable by through tech reps
# (3.0) Identify each batch of HPLC completed by Standard
# Each run has a correction based on the Standard 
# This can be used to identify unique HPLC runs. 
Raw.data.2$Standard <- format(Raw.data.2$Standard, nsmall =2)
Batches <- unique(Raw.data.2$Standard)
# Fifteen "batches' of HPLC completed! 
# (3.1) Create New Column to differentiate batches with "categorical variable"
Raw.Data.2 <- Raw.data.2 %>%
  mutate(Batches=Raw.data.2$Standard)
# (3.2) Create list replace Numeric Raw.Data$Standard 
HPLC.Runs <- c("Uno", "Dos", "Tres", "Cuatro", "Cinco",
               "Sies", "Siete", "Ocho", "Nueve", "Diez",
               "Once", "Doce", "Trece", "Catorce", 
               "Quince")
# (3.3) Make new column to store cateogrical batch data
Raw.Data.2$Batches <- as.character(Raw.Data.2$Batches)
# (3.4) Replace Numeric batch information with categories  
for (i in 1:length(Batches)){
  Named <- which(Raw.Data.2$Batches == Batches[i])
  Raw.Data.2$Batches[Named] <- HPLC.Runs[i]}
# (3.5) Check for differences in average, minumum, and maximum between batches. 
Batch.Check <- Raw.Data.2 %>%
  group_by(Batches) %>%
  summarize(Alpha.Avg = mean(`ug/g dry Alpha`, na.rm=TRUE), 
            Alpha.Min = min(`ug/g dry Alpha`, na.rm=TRUE),
            Alpha.Max = max(`ug/g dry Alpha`, na.rm=TRUE),
            Beta.Avg = mean(`ug/g dry Beta`, na.rm=TRUE), 
            Beta.Min = min(`ug/g dry Beta`, na.rm=TRUE),
            Beta.Max = max(`ug/g dry Beta`, na.rm=TRUE),
            Lut.Avg = mean(`ug/g dry Lutein`, na.rm=TRUE), 
            Lut.Min = min(`ug/g dry Lutein`, na.rm=TRUE),
            Lut.Max = max(`ug/g dry Lutein`, na.rm=TRUE),
            Lyco.Avg = mean(`ug/g dry Lycopene`, na.rm=TRUE), 
            Lyco.Min = min(`ug/g dry Lycopene`, na.rm=TRUE),
            Lyco.Max = max(`ug/g dry Lycopene`, na.rm=TRUE),
            Phyto.Avg = mean(`ug/g dry Phytoene`, na.rm=TRUE), 
            Phyto.Min = min(`ug/g dry Phytoene`, na.rm=TRUE),
            Phyto.Max = max(`ug/g dry Phytoene`, na.rm=TRUE),
            Zeta.Avg = mean(`ug/g dry Zeta`, na.rm=TRUE), 
            Zeta.Min = min(`ug/g dry Zeta`, na.rm=TRUE),
            Zeta.Max = max(`ug/g dry Zeta`, na.rm=TRUE)) 
# (3.6) Look at results
View(Batch.Check)
# There are differences in batchs, though nothing too problematic

# (3.7) Check distribution of samples!  
Alpha.Batch <- ggplot(data=Raw.Data.2,aes(x=Batches, 
                                          y=`ug/g dry Alpha`, fill=Batches)) +
  geom_boxplot() + # Define boxplot
  theme_classic() + # keep it simple
  theme(legend.position = "none") # no legend needed

Beta.Batch <- ggplot(data=Raw.Data.2,aes(x=Batches, 
                                         y=`ug/g dry Beta`, fill=Batches)) +
  geom_boxplot() +
  theme_classic() + 
  theme(legend.position = "none")

Lut.Batch <- ggplot(data=Raw.Data.2,aes(x=Batches, 
                                        y=`ug/g dry Lutein`, fill=Batches)) +
  geom_boxplot() +
  theme_classic() + 
  theme(legend.position = "none")

Lyco.Batch <- ggplot(data=Raw.Data.2,aes(x=Batches, 
                                         y=`ug/g dry Lycopene`, fill=Batches)) +
  geom_boxplot() +
  theme_classic() + 
  theme(legend.position = "none")

Phyt.Batch <- ggplot(data=Raw.Data.2,aes(x=Batches, 
                                         y=`ug/g dry Phytoene`, fill=Batches)) +
  geom_boxplot() +
  theme_classic() + 
  theme(legend.position = "none")

Zeta.Batch <- ggplot(data=Raw.Data.2,aes(x=Batches, 
                                         y=`ug/g dry Zeta`, fill=Batches)) +
  geom_boxplot() +
  theme_classic() + 
  theme(legend.position = "none")

# (3.8) Smoosh Plots together for quick viewing
Batch.Check.Boxes <- ggarrange(Alpha.Batch,Beta.Batch,Lut.Batch, 
                               Lyco.Batch, Phyt.Batch, Zeta.Batch,
                               ncol=1, nrow=6)
# (3.9) Clean Up the environment
rm(Alpha.Batch,Beta.Batch,Lutein.Batch,Lycopene.Batch, Phyto.Batch, Zeta.Batch) 
# (3.10) Write Plot to file 
png(filename="BatchPhenoDis.png",
    height=1200, width = 600)
Batch.Check.Boxes  
dev.off()

####### Part Four: Calculate Average Concentration ########
Averaged.Data <- Raw.Data.2 %>%
  group_by(`Sample.Name`) %>%
  dplyr::summarize(`Alpha (ug/g)` = mean(`ug/g dry Alpha`, na.rm=TRUE), 
                   Alpha.sd = sd(`ug/g dry Alpha`, na.rm=TRUE),
                   `Beta (ug/g)` = mean(`ug/g dry Beta`, na.rm=TRUE), 
                   Beta.sd = sd(`ug/g dry Beta`, na.rm=TRUE),
                   `Lut (ug/g)` = mean(`ug/g dry Lutein`, na.rm=TRUE), 
                   Lut.sd = sd(`ug/g dry Lutein`, na.rm=TRUE),
                   `Lycopene (ug/g)` = mean(`ug/g dry Lycopene`, na.rm=TRUE), 
                   Lycopene.sd = sd(`ug/g dry Lutein`, na.rm=TRUE),
                   `Phyto (ug/g)` = mean(`ug/g dry Phytoene`, na.rm=TRUE), 
                   Phyto.sd = sd(`ug/g dry Phytoene`, na.rm=TRUE),
                   `Zeta (ug/g)` = mean(`ug/g dry Zeta`, na.rm=TRUE), 
                   Zeta.sd = sd(`ug/g dry Zeta`, na.rm=TRUE))

# (4.2) Add Batch Information to compare the variation between tech reps. 
for (i in 1:length(Averaged.Data$`Sample.Name`)){
  where <- match(Averaged.Data$`Sample.Name`[i], Raw.Data.2$`Sample.Name`)
  Averaged.Data[i,14] <- Raw.Data.2$Batches[where]
}
colnames(Averaged.Data)[14] <- "Batch.Eff"  

####### Part Five: Identify technical replicate issues ####### 
# 5.1 Calculated CV for tech replicates. 
# 5.1 Calculated CV for tech replicates.   
TechRepVar <- Raw.Data.2 %>%
  group_by(`Sample.Name`) %>%
  dplyr::summarize(`AlphaVar` = sd(`ug/g dry Alpha`, na.rm=TRUE)/
                     mean(`ug/g dry Alpha`, na.rm=TRUE), 
                   BetaVar = sd(`ug/g dry Beta`, na.rm=TRUE) /
                     mean(`ug/g dry Beta`, na.rm=TRUE), 
                   LutVar = sd(`ug/g dry Lutein`, na.rm=TRUE)/
                     mean(`ug/g dry Lutein`, na.rm=TRUE),
                   LycoVar = sd(`ug/g dry Lutein`, na.rm=TRUE)/ 
                     mean(`ug/g dry Lycopene`, na.rm=TRUE), 
                   PhytoVar = sd(`ug/g dry Phytoene`, na.rm=TRUE)/
                     mean(`ug/g dry Phytoene`, na.rm=TRUE), 
                   ZetaVar =sd(`ug/g dry Zeta`, na.rm=TRUE) /          
                     mean(`ug/g dry Zeta`, na.rm=TRUE))

# Notes: Threshold for CV = 0.1 (Coefficient of Variance)
# Either Alpha/Beta fail CV test - Remove Sample
# Other compounds remove just remove for that measurement
# Calculate Total from Alpha, Beta, Zeta & Phytoene
# Any failure of CV for the above, remove sample
# Lutein = Low confidence - use for ratio or for GWA w/ Yellow alone
# Lycopene = Low confidence 
# 5.2 OK Identify all samples that did not pass CV     
Alpha.Failed <- TechRepVar$`Sample.Name`[which(TechRepVar$AlphaVar > 0.1)]
Beta.Failed <-  TechRepVar$`Sample.Name`[which(TechRepVar$BetaVar > 0.1)] 
Phyto.Failed <-  TechRepVar$`Sample.Name`[which(TechRepVar$PhytoVar > 0.1)] 
Zeta.Failed <-   TechRepVar$`Sample.Name`[which(TechRepVar$ZetaVar > 0.1)]
Lyco.Failed <- TechRepVar$`Sample.Name`[which(TechRepVar$LycoVar > 0.1)]
Lut.Failed <- TechRepVar$`Sample.Name`[which(TechRepVar$LutVar > 0.1)] 


####### Part Six: Filter & Create Data Frame####### 
# 6.1 Set NAs in Averaged Data  
Output.data <- Averaged.Data %>%
  select(c(`Sample.Name`, `Alpha (ug/g)`, `Beta (ug/g)`,
           `Lut (ug/g)`,`Lycopene (ug/g)`,`Phyto (ug/g)`,
           `Zeta (ug/g)`))
# 6.2 Remove those samples that have large issues with Alpha/Beta
Output.data$`Alpha (ug/g)`[which(Output.data$`Sample.Name` %in% Alpha.Failed)] <- NA
Output.data$`Beta (ug/g)`[which(Output.data$`Sample.Name` %in% Beta.Failed)] <-NA
Output.data$`Lut (ug/g)`[which(Output.data$`Sample.Name` %in% Lut.Failed)] <- NA
Output.data$`Lycopene (ug/g)`[which(Output.data$`Sample.Name` %in% Lyco.Failed)] <- NA
Output.data$`Phyto (ug/g)`[which(Output.data$`Sample.Name`%in% Phyto.Failed)] <- NA
Output.data$`Zeta (ug/g)`[which(Output.data$`Sample.Name` %in% Zeta.Failed)]<- NA


# 6.3 Calculate a Total of alpha, beta, phytoene, and zeta
Output.data[,8] <- Output.data$`Alpha (ug/g)` + Output.data$`Beta (ug/g)`+
  Output.data$`Phyto (ug/g)`+ Output.data$`Zeta (ug/g)`
colnames(Output.data)[8] <- "Total (ug/g)"

# 6.4 Calculate a ratio of alpha & Beta to Lutein
Output.data[,9] <- (Output.data$`Alpha (ug/g)` + Output.data$`Beta (ug/g)`)/
  Output.data$`Lut (ug/g)`
colnames(Output.data)[9] <- "AB.Ratio"

# 6.5 Find rows that have only missing data
Missingdat <- rowSums(is.na(Output.data))
HighMissingValues <- Output.data[which(Missingdat > 3),]

# 6.6 Remove Samples with a high amount of missing data
Output.data.1 <- Output.data[-which(Missingdat >3),]
write.csv(x=Output.data.1, file="HPLC_dataWI18.csv", row.names = F)



