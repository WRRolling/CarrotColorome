###################################
# Load tidyverse to organize data
library(tidyverse)

##################################
# determine if there is a  level of correlation in the data!
# load random results
ran1 <- read.csv("../Analysis_14/SecondRound/CA19/Mark81/RandomAccuracy.csv", 
                 head=T)
ran2 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark161/RandomAccuracy.csv", 
                 head=T)
ran3 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark323/RandomAccuracy.csv", 
               head=T)
ran4 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark646/RandomAccuracy.csv", 
         head=T)
ran5 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark1291/RandomAccuracy.csv", 
         head=T)
ran6 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark2582/RandomAccuracy.csv", 
                 head=T)
ran7 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark5164/RandomAccuracy.csv", 
                 head=T)
ran8 <- read.csv(file="../Analysis_14/SecondRound/CA19/Mark10328/RandomAccuracy.csv", 
                 head=T)

# combine random results
random.dat <- rbind(ran1,ran2,ran3,ran4,ran5,ran6,ran7,ran8)
colnames(random.dat) <- c("Trait", "Correlation")
# calculate averages
Baseline <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 
# write output
write.csv(x=Baseline, file="RandomAccuracy", row.names = F)

###################################
# determine if there is population alone can provide accurate predictions!
# load population alone results
pop1 <- read.csv("../Analysis_14/SecondRound/CA19/Mark81/PopulationAlone.csv", 
                 head=T)
#pop2 <- read.csv("../Analysis_14/SecondRound/CA19/Mark161/PopulationAlone.csv", 
                # head=T) # Forgot to run in Mark161
#pop3 <- read.csv("../Analysis_14/SecondRound/CA19/Mark323/PopulationAlone.csv", 
                 # head=T) for to run in Mark323
pop4 <- read.csv("../Analysis_14/SecondRound/CA19/Mark646/PopulationAlone.csv", 
                 head=T)
pop5 <- read.csv("../Analysis_14/SecondRound/CA19/Mark1291/PopulationAlone.csv", 
                 head=T)
pop6 <- read.csv("../Analysis_14/SecondRound/CA19/Mark2582/PopulationAlone.csv", 
                 head=T)
pop7 <- read.csv("../Analysis_14/SecondRound/CA19/Mark5164/PopulationAlone.csv", 
                 head=T)
pop8 <- read.csv("../Analysis_14/SecondRound/CA19/Mark10328/PopulationAlone.csv", 
                 head=T)
# combine population alone
random.dat <- rbind(pop1,pop4,pop5,pop6,pop7,pop8)
colnames(random.dat) <- c("Trait", "Correlation", 
                          "pvalue",
                          "H2", 
                          "Top25",
                          "Bot25", 
                          "Time")
# Summarize
Baseline.2 <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 
# write output
write.csv(x=Baseline.2, file="PopulationAccuracy.csv", row.names = F)
###################################
# determine how accurate you would be with core color alone!
dat1 <- read.csv("../Analysis_14/SecondRound/CA19/Mark81/CoreAlone.csv", 
                 head=T)
#dat2 <- read.csv("../Analysis_14/SecondRound/CA19/Mark161/CoreAlone.csv", 
                # head=T) # forgot 
#dat3 <- read.csv("../Analysis_14/SecondRound/CA19/Mark323/CoreAlone.csv", 
               #  head=T) # forgot
dat4 <- read.csv("../Analysis_14/SecondRound/CA19/Mark646/CoreAlone.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/CA19/Mark1291/CoreAlone.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/CA19/Mark2582/CoreAlone.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/CA19/Mark5164/CoreAlone.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/CA19/Mark10328/CoreAlone.csv", 
                 head=T)
# combine population alone
random.dat <- rbind(dat1,dat4,dat5,dat6,dat7,dat8)
colnames(random.dat) <- c("Trait", "Correlation", 
                          "pvalue",
                          "H2", 
                          "Top25",
                          "Bot25", 
                          "Time")
# Summarize
Baseline.3 <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 

# write output
write.csv(x=Baseline.3, file="CoreAccuracy.csv", row.names = F)

###################################
# How well will kinship alone perform? 
dat1 <- read.csv("../Analysis_14/SecondRound/CA19/Mark81/KinshipAlone.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/CA19/Mark161/KinshipAlone.csv", 
               head=T)
dat3 <- read.csv("../Analysis_14/SecondRound/CA19/Mark323/KinshipAlone.csv", 
                 head=T)
dat4 <- read.csv("../Analysis_14/SecondRound/CA19/Mark646/KinshipAlone.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/CA19/Mark1291/KinshipAlone.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/CA19/Mark2582/KinshipAlone.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/CA19/Mark5164/KinshipAlone.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/CA19/Mark10328/KinshipAlone.csv", 
                 head=T)

# Rename Traits to match
Total <- which(str_detect(string=dat1$X, "Total")==T)
Lutein <-  which(str_detect(string=dat1$X, "Lutein")==T)
Alpha <-  which(str_detect(string=dat1$X, "Alpha")==T)
Beta <-  which(str_detect(string=dat1$X, "Beta")==T)
Lycopene <-  which(str_detect(string=dat1$X, "Lycopene")==T)
Phytoene <-  which(str_detect(string=dat1$X, "Phytoene")==T)
Zeta <-  which(str_detect(string=dat1$X, "Zeta")==T)

rename <- function(input){
  modified_input <- input
  modified_input$X[Total] <- "Total"
  modified_input$X[Lutein] <- "Lutein"
  modified_input$X[Alpha] <- "Alpha"
  modified_input$X[Beta] <- "Beta"
  modified_input$X[Lycopene] <- "Lycopene"
  modified_input$X[Phytoene] <- "Phytoene"
  modified_input$X[Zeta] <- "Zeta"
  return(modified_input)
}


# Rename Traits to match need to update reference for 8th only 8 iterations
#Total <- which(str_detect(string=dat8$X, "Total")==T)
#Lutein <-  which(str_detect(string=dat8$X, "Lutein")==T)
#Alpha <-  which(str_detect(string=dat8$X, "Alpha")==T)
##Beta <-  which(str_detect(string=dat8$X, "Beta")==T)
#Lycopene <-  which(str_detect(string=dat8$X, "Lycopene")==T)
#Phytoene <-  which(str_detect(string=dat8$X, "Phytoene")==T)
#Zeta <-  which(str_detect(string=dat8$X, "Zeta")==T)

dat1 <- rename(dat1)
dat2 <- rename(dat2)
dat3 <- rename(dat3)
dat4<- rename(dat4)
dat5 <- rename(dat5)
dat6 <- rename(dat6)
dat7 <- rename(dat7)
dat8 <- rename(dat8)


# check results for each trait
dat1.1 <- dat1 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat1.1[,14] <- 81

dat2.1 <- dat2 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat2.1[,14] <- 161

dat3.1 <- dat3 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat3.1[,14] <- 323

dat4.1 <- dat4 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat4.1[,14] <- 646

dat5.1 <- dat5 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat5.1[,14] <- 1291


dat6.1 <- dat6 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat6.1[,14] <- 2582


dat7.1 <- dat7 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat7.1[,14] <- 5164

dat8.1 <- dat8 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat8.1[,14] <- 10323


combinedKin <- rbind(dat1.1, dat2.1,dat3.1,dat4.1,dat5.1,
      dat6.1,dat7.1,dat8.1)

write.csv(x=combinedKin, file="KinshipAlone.csv", row.names = F)

###################################
# Test how well a prediction models + SNPs do alone. 
dat1 <- read.csv("../Analysis_14/SecondRound/CA19/Mark81/ModSummary.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/CA19/Mark161/ModSummary.csv", 
                 head=T)
dat3 <- read.csv("../Analysis_14/SecondRound/CA19/Mark323/ModSummary.csv", 
                 head=T)
dat4 <- read.csv("../Analysis_14/SecondRound/CA19/Mark646/ModSummary.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/CA19/Mark1291/ModSummary.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/CA19/Mark2582/ModSummary.csv", 
                head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/CA19/Mark5164/ModSummary.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/CA19/Mark10328/ModSummary.csv", 
                 head=T)
dat1.1 <- as.data.frame(t(dat1))
dat2.1 <- as.data.frame(t(dat2))
dat3.1 <- as.data.frame(t(dat3))
dat4.1 <- as.data.frame(t(dat4))
dat5.1 <- as.data.frame(t(dat5))
dat6.1 <- as.data.frame(t(dat6))
dat7.1 <- as.data.frame(t(dat7))
dat8.1 <- as.data.frame(t(dat8))

format.mod.data <- function(input){
  modified.input <- input
  colnames(modified.input) <- modified.input[1,]
  modified.input <- modified.input[-1,]
  modified.input[,13] <- row.names(modified.input)
  modified.input <- separate_wider_delim(modified.input , cols = V13, delim = ".", 
                                 names = c("Model","Waste","Trait"))
  modified.input <- modified.input[,-14]
  return(modified.input)
}

dat1.2 <- format.mod.data(dat1.1)
dat2.2 <- format.mod.data(dat2.1)
dat3.2 <- format.mod.data(dat3.1)
dat4.2 <- format.mod.data(dat4.1)
dat5.2 <- format.mod.data(dat5.1)
dat6.2 <- format.mod.data(dat6.1)
dat7.2 <- format.mod.data(dat7.1)
dat8.2 <- format.mod.data(dat8.1)

dat1.2[,15] <-  81
dat2.2[,15] <-  161
dat3.2[,15] <-  323
dat4.2[,15] <-  646
dat5.2[,15] <-  1291
dat6.2[,15] <-  2582
dat7.2[,15] <-  5164
dat8.2[,15] <-  10323

combinedMods <- rbind(dat1.2, dat2.2,dat3.2,dat4.2, dat5.2,
                      dat6.2,dat7.2,dat8.2)
write.csv(x=combinedMods, file="ModelOutputs.csv", row.names = F)

###################################
# Consider how adding core, population, kinship to models!
dat1 <- read.csv("../Analysis_14/SecondRound/CA19/Mark81/ModCombinedSummary.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/CA19/Mark161/ModCombinedSummary.csv", 
                 head=T)
dat3 <- read.csv("../Analysis_14/SecondRound/CA19/Mark323/ModCombinedSummary.csv", 
                 head=T)
dat4 <- read.csv("../Analysis_14/SecondRound/CA19/Mark646/ModCombinedSummary.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/CA19/Mark1291/ModCombinedSummary.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/CA19/Mark2582/ModCombinedSummary.csv", 
                head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/CA19/Mark5164/ModCombinedSummary.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/CA19/Mark10328/ModCombinedSummary.csv", 
                 head=T)
dat1.1 <- as.data.frame(t(dat1))
dat2.1 <- as.data.frame(t(dat2))
dat3.1 <- as.data.frame(t(dat3))
dat4.1 <- as.data.frame(t(dat4))
dat5.1 <- as.data.frame(t(dat5))
dat6.1 <- as.data.frame(t(dat6))
dat7.1 <- as.data.frame(t(dat7))
dat8.1 <- as.data.frame(t(dat8))

dat1.2 <- format.mod.data(dat1.1)
dat2.2 <- format.mod.data(dat2.1)
dat3.2 <- format.mod.data(dat3.1)
dat4.2 <- format.mod.data(dat4.1)
dat5.2 <- format.mod.data(dat5.1)
dat6.2 <- format.mod.data(dat6.1)
dat7.2 <- format.mod.data(dat7.1)
dat8.2 <- format.mod.data(dat8.1)
dat1.2[,15] <-  81
dat2.2[,15] <-  161
dat3.2[,15] <-  323
dat4.2[,15] <-  646
dat5.2[,15] <-  1291
dat6.2[,15] <-  2582
dat7.2[,15] <-  5164
dat8.2[,15] <-  10328

combinedMods <- rbind(dat1.2, dat2.2,dat3.2,dat4.2, dat5.2, dat6.2,dat7.2,dat8.2)
write.csv(x=combinedMods, file="FullModelOutputs.csv", row.names = F)

#########################################
# Wisconsin Dataset!
########################################
ran1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/RandomAccuracy.csv", 
                 head=T)
#ran2 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark161/RandomAccuracy.csv", 
                 #head=T) Forgot
ran3 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark323/RandomAccuracy.csv", 
                 head=T)
#ran4 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark646/RandomAccuracy.csv", 
               #  head=T)
ran5 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark1291/RandomAccuracy.csv", 
                 head=T)
ran6 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark2582/RandomAccuracy.csv", 
                 head=T)
ran7 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark5164/RandomAccuracy.csv", 
                 head=T)
ran8 <- read.csv(file="../Analysis_14/SecondRound/WI18/Mark10328/RandomAccuracy.csv", 
                 head=T)

# combine random results
random.dat <- rbind(ran1,ran3,ran5,ran6,ran7,ran8)
colnames(random.dat) <- c("Trait", "Correlation")
# calculate averages
Baseline <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 
# write output
write.csv(x=Baseline, file="WIoutput/RandomAccuracy", row.names = F)
# determine if there is population alone can provide accurate predictions!
# load population alone results

##################
pop1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/PopulationAlone.csv", 
                 head=T)
#pop2 <- read.csv("../Analysis_14/SecondRound/WI18/Mark161/PopulationAlone.csv", 
 #head=T) # Forgot to run in Mark161
pop3 <- read.csv("../Analysis_14/SecondRound/WI18/Mark323/PopulationAlone.csv", 
 head=T) 
pop4 <- read.csv("../Analysis_14/SecondRound/WI18/Mark646/PopulationAlone.csv", 
                 head=T)
pop5 <- read.csv("../Analysis_14/SecondRound/WI18/Mark1291/PopulationAlone.csv", 
                 head=T)
pop6 <- read.csv("../Analysis_14/SecondRound/WI18/Mark2582/PopulationAlone.csv", 
                 head=T)
pop7 <- read.csv("../Analysis_14/SecondRound/WI18/Mark5164/PopulationAlone.csv", 
                 head=T)
pop8 <- read.csv("../Analysis_14/SecondRound/WI18/Mark10328/PopulationAlone.csv", 
                 head=T)
# combine population alone
random.dat <- rbind(pop1,pop3, pop4,pop5,pop6,pop7,pop8)
colnames(random.dat) <- c("Trait", "Correlation", 
                          "pvalue",
                          "H2", 
                          "Top25",
                          "Bot25", 
                          "Time")
# Summarize
Baseline.2 <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 
# write output
write.csv(x=Baseline.2, file="WIoutput/PopulationAccuracy.csv", row.names = F)
###################################
# determine how accurate you would be with core color alone!
dat1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/CoreAlone.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/WI18/Mark161/CoreAlone.csv", 
 head=T) # forgot 
dat3 <- read.csv("../Analysis_14/SecondRound/WI18/Mark323/CoreAlone.csv", 
  head=T) # forgot
dat4 <- read.csv("../Analysis_14/SecondRound/WI18/Mark646/CoreAlone.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/WI18/Mark1291/CoreAlone.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/WI18/Mark2582/CoreAlone.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/WI18/Mark5164/CoreAlone.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/WI18/Mark10328/CoreAlone.csv", 
                 head=T)
# combine population alone
random.dat <- rbind(dat1,dat4,dat5,dat6,dat7,dat8)
colnames(random.dat) <- c("Trait", "Correlation", 
                          "pvalue",
                          "H2", 
                          "Top25",
                          "Bot25", 
                          "Time")
# Summarize
Baseline.3 <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 

# write output
write.csv(x=Baseline.3, file="CoreAccuracy.csv", row.names = F)

# determine how accurate you would be with core color alone!
dat1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/CorePop.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/WI18/Mark161/CorePop.csv", 
                 head=T) # forgot 
dat3 <- read.csv("../Analysis_14/SecondRound/WI18/Mark323/CorePop.csv", 
                 head=T) # forgot
dat4 <- read.csv("../Analysis_14/SecondRound/WI18/Mark646/CorePop.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/WI18/Mark1291/CorePop.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/WI18/Mark2582/CorePop.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/WI18/Mark5164/CorePop.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/WI18/Mark10328/CorePop.csv", 
                 head=T)
# combine population alone
random.dat <- rbind(dat1,dat3,dat4,dat5,dat6,dat7,dat8)
colnames(random.dat) <- c("Trait", "Correlation", 
                          "pvalue",
                          "H2", 
                          "Top25",
                          "Bot25", 
                          "Time")
# Summarize
Baseline.3 <- random.dat %>%
  group_by(Trait) %>%
  summarise_at(vars(Correlation), list(Average = mean, Stdev = sd)) 

# write output
write.csv(x=Baseline.3, file="CorePop.csv", row.names = F)


###################################
# How well will kinship alone perform? 
dat1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/KinshipAlone.csv", 
                 head=T)
#dat2 <- read.csv("../Analysis_14/SecondRound/WI18/Mark161/KinshipAlone.csv", 
                 #head=T) #forgot
dat3 <- read.csv("../Analysis_14/SecondRound/WI18/Mark323/KinshipAlone.csv", 
                 head=T)
dat4 <- read.csv("../Analysis_14/SecondRound/WI18/Mark646/KinshipAlone.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/WI18/Mark1291/KinshipAlone.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/WI18/Mark2582/KinshipAlone.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/WI18/Mark5164/KinshipAlone.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/WI18/Mark10328/KinshipAlone.csv", 
                 head=T)

# Rename Traits to match
Total <- which(str_detect(string=dat1$X, "Total")==T)
Lutein <-  which(str_detect(string=dat1$X, "Lutein")==T)
Alpha <-  which(str_detect(string=dat1$X, "Alpha")==T)
Beta <-  which(str_detect(string=dat1$X, "Beta")==T)
Lycopene <-  which(str_detect(string=dat1$X, "Lycopene")==T)
Phytoene <-  which(str_detect(string=dat1$X, "Phytoene")==T)
Zeta <-  which(str_detect(string=dat1$X, "Zeta")==T)

rename <- function(input){
  modified_input <- input
  modified_input$X[Total] <- "Total"
  modified_input$X[Lutein] <- "Lutein"
  modified_input$X[Alpha] <- "Alpha"
  modified_input$X[Beta] <- "Beta"
  modified_input$X[Lycopene] <- "Lycopene"
  modified_input$X[Phytoene] <- "Phytoene"
  modified_input$X[Zeta] <- "Zeta"
  return(modified_input)
}



dat1 <- rename(dat1)
#dat2 <- rename(dat2)
dat3 <- rename(dat3)
dat4<- rename(dat4)
dat5 <- rename(dat5)
dat6 <- rename(dat6)
dat7 <- rename(dat7)
# Rename Traits to match need to update reference for 8th only 8 iterations
Total <- which(str_detect(string=dat8$X, "Total")==T)
Lutein <-  which(str_detect(string=dat8$X, "Lutein")==T)
Alpha <-  which(str_detect(string=dat8$X, "Alpha")==T)
Beta <-  which(str_detect(string=dat8$X, "Beta")==T)
Lycopene <-  which(str_detect(string=dat8$X, "Lycopene")==T)
Phytoene <-  which(str_detect(string=dat8$X, "Phytoene")==T)
Zeta <-  which(str_detect(string=dat8$X, "Zeta")==T)

dat8 <- rename(dat8)


# check results for each trait
dat1.1 <- dat1 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat1.1[,14] <- 81

#dat2.1 <- dat2 %>%
#  group_by(X) %>%
#  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
#dat2.1[,14] <- 161

dat3.1 <- dat3 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat3.1[,14] <- 323

dat4.1 <- dat4 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat4.1[,14] <- 646

dat5.1 <- dat5 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat5.1[,14] <- 1291


dat6.1 <- dat6 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat6.1[,14] <- 2582


dat7.1 <- dat7 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat7.1[,14] <- 5164

dat8.1 <- dat8 %>%
  group_by(X) %>%
  summarise_at(vars(Cor, Pval, H2, Top25,Bot50, Time),
               list(Average = mean, Stdev = sd))
dat8.1[,14] <- 10323


combinedKin <- rbind(dat1.1, dat3.1,dat4.1,dat5.1,
                     dat6.1,dat7.1,dat8.1)

write.csv(x=combinedKin, file="WIoutput/KinshipAlone.csv", row.names = F)

###################################
# Test how well a prediction models + SNPs do alone. 
dat1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/ModSummary.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/WI18/Mark162/ModSummary.csv", 
                 head=T)
dat3 <- read.csv("../Analysis_14/SecondRound/WI18/Mark323/ModSummary.csv", 
                 head=T)
dat4 <- read.csv("../Analysis_14/SecondRound/WI18/Mark646/ModSummary.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/WI18/Mark1291/ModSummary.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/WI18/Mark2582/ModSummary.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/WI18/Mark5164/ModSummary.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/WI18/Mark10328/ModSummary.csv", 
                 head=T)
dat1.1 <- as.data.frame(t(dat1))
dat2.1 <- as.data.frame(t(dat2))
dat3.1 <- as.data.frame(t(dat3))
dat4.1 <- as.data.frame(t(dat4))
dat5.1 <- as.data.frame(t(dat5))
dat6.1 <- as.data.frame(t(dat6))
dat7.1 <- as.data.frame(t(dat7))
dat8.1 <- as.data.frame(t(dat8))

format.mod.data <- function(input){
  modified.input <- input
  colnames(modified.input) <- modified.input[1,]
  modified.input <- modified.input[-1,]
  modified.input[,13] <- row.names(modified.input)
  modified.input <- separate_wider_delim(modified.input , cols = V13, delim = ".", 
                                         names = c("Model","Waste","Trait"))
  modified.input <- modified.input[,-14]
  return(modified.input)
}

dat1.2 <- format.mod.data(dat1.1)
dat2.2 <- format.mod.data(dat2.1)
dat3.2 <- format.mod.data(dat3.1)
dat4.2 <- format.mod.data(dat4.1)
dat5.2 <- format.mod.data(dat5.1)
dat6.2 <- format.mod.data(dat6.1)
dat7.2 <- format.mod.data(dat7.1)
dat8.2 <- format.mod.data(dat8.1)

dat1.2[,15] <-  81
dat2.2[,15] <-  161
dat3.2[,15] <-  323
dat4.2[,15] <-  646
dat5.2[,15] <-  1291
dat6.2[,15] <-  2582
dat7.2[,15] <-  5164
dat8.2[,15] <-  10323

combinedMods <- rbind(dat1.2, dat2.2,dat3.2,dat4.2, dat5.2,
                      dat6.2,dat7.2,dat8.2)
write.csv(x=combinedMods, file="WIoutput/ModelOutputs.csv", row.names = F)
###################################
# Consider how adding core, population, kinship to models!
dat1 <- read.csv("../Analysis_14/SecondRound/WI18/Mark81/ModCombinedSummary.csv", 
                 head=T)
dat2 <- read.csv("../Analysis_14/SecondRound/WI18/Mark162/ModCombinedSummary.csv", 
                 head=T)
dat3 <- read.csv("../Analysis_14/SecondRound/WI18/Mark323/ModCombinedSummary.csv", 
                 head=T)
dat4 <- read.csv("../Analysis_14/SecondRound/WI18/Mark646/ModCombinedSummary.csv", 
                 head=T)
dat5 <- read.csv("../Analysis_14/SecondRound/WI18/Mark1291/ModCombinedSummary.csv", 
                 head=T)
dat6 <- read.csv("../Analysis_14/SecondRound/WI18/Mark2582/ModCombinedSummary.csv", 
                 head=T)
dat7 <- read.csv("../Analysis_14/SecondRound/WI18/Mark5164/ModCombinedSummary.csv", 
                 head=T)
dat8 <- read.csv("../Analysis_14/SecondRound/WI18/Mark10328/ModCombinedSummary.csv", 
                 head=T)
dat1.1 <- as.data.frame(t(dat1))
dat2.1 <- as.data.frame(t(dat2))
dat3.1 <- as.data.frame(t(dat3))
dat4.1 <- as.data.frame(t(dat4))
dat5.1 <- as.data.frame(t(dat5))
dat6.1 <- as.data.frame(t(dat6))
dat7.1 <- as.data.frame(t(dat7))
dat8.1 <- as.data.frame(t(dat8))

dat1.2 <- format.mod.data(dat1.1)
dat2.2 <- format.mod.data(dat2.1)
dat3.2 <- format.mod.data(dat3.1)
dat4.2 <- format.mod.data(dat4.1)
dat5.2 <- format.mod.data(dat5.1)
dat6.2 <- format.mod.data(dat6.1)
dat7.2 <- format.mod.data(dat7.1)
dat8.2 <- format.mod.data(dat8.1)
dat1.2[,15] <-  81
dat2.2[,15] <-  161
dat3.2[,15] <-  323
dat4.2[,15] <-  646
dat5.2[,15] <-  1291
dat6.2[,15] <-  2582
dat7.2[,15] <-  5164
dat8.2[,15] <-  10328

combinedMods <- rbind(dat1.2, dat2.2,dat3.2,dat4.2, dat5.2, dat6.2,dat7.2,dat8.2)
write.csv(x=combinedMods, file="WIoutput/FullModelOutputs.csv", row.names = F)