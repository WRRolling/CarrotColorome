########
# Step One: Load Tidyverse
library('tidyverse')
library('matrixStats')

# Step Two: Load plink distance data
Dist.dat <- read_table("CA19DistMat.txt", 
                     col_names = F)

# Step Three: Format distance data
Dist.dat[Dist.dat == 0 ] <- NA
Dist.dat <- as.data.frame(Dist.dat)
row.names(Dist.dat) <- Dist.dat[,1]
colnames(Dist.dat) <- Dist.dat[,1]
Dist.dat <- Dist.dat[,-1]
Dist.dat <- Dist.dat[-1,]
Dist.dat <- as.matrix(Dist.dat)

# Step Four: Calculate the median identity by by state value for this population
All.Means <- colMedians(Dist.dat, na.rm = T)
Value.of.All.Means <- median(All.Means)
# 53381.45 

# Step Five: Load Population File
Population <- read_csv("CA19pops.csv", col_names = T)
Population <- Population[order(Population$`Root Number`),] 

# Step Seven Find Accessions in P1 - Eastern Germplasm 
P1 <- which(Population$Population == "Western")

P1.dat <- Population[P1,]
P1.IDs <- which(row.names(Dist.dat) %in% P1.dat$`Root Number`)
P1.dist.dat <-Dist.dat[P1.IDs,P1.IDs]

# Step Eight: Calculate Median
P1.Means <- colMedians(P1.dist.dat, na.rm = T)
P1.Mean <- median(P1.Means)
# 49642.5

# Step Nine Find Accessions in P2 - Western Germplasm 
P2 <- which(Population$Population == "Eastern")
P2.dat <- Population[P2,]
P2.IDs <- which(row.names(Dist.dat) %in% P2.dat$`Root Number`)
P2.dist.dat <-Dist.dat[P2.IDs,P2.IDs]

# Step Eight: Calculate Median
P2.Means <- colMedians(P2.dist.dat, na.rm = T)
P2.Mean <- median(P2.Means)
# 45184.3

# Step Nine: Find All P1 & P2 individuals
P1.P2.IDs <- which(Population$Population %in% c("Eastern","Western"))
P1.P2.dat <- Population[P1.P2.IDs,]      
P1.P2 <- which(row.names(Dist.dat) %in% P1.P2.dat$`Root Number`)
P1.P2.dist.dat <- Dist.dat[P1.P2,P1.P2]


# Step Ten: If else for loop to calculate distances
  
P1.P2.means <- list()

for (i in 1:length(P1.P2.dist.dat)) {
  if (P1.P2.dat$Population[i] == "Western") {
    wheres <- which(P1.P2.dat$Population == "Eastern")
    temp.dat <- P1.P2.dist.dat[c(i,wheres),c(i,wheres)]
    P1.P2.means[i] <- median(temp.dat[,1], na.rm = T) 
  } else if (P1.P2.dat$Population[i] == "Eastern") {
      wheres <- which(P1.P2.dat$Population == "Western")
      temp.dat <- P1.P2.dist.dat[c(i,wheres),c(i,wheres)]
      P1.P2.means[i]  <- median(temp.dat[,1], na.rm = T)
  }
  }

# Step Eleven: Calculate median
P1.P2.mean <- median(as.vector(unlist(P1.P2.means))) 
#58792.8

# Step Twelve: Format to write to csv file 
One <- c("All", "Within P1", "Within P2", "Between P1 and P2")
Two <- list()
Two[1] <- Value.of.All.Means
Two[2] <- P1.Mean
Two[3] <- P2.Mean
Two[4] <- P1.P2.mean
Two <- as.vector(unlist(Two))

# Step Thirteen: Write to csv file 
Output <- rbind(One, Two)
write.csv(x=Output, file="CAdistmeans.csv")

########
# Repeat for WI18 data

# Step Two: Load plink distance data
Dist.dat <- read_csv("WI18Distmat.csv", 
                       col_names = F)

# Step Three: Format distance data
Dist.dat[Dist.dat == 0 ] <- NA
Dist.dat <- as.data.frame(Dist.dat)
row.names(Dist.dat) <- Dist.dat[,1]
colnames(Dist.dat) <- Dist.dat[,1]
Dist.dat <- Dist.dat[,-1]
Dist.dat <- Dist.dat[-1,]
Dist.dat <- as.matrix(Dist.dat)

# Step Four: Calculate the median identity by by state value for this population
All.Means <- colMedians(Dist.dat, na.rm = T)
Value.of.All.Means <- median(All.Means)
# 52532.95 

# Step Five: Load Population File
Population <- read_csv("WI18Pops.csv", col_names = T)
Population <- Population[order(Population$Individual),] 

# Step Seven Find Accessions in P1 - Eastern Germplasm 
P1 <- which(Population$Population == "Western")

P1.dat <- Population[P1,]
P1.IDs <- which(row.names(Dist.dat) %in% P1.dat$Individual)
P1.dist.dat <-Dist.dat[P1.IDs,P1.IDs]

# Step Eight: Calculate Median
P1.Means <- colMedians(P1.dist.dat, na.rm = T)
P1.Mean <- median(P1.Means)
# 51795.25

# Step Nine Find Accessions in P2 - Western Germplasm 
P2 <- which(Population$Population == "Eastern")
P2.dat <- Population[P2,]
P2.IDs <- which(row.names(Dist.dat) %in% P2.dat$Individual)
P2.dist.dat <-Dist.dat[P2.IDs,P2.IDs]

# Step Eight: Calculate Median
P2.Means <- colMedians(P2.dist.dat, na.rm = T)
P2.Mean <- median(P2.Means)
# 47796.5

# Step Nine: Find All P1 & P2 individuals
P1.P2.IDs <- which(Population$Population %in% c("Eastern","Western"))
P1.P2.dat <- Population[P1.P2.IDs,]      
P1.P2 <- which(row.names(Dist.dat) %in% P1.P2.dat$Individual)
P1.P2.dist.dat <- Dist.dat[P1.P2,P1.P2]


# Step Ten: If else for loop to calculate distances

P1.P2.means <- list()

for (i in 1:length(P1.P2.dist.dat)) {
  if (P1.P2.dat$Population[i] == "Western") {
    wheres <- which(P1.P2.dat$Population == "Eastern")
    temp.dat <- P1.P2.dist.dat[c(i,wheres),c(i,wheres)]
    P1.P2.means[i] <- median(temp.dat[,1], na.rm = T) 
  } else if (P1.P2.dat$Population[i] == "Eastern") {
    wheres <- which(P1.P2.dat$Population == "Western")
    temp.dat <- P1.P2.dist.dat[c(i,wheres),c(i,wheres)]
    P1.P2.means[i]  <- median(temp.dat[,1], na.rm = T)
  }
}

# Step Eleven: Calculate median
P1.P2.mean <- median(as.vector(unlist(P1.P2.means))) 
#61821.1

# Step Twelve: Format to write to csv file 
One <- c("All", "Within P1", "Within P2", "Between P1 and P2")
Two <- list()
Two[1] <- Value.of.All.Means
Two[2] <- P1.Mean
Two[3] <- P2.Mean
Two[4] <- P1.P2.mean
Two <- as.vector(unlist(Two))

# Step Thirteen: Write to csv file 
Output <- rbind(One, Two)
write.csv(x=Output, file="WIdistmeans.csv")
