# Step Zero 
rm(list=ls()) 

# Step One: Load libraries
library('tidyverse')
library('factoextra')
library('rgl')


# Step Two: Load Data
dat <- read.csv("PCA.Pop.csv", head=T)

FigV0 <- ggplot(data=dat, 
                mapping = (aes(x=PC1, y=PC2)))+
                geom_point(size =2) + annotate(geom="text", x = 20, y=(-25), label="PC1 = 8.5% of Explained Variance", color="black") +
  annotate(geom="text", x = -25, y=(45), label="PC2 = 2.0% of Explained Variance", color="black") +
  theme_classic() 


