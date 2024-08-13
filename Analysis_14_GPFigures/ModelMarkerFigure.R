# Prepare R environment: Load Data and Library
library(tidyverse)
#gp.dat <- read.csv("CAModelMarkerData.csv", head=T)
gp.dat <- read.csv("CAModelMarkerData.csv", head=T)

gp.dat$Marker.Number <- as.numeric(gp.dat$Marker.Number)
#Separate Figure By Traits
Traits <- unique(gp.dat$Trait)

# First Trait
Traits[1]
#Alpha
al.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[1]),]

alpha <- ggplot(data=al.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="α-carotene", x="Marker Number")+
  guides(shape="none")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))


# Second Trait
Traits[2]
#Beta
beta.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[2]),]

beta <- ggplot(data=beta.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="β-carotene",  x="Marker Number")+
  guides(shape="none")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))

# Third Trait
Traits[3]
#Lutein
lut.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[3]),]

lut <- ggplot(data=lut.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="Lutein",  x="Marker Number")+
  guides(shape="none")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))


# Fourth Trait
Traits[4]
#Lycopene
lyco.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[4]),]

lyco <- ggplot(data=lyco.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="Lycopene",  x="Marker Number")+
  guides(shape="none")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))



# Fifth Trait
Traits[5]
#phytoene
phyto.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[5]),]

phyto <- ggplot(data=phyto.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="Phytoene", x="Marker Number")+
  guides(shape="none")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))
# Sixth Trait
Traits[6]
#Total
Total.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[6]),]

Total <- ggplot(data=Total.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="Total Carotenoids", x="Marker Number")+
  guides(shape="none")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))

# Seventh Trait
Traits[7]
#Total
Zeta.gp.dat <-gp.dat[which(gp.dat$Trait == Traits[7]),]

Zeta <- ggplot(data=Zeta.gp.dat, aes(x=Marker.Number, y=Correlation, group=Model))+
  geom_point(size=4,aes(shape=Model)) +
  labs(title="ζ-carotene", x="Marker Number")+
  theme(text=element_text(size=21),plot.title = element_text(hjust = 0.5))

library(ggpubr)
png(filename="CAModelMarker.png", height = 1800, width=1500)
ggarrange(alpha, beta, lut, lyco, phyto, Total, Zeta, 
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 2, nrow = 4)
dev.off()
