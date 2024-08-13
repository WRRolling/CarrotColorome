####### Part One:Prepare R Environment ####### 
# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Set working directory
setwd("../Desktop/") # Update with new directory
# (1.2) Load libraries
library('tidyverse')
library(gridExtra)
# (1.3) Load Data!
data <- read.csv("MarkerNumberData.csv", head=T)
####### Part Two:Create Figure to shower marker number plateau ####### 
# (2.0) Identify average accuracy based on marker number 
data_average <- data %>%
  group_by(Marker.Number, Dataset) %>%
  summarize(AvgCorrelation = mean(Correlation))
# (2.1) Identify time increase with more markers
data_average_time <- data %>%
  group_by(Marker.Number) %>%
  summarize(AvgTime = mean(Time), StandardDev = sd(Time)) 
# (2.2) Create Figure plotting all traits for all marker numbers and all models
Fig1 <- ggplot(data=data, aes(x=Marker.Number,y=Correlation, color=Dataset))+
  geom_point() +
  geom_line(data = data_average, aes(x = Marker.Number, y = AvgCorrelation, color = Dataset), linetype = "dashed") +  # Average line with dashed style
  ylab("Accuracy")+
  xlab("Number of Markers")+
  scale_color_manual(values = c("orange", "grey")) +
  theme(legend.position = "none")
# (2.3) Add annotations for the seconds required at approximately each marker number
Fig1.2 <- Fig1 + annotate("text", x=81, y=0.30, label="33 sec") + 
  annotate("text", x=323, y=0.37, label="115 sec") + 
  annotate("text", x=1291, y=0.385, label="170 sec") + 
  annotate("text", x=2582, y=0.39, label="331 sec") +
  annotate("text", x=5164, y=0.0, label="602 sec") + 
  annotate("text", x=10323, y=0.0, label="1135 sec") 

####### Part Two:Create Figure to Model Averages ####### 
# (3.0) SUbset to just 2500 markers
data.sub <- data[which(data$Marker.Number == 2582),]
# (3.1) Summarize the results of the marker number
data_average <- data.sub %>%
  group_by(Model, Dataset) %>%
  summarize(AvgCorrelation = mean(Correlation), StandardEv=sd(Correlation), Time=mean(Time))

Fig2 <- ggplot(data=data_average, aes(x=Model,y=AvgCorrelation, fill=Dataset))+
  geom_bar(stat='identity', position='dodge') +
  ylab("Accuracy")+
  xlab("BGLR Genomic Prediction Model") +
  scale_fill_manual(values = c("orange", "grey")) +
  geom_errorbar(aes(x=Model, ymin=AvgCorrelation-StandardEv, ymax=AvgCorrelation+StandardEv), 
                width = 0.3, position = position_dodge(0.9)) +
  #, width=0.4, colour="black", alpha=0.9, size=1.3) +
  annotate("text", x=.75, y=.255, label="470 sec") + 
  annotate("text", x=1.25, y=.22, label="403 sec") +
  annotate("text", x=1.75, y=.255, label="514 sec") + 
  annotate("text", x=2.25, y=.22, label="256 sec") +
  annotate("text", x=2.75, y=.255, label="483 sec") + 
  annotate("text", x=3.25, y=.22, label="222 sec") +
  annotate("text", x=3.75, y=.255, label="535 sec") + 
  annotate("text", x=4.25, y=.22, label="271 sec") +
  annotate("text", x=4.75, y=.255, label="433 sec") + 
  annotate("text", x=5.25, y=.22, label="201 sec") +
  annotate("text", x=5.75, y=.255, label="131 sec") + 
  annotate("text", x=6.25, y=.22, label="47 sec") 
  
####### Part Two:Create Figure to trait differences ####### 
# (3.0) SUbset to just 2500 markers # Need new data here. 

data.sub.2 <- data.sub[which(data.sub$Model == "BRR"),]
data.sub.2[,7] <- c(2.44E-02,
2.04E-02,
1.52E-02,
4.93E-02,
4.06E-02,
2.68E-02,
3.51E-02,
5.15E-02,
4.13E-02,
0.02318781,
3.20E-02,
2.79E-02,
0.02276592,
1.51E-02)

colnames(data.sub.2)[7] <- "StDevCor"


Fig3 <- ggplot(data=data.sub.2, aes(x=Trait, y=Correlation, fill=Dataset))+
  geom_bar(stat='identity', position='dodge') +
  ylab("Accuracy")+
  xlab("Carotenoid Predicted") +
  scale_fill_manual(values = c("orange", "grey")) +
  geom_errorbar(aes(x=Trait, ymin=Correlation-StDevCor, ymax=Correlation+StDevCor), 
                width = 0.3, position = position_dodge(0.9)) 


#################
data.3 <- read.csv(file="NextDat2.csv", head=T)
rownames(data.3) <- data.3[,1]
data.4 <- data.3
data.4 <- data.4[,-1]
ggradar(data.3) 

data4.5 <- data.frame(
  Total = c(0.5,0.0),
  Lutein = c(0.5, 0.0),
  Alpha = c(0.5, 0.0),
  Beta = c(0.5, 0.0),
  Lycopene = c(0.5, 0.0),
  Phytoene = c(0.5, 0.0),
  Zeta=c(0.5, 0.0)
)
rownames(data4.5) <- c("Max", "Mix")
data.5 <- rbind(data4.5, data.4)



create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = NA, plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

create_beautiful_radarchart(
   data = data.5, caxislabels = c(0, 0.125, 0.25, 0.375, 0.5),
   vlcex=2, calcex = 1.5,
   color = c("Orange", "Grey", "Red", "Yellow","Black", "Blue")
  )
legend(
  x = "bottom", legend = rownames(data.5[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20, col =c("Orange", "Grey", "Red", "Yellow","Black", "Blue") ,
  text.col = "black", cex = 2, pt.cex = 2, x.intersp = 0.5, y.intersp=-0.25,text.width=1
  
)


data.6 <- read.csv(file="NextDat3.csv", head=T)
rownames(data.6) <- data.6[,1]
data.7 <- data.6
data.7 <- data.7[,-1]
ggradar(data.7) 

data7.5 <- data.frame(
  Total = c(0.5,0.0),
  Lutein = c(0.5, 0.0),
  Alpha = c(0.5, 0.0),
  Beta = c(0.5, 0.0),
  Lycopene = c(0.5, 0.0),
  Phytoene = c(0.5, 0.0),
  Zeta=c(0.5, 0.0)
)
rownames(data4.5) <- c("Max", "Mix")
data.8 <- rbind(data7.5, data.7)




create_beautiful_radarchart(
  data = data.8, caxislabels = c(0, 0.125, 0.25, 0.375, 0.5),
  vlcex=2, calcex = 1.5,
  color = c("Orange", "Grey", "Red", "Yellow","Black", "Blue")
)
legend(
  x = "bottom", legend = rownames(data.5[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20, col =c("Orange", "Grey", "Red", "Yellow","Black", "Blue") ,
  text.col = "black", cex = 2, pt.cex = 2, x.intersp = 0.5, y.intersp=-0.25,text.width=1
  
)
