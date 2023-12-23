#1. assign methylation level to different methylation group rank

#methylation level	plot rank
#<5%	0
#5-15%	10
#15-25%	20
#25-35%	30
#35-45%	40
#45-55%	50
#55-65%	60
#65-75%	70
#75-85%	80
#85-95%	90
#>95% 100


#sample methylation level visualization with different methylation group rank
##line plot
library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)
library(reshape2)
library(ggthemes)
library(ggtext)

path <- "C:\\Users\\qjp\\Desktop\\Figures\\TissueCpGMethylationLevel\\"
data <- read.csv(paste0(path,"MethylationLevelOfDifferentRank.csv"),header = T)

colnames(data) <- str_sub(colnames(data),start = 2, end=-1)
colnames(data)[1] <- "Name"
colnames(data)[2] <- "Tissue"

data2 <- melt(data, id.vars = c("Name", "Tissue"), variable.name = "Level", value.name = "Percentage")

data2$Level <- factor(data2$Level, levels = c(0,10,20,30,40,50,60,70,80,90,100))

##
p1 <- ggplot(data2,aes(x=Level,y=Percentage,color = Tissue))+
  geom_line(aes(group = Name),size=0.5)+
  geom_point(aes(shape=Tissue),alpha = 0.8,size = 1)+
  scale_color_aaas()+xlab("Methylation Level")+
  ylab("CpG Sites Percentage")+
  scale_shape_manual(values=c(15,16,17,18))+
  theme_bw()+
  geom_richtext(aes(x= 9.7, y=30, label = "Brain"),  color = "#3B4992",show.legend = FALSE)+
  geom_richtext(aes(x= 8.2, y=20, label = "Liver"),   color = "#EE2200",show.legend = FALSE)+
  geom_richtext(aes(x=10.9,y=20, label = "Muscle"),  color = "#008B45",show.legend = FALSE)+ 
  geom_richtext(aes(x= 3, y=12, label = "Placenta"),color = "#631779",show.legend = FALSE)+
  theme(legend.position="none", axis.text=element_text(size=12), axis.title = element_text(size=12))

pdf("MethylationLevelOfDifferentRank.pdf",width=6,height=4)
p1      
dev.off()