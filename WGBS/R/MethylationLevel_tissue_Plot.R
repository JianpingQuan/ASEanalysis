
#he global methylation level of each sample was calculated according to bismark methylation extractor module

library(ggplot2)
library(ggsci)
library(ggtext)

path <- "C:\\Users\\qjp\\Desktop\\Figures\\TissueCpGMethylationLevel\\"
dat <- read.table(paste0(path,"MethylationLevelOfTissues.txt"), header = T, row.names = 1)

ggplot(dat, aes(x=Tissue, y=MCpG*100, color=Tissue))+
  geom_boxplot()+
  scale_color_aaas()+
  theme_bw()+
  ylab(label = "Methylation level of CpG context")+
  theme(axis.ticks.x=element_blank(),axis.text.x = element_blank())+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 13))+
  geom_richtext(aes(x= 1, y=72, label = "Brain"),   color = "#3B4992",show.legend = FALSE)+
  geom_richtext(aes(x= 2, y=55, label = "Liver"),   color = "#EE2200",show.legend = FALSE)+
  geom_richtext(aes(x= 3, y=67, label = "Muscle"),  color = "#008B45",show.legend = FALSE)+ 
  geom_richtext(aes(x= 4, y=62, label = "Placenta"),color = "#631779",show.legend = FALSE)