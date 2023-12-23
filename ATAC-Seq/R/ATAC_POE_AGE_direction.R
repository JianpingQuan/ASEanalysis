##------------------generated data ------------------
cd /f/Manuscript/ATAC/figure/POE05112023/ATAC_POE/method2/FilterDistal
ls *nodistal.csv|while read file
do
paste <(echo $file) <(cat $file|tr ',' '\t'|awk 'NR>1'|awk '$6<=0.05'|awk '$2>0'|wc -l) <(cat $file|tr ',' '\t'|awk 'NR>1'|awk '$6<=0.05'|awk '$2<0'|wc -l)
done

cd /f/Manuscript/ATAC/figure/POE05112023/ATAC_AGE/method2
ls *nodistal.csv|while read file
do
paste <(echo $file) <(cat $file|tr ',' '\t'|awk 'NR>1'|awk '$2>0'|wc -l) <(cat $file|tr ',' '\t'|awk 'NR>1'|awk '$2<0'|wc -l)
done
##------------------------------


##Plot1
library(dplyr)
library(ggplot2)
setwd("F:\\Manuscript\\ATAC\\figure\\mg_POE_AGE_biasDirecitionNum\\")
dat <- read.table("ATAC_POE_AGE_num_direction.txt", header = T)

dat$Group <- factor(dat$Group, levels=c("Pl_F70","Pl_F40","Mu_D168","Mu_D1","Mu_F70","Mu_F40","Li_D168","Li_D1","Li_F70","Li_F40","Br_D168","Br_D1","Br_F70","Br_F40"))


dat1 <- filter(dat, Class=="Maternal"|Class=="Duroc")
dat2 <- filter(dat, Class=="Paternal"|Class=="Lulai")
dat2$Number <- 0 - dat2$Number


##Plot2
p=ggplot()+geom_col(data=dat1, aes(x=Number,y=Group,fill=Class),position = "dodge")+
geom_col(data=dat2,aes(x=Number,y=Group,fill=Class),position = "dodge")+ 
geom_text(data = dat1, aes(x=Number,y=Group,fill=Class, label = Number),
            position = position_dodge(0.9), hjust = -0.01) + 
geom_text(data = dat2, 
            aes(x=Number,y=Group,fill=Class, label = (0-Number)),
            position = position_dodge(0.9), hjust = 1.01) + theme_bw()+
			scale_fill_manual(values = c("#42B540FF","#925E9FFF","#ED0000FF", "#00468BFF"))+ theme(legend.position = "top", legend.direction = "horizontal", legend.key.size = unit(0.75, "lines"), legend.box.spacing = unit(0.25, "lines"))+labs(x=NULL,y=NULL)+scale_x_continuous(limits = c(-250,250),position="top",breaks=c(-250,-150,-50,0,50,150,250),labels = c("250","150","50","0","50","150","250"))
pdf("ATAC_POE_AGE_number_direction.pdf", height=7,width=6)
p
dev.off()
			