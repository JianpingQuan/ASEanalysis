library(dplyr)

dat <- read.table("C:\\Users\\qjp\\Desktop\\POE_AGE_num_direction.txt", header = T)

dat$Group <- factor(dat$Group, levels=c("Pl_F70","Pl_F40","Mu_D168","Mu_D1","Mu_F70","Mu_F40","Li_D168","Li_D1","Li_F70","Li_F40","Br_D168","Br_D1","Br_F70","Br_F40"))


dat1 <- filter(dat, Class=="Maternal"|Class=="Duroc")
dat2 <- filter(dat, Class=="Paternal"|Class=="Lulai")
dat2$Number2 <- dat1$Number+dat2$Number

##plot
P <- ggplot()+
geom_col(data=dat2,aes(x=Number2,y=Group,fill=Class),position = "dodge")+
geom_col(data=dat1, aes(x=Number,y=Group,fill=Class),position = "dodge")+
geom_text(data = dat1, aes(x=Number,y=Group,fill=Class, label = Number),
position = position_dodge(0.9), hjust = 1.1) + 
geom_text(data = dat2, aes(x=Number2,y=Group,fill=Class, label = Number2),position = position_dodge(0.9), hjust = 1.1) + theme_bw()+scale_fill_manual(values = c("#42B540FF","#925E9FFF","#ED0000FF", "#00468BFF")) + scale_x_continuous(breaks=c(-45,-25,0,25,50),labels = c("45", "25", "0", "25", "50"))+ theme(legend.position = "top", legend.direction = "horizontal", legend.key.size = unit(0.75, "lines"), legend.box.spacing = unit(0.25, "lines")) + scale_x_continuous(position="top") + 
labs(x=NULL,y = NULL,title =NULL)

pdf("POE_AGE_number_dis2.pdf", height=7,width=5)
P
dev.off()