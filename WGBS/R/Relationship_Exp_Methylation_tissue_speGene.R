#relationship between methylation and expression of tissues-specific gene

library(dplyr)
library(ggplot2)
library(magrittr)
library(ggtext)
library(ggpubr)

dataExp <- read.csv("F:\\RNA\\RawTables\\OriginalRef\\GeneCountMatrix.tpmUnit.csv", header=T, row.names=1) %>% select(., contains("70"))

feature <- c("geneBody", "gene1stExon", "geneIntron", "genePromoter","geneUTR5","geneUTR3")

## choose which feature to analysis
i= 2

methdata <- paste0("C:\\Users\\qjp\\Desktop\\genomeFeature20220823\\data\\AllSample.", feature[i],".methyLevel.txt")

Tis_meth <- read.table(methdata, header = F,row.names = 1)

sampleName <- read.table("C:\\Users\\qjp\\Desktop\\genomeFeature20220823\\data\\sampleName.txt", header = F)

colnames(Tis_meth) <- sampleName$V1

Tis_spe <- read.csv("C:\\Users\\qjp\\Desktop\\RNA-Tissue_specific.genes.csv")

gene_com <- intersect(rownames(Tis_meth), rownames(dataExp))

gene_com_exp  <- filter(dataExp, rownames(dataExp) %in% gene_com) %>% arrange(.,rownames(.))

gene_com_meth <- filter(Tis_meth, rownames(Tis_meth) %in% gene_com) %>% arrange(.,rownames(.))

gene_com_exp_tis_mean <- data.frame(Br=rowMeans(select(gene_com_exp, contains("Br"))), Li=rowMeans(select(gene_com_exp, contains("Li"))), Mu=rowMeans(select(gene_com_exp, contains("Mu"))), Pl=rowMeans(select(gene_com_exp, contains("Pl"))))

gene_com_meth_tis_mean <- data.frame(Br=rowMeans(select(gene_com_meth, contains("_B"))), Li=rowMeans(select(gene_com_meth, contains("_L"))), Mu=rowMeans(select(gene_com_meth, contains("_M"))), Pl=rowMeans(select(gene_com_meth, contains("_P"))))

gene_com_exp_tis_mean  <- arrange(gene_com_exp_tis_mean,  rownames(gene_com_exp_tis_mean))

gene_com_meth_tis_mean <- arrange(gene_com_meth_tis_mean, rownames(gene_com_meth_tis_mean))

Br <- data.frame(Exp=gene_com_exp_tis_mean$Br, Meth=gene_com_meth_tis_mean$Br, Tissue = rep("Br", length(gene_com_exp_tis_mean$Br))) %>% set_rownames(rownames(gene_com_exp_tis_mean))

Li <- data.frame(Exp=gene_com_exp_tis_mean$Li, Meth=gene_com_meth_tis_mean$Li, Tissue = rep("Li", length(gene_com_exp_tis_mean$Li))) %>% set_rownames(rownames(gene_com_exp_tis_mean))

Mu <- data.frame(Exp=gene_com_exp_tis_mean$Mu, Meth=gene_com_meth_tis_mean$Mu, Tissue = rep("Mu", length(gene_com_exp_tis_mean$Mu))) %>% set_rownames(rownames(gene_com_exp_tis_mean))

Pl <- data.frame(Exp=gene_com_exp_tis_mean$Pl, Meth=gene_com_meth_tis_mean$Pl, Tissue = rep("Pl", length(gene_com_exp_tis_mean$Pl))) %>% set_rownames(rownames(gene_com_exp_tis_mean))

BrSpe <- filter(Tis_spe, Tissue=="Br")
LiSpe <- filter(Tis_spe, Tissue=="Li")
MuSpe <- filter(Tis_spe, Tissue=="Mu")
PlSpe <- filter(Tis_spe, Tissue=="Pl")

allTisSpe <- rbind(BrSpe, LiSpe, MuSpe, PlSpe)

###Brain
Specific <- c()
for(i in 1:nrow(Br)){
  if(rownames(Br)[i] %in% BrSpe$Gene){
    Specific[i] <- c("Spe")
  }else if(!(rownames(Br)[i] %in% BrSpe$Gene)&&((rownames(Br)[i] %in% allTisSpe$Gene))){
    Specific[i] <- c("Other_Spe")
  }else{
    Specific[i] <- c("NonSpe")
  }
}

Br$Specific <- Specific


##Liver
Specific <- c()
for(i in 1:nrow(Li)){
  if(rownames(Li)[i] %in% LiSpe$Gene){
    Specific[i] <- c("Spe")
  }else if(!(rownames(Li)[i] %in% LiSpe$Gene)&&((rownames(Li)[i] %in% allTisSpe$Gene))){
    Specific[i] <- c("Other_Spe")
  }else{
    Specific[i] <- c("NonSpe")
  }
}
Li$Specific <- Specific


##Muscle
Specific <- c()
for(i in 1:nrow(Mu)){
  if(rownames(Mu)[i] %in% MuSpe$Gene){
    Specific[i] <- c("Spe")
  }else if(!(rownames(Mu)[i] %in% MuSpe$Gene)&&((rownames(Mu)[i] %in% allTisSpe$Gene))){
    Specific[i] <- c("Other_Spe")
  }else{
    Specific[i] <- c("NonSpe")
  }
}
Mu$Specific <- Specific


##Placenta
Specific <- c()
for(i in 1:nrow(Pl)){
  if(rownames(Pl)[i] %in% PlSpe$Gene){
    Specific[i] <- c("Spe")
  }else if(!(rownames(Pl)[i] %in% PlSpe$Gene)&&((rownames(Pl)[i] %in% allTisSpe$Gene))){
    Specific[i] <- c("Other_Spe")
  }else{
    Specific[i] <- c("NonSpe")
  }
}
Pl$Specific <- Specific

allTis <- rbind(Br,Li,Mu,Pl)


#1st exon
data_text <- data.frame(label=c("rho=-0.546","rho=-0.471","rho=-0.524","rho=-0.467"),Tissue=c("Br","Li","Mu","Pl"), x=c(10,10,10,10),y=c(1.0,1.0,1.0,1.0))
#promoter
data_text <- data.frame(label=c("rho=-0.463","rho=-0.282","rho=-0.402","rho=-0.251"),Tissue=c("Br","Li","Mu","Pl"), x=c(10,10,10,10),y=c(1.0,1.0,1.0,1.0))


##plot
Br_plot <- ggplot(data = filter(Br, Specific=="NonSpe"),aes(log2(Exp), Meth)) + 
  geom_point(alpha = 0.6, shape = 21, color="gray") +
  labs(x = NULL, y = "Methylation Level", title = "Brain") + 
  scale_y_continuous(breaks = seq(0,1,length.out = 6)) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10)) +
  geom_smooth(data=Br, aes(x=log2(Exp), y=Meth), color = "blue") + 
  geom_point(data=filter(Br, Specific=="Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha=0.6, shape = 16, color = "#b30000") +
  geom_point(data=filter(Br, Specific=="Other_Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha = 0.6, shape = 18, color = "#fe9929") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=9.5, y=1, label = "rho = -0.546")
  

Li_plot <- ggplot(data = filter(Li, Specific=="NonSpe"),aes(log2(Exp), Meth)) + 
  geom_point(alpha = 0.6, shape = 21, color="gray") +
    labs(x = NULL, y = NULL, title = "Liver") + 
  scale_y_continuous(breaks = seq(0,1,length.out = 6)) +
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10, 15))+
  geom_smooth(data=Li, aes(x=log2(Exp), y=Meth), color = "blue") + 
  geom_point(data=filter(Li, Specific=="Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha=0.6, shape = 16, color = "#b30000") +
  geom_point(data=filter(Li, Specific=="Other_Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha = 0.6, shape = 18, color = "#fe9929") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))  +
  annotate("text",x=12.6, y=1, label = "rho = -0.471")


Mu_plot <- ggplot(data = filter(Mu, Specific=="NonSpe"),aes(log2(Exp), Meth)) + 
  geom_point(alpha = 0.6, shape = 21, color="gray") +
  labs(x = "Genes Expression log2(TPM)", y = "Methylation Level", title = "Muscle") + 
  scale_y_continuous(breaks = seq(0,1,length.out = 6)) +
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10))+
  geom_smooth(data=Mu, aes(x=log2(Exp), y=Meth), color = "blue") + 
  geom_point(data=filter(Mu, Specific=="Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha=0.6, shape = 16, color = "#b30000") +
  geom_point(data=filter(Mu, Specific=="Other_Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha = 0.6, shape = 18, color = "#fe9929") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))  +
  annotate("text",x=10, y=1, label = "rho = -0.524")


Pl_plot <- ggplot(data = filter(Pl, Specific=="NonSpe"),aes(log2(Exp), Meth)) + 
  geom_point(alpha = 0.6, shape = 21, color="gray") +
  labs(x = "Genes Expression log2(TPM)", y = NULL, title = "Placenta") + 
  scale_y_continuous(breaks = seq(0,1,length.out = 6)) + 
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10, 15))+
  geom_smooth(data=Pl, aes(x=log2(Exp), y=Meth), color = "blue")+
  geom_point(data=filter(Pl, Specific=="Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha=0.6, shape = 16, color = "#b30000") +
  geom_point(data=filter(Pl, Specific=="Other_Spe"), mapping=aes(x=log2(Exp), y=Meth), 
             alpha = 0.6, shape = 18, color = "#fe9929") + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text",x=11.2, y=1, label = "rho = -0.467")

ALL <- ggarrange(Br_plot,Li_plot,Mu_plot,Pl_plot, ncol = 2, nrow = 2)

pdf("ExpAndMethWithHlight_Promoter.pdf",width = 6, height = 6)
print(ALL)
dev.off()