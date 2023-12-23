
#Smooth line plot to show the relationship between methylation and expression of global gene
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(magrittr)

dataExp <- read.csv("F:\\RNA\\RawTables\\OriginalRef\\GeneCountMatrix.tpmUnit.csv", header=T, row.names=1) %>% select(., contains("70"))

feature <- c("geneBody", "gene1stExon", "geneIntron", "genePromoter","geneUTR5","geneUTR3")
i= 2
path <- paste0("C:\\Users\\qjp\\Desktop\\genomeFeature20220823\\data\\AllSample.", feature[i],".methyLevel.txt")
datMeth <- read.table(path, header = F,row.names = 1)
sampleName <- read.table("C:\\Users\\qjp\\Desktop\\genomeFeature20220823\\data\\sampleName.txt", header = F)
colnames(datMeth) <- sampleName$V1

## get the commone gene between exp and meth
geneCom <- intersect(rownames(dataExp), rownames(datMeth))

gene_exp_mean <- data.frame(Br=rowMeans(select(dataExp, contains("Br"))), Li=rowMeans(select(dataExp, contains("Li"))), Mu=rowMeans(select(dataExp, contains("Mu"))), Pl=rowMeans(select(dataExp, contains("Pl")))) %>% set_rownames(rownames(dataExp)) %>% filter(.,rownames(.) %in% geneCom) %>% arrange(., rownames(.))

gene_meth_mean <- data.frame(Br=rowMeans(select(datMeth, contains("_B"))), Li=rowMeans(select(datMeth, contains("_L"))), Mu=rowMeans(select(datMeth, contains("_M"))), Pl=rowMeans(select(datMeth, contains("_P")))) %>% set_rownames(rownames(datMeth)) %>% filter(.,rownames(.) %in% geneCom) %>% arrange(., rownames(.))

## rank gene according to gene expressione value
gene_exp_mean_rank <- data.frame(BrO=rank(gene_exp_mean$Br), LiO=rank(gene_exp_mean$Li), MuO = rank(gene_exp_mean$Mu), PlO=rank(gene_exp_mean$Pl))

## get each tissue exp and meth arrange by exp rank
Br_Exp_Meth <- data.frame(Exp=gene_exp_mean$Br, Rank=gene_exp_mean_rank$BrO, Meth=gene_meth_mean$Br) %>% set_rownames(rownames(gene_exp_mean)) %>% arrange(.,Rank)

Li_Exp_Meth <- data.frame(Exp=gene_exp_mean$Li, Rank=gene_exp_mean_rank$LiO, Meth=gene_meth_mean$Li) %>% set_rownames(rownames(gene_exp_mean)) %>% arrange(.,Rank)

Mu_Exp_Meth <- data.frame(Exp=gene_exp_mean$Mu, Rank=gene_exp_mean_rank$MuO, Meth=gene_meth_mean$Mu) %>% set_rownames(rownames(gene_exp_mean)) %>% arrange(.,Rank)

Pl_Exp_Meth <- data.frame(Exp=gene_exp_mean$Pl, Rank=gene_exp_mean_rank$PlO, Meth=gene_meth_mean$Pl) %>% set_rownames(rownames(gene_exp_mean)) %>% arrange(.,Rank)

##correlation analysis
cor.test(Br_Exp_Meth$Exp,Br_Exp_Meth$Meth, method = "spearman")
cor.test(Li_Exp_Meth$Exp,Li_Exp_Meth$Meth, method = "spearman")
cor.test(Mu_Exp_Meth$Exp,Mu_Exp_Meth$Meth, method = "spearman")
cor.test(Pl_Exp_Meth$Exp,Pl_Exp_Meth$Meth, method = "spearman")

Br_Exp_Meth$Rank <- seq(1,nrow(Br_Exp_Meth))
Li_Exp_Meth$Rank <- seq(1,nrow(Li_Exp_Meth))
Mu_Exp_Meth$Rank <- seq(1,nrow(Mu_Exp_Meth))
Pl_Exp_Meth$Rank <- seq(1,nrow(Pl_Exp_Meth))

exp <- c()
for(i in 1:nrow(Br_Exp_Meth)){
  if(Br_Exp_Meth$Exp[i]==0){
    exp[i] <- 0.0001
  }else{
    exp[i] <- Br_Exp_Meth$Exp[i]
  }
}
Br_Exp_Meth$Exp <- exp

exp <- c()
for(i in 1:nrow(Li_Exp_Meth)){
  if(Li_Exp_Meth$Exp[i]==0){
    exp[i] <- 0.0001
  }else{
    exp[i] <- Li_Exp_Meth$Exp[i]
  }
}
Li_Exp_Meth$Exp <- exp

exp <- c()
for(i in 1:nrow(Mu_Exp_Meth)){
  if(Mu_Exp_Meth$Exp[i]==0){
    exp[i] <- 0.0001
  }else{
    exp[i] <- Mu_Exp_Meth$Exp[i]
  }
}
Mu_Exp_Meth$Exp <- exp


exp <- c()
for(i in 1:nrow(Pl_Exp_Meth)){
  if(Pl_Exp_Meth$Exp[i]==0){
    exp[i] <- 0.0001
  }else{
    exp[i] <- Pl_Exp_Meth$Exp[i]
  }
}
Pl_Exp_Meth$Exp <- exp


###
ylim.prim <- c(0, 100)   # in this example, precipitation
ylim.sec <- c(-10, 6)

Br <- ggplot(data = Br_Exp_Meth,aes(x=Rank))+ 
  geom_smooth(aes(y=Meth*100), color = "#3B4992")+
  geom_line(aes(y=log10(Exp)*10+40),color="#EE2200", size = 1)+ 
  theme_bw()+ xlab(NULL) + 
  scale_y_continuous(name="Gene Promoter Methylation level",limits=c(0,100), sec.axis = sec_axis(~(.*0.1)-4, name = NULL))+
  theme(axis.title.y = element_text(color ="#3B4992", size=13),axis.title.y.right = element_text(color = '#EE2200', size=13)) + 
  geom_vline(aes(xintercept=nrow(Br_Exp_Meth)/2),color = "#631779",linetype="dashed",size = 0.5) +
  annotate("text", x = 16000, y = 95, label = "r1 = -0.546",color="black",size = 4, fontface="bold" ) + labs(title = "Brain") + theme(plot.title = element_text(hjust = 0.5))
  
Li <- ggplot(data = Li_Exp_Meth,aes(x=Rank))+ 
  geom_smooth(aes(y=Meth*100), color = "#3B4992")+
  geom_line(aes(y=log10(Exp)*10+40),color="#EE2200", size = 1)+ 
  theme_bw()+ xlab(NULL) + 
  scale_y_continuous(name=NULL,limits=c(0,100), sec.axis = sec_axis(~(.*0.1)-4, name = "Gene Expression log10(TPM)"))+
  theme(axis.title.y = element_text(color ="#3B4992", size=13),axis.title.y.right = element_text(color = '#EE2200', size=13)) + 
  geom_vline(aes(xintercept=nrow(Li_Exp_Meth)/2),color = "#631779",linetype="dashed",size = 0.5) +
  annotate("text", x = 16000, y = 95, label = "r1 = -0.471",color="black",size = 4, fontface="bold" ) + labs(title = "Liver") + theme(plot.title = element_text(hjust = 0.5))

Mu <- ggplot(data = Mu_Exp_Meth,aes(x=Rank))+ 
  geom_smooth(aes(y=Meth*100), color = "#3B4992")+
  geom_line(aes(y=log10(Exp)*10+40),color="#EE2200", size = 1)+ 
  theme_bw()+ xlab("Rank Gene by Gene Expression Level") + 
  scale_y_continuous(name="Gene Promoter Methylation level",limits=c(0,100), sec.axis = sec_axis(~(.*0.1)-4, name = NULL))+
  theme(axis.title.y = element_text(color ="#3B4992", size=13),axis.title.y.right = element_text(color = '#EE2200', size=13)) + 
  geom_vline(aes(xintercept=nrow(Mu_Exp_Meth)/2),color = "#631779",linetype="dashed",size = 0.5) +
  annotate("text", x = 16000, y = 95, label = "r1 = -0.524",color="black",size = 4, fontface="bold" ) + labs(title = "Muscle") + theme(plot.title = element_text(hjust = 0.5))

Pl <- ggplot(data = Pl_Exp_Meth,aes(x=Rank))+ 
  geom_smooth(aes(y=Meth*100), color = "#3B4992")+
  geom_line(aes(y=log10(Exp)*10+40),color="#EE2200", size = 1)+ 
  theme_bw()+ xlab("Rank Gene by Gene Expression Level") + 
  scale_y_continuous(name=NULL,limits=c(0,100), sec.axis = sec_axis(~(.*0.1)-4, name = "Gene Expression log10(TPM)"))+
  theme(axis.title.y = element_text(color ="#3B4992", size=13),axis.title.y.right = element_text(color = '#EE2200', size=13)) + 
  geom_vline(aes(xintercept=nrow(Pl_Exp_Meth)/2),color = "#631779",linetype="dashed",size = 0.5) +
  annotate("text", x = 16000, y = 95, label = "r1 = -0.467",color="black",size = 4, fontface="bold" ) + labs(title = "Placenta") + theme(plot.title = element_text(hjust = 0.5))

ggarrange(Br,Li,Mu,Pl,ncol = 2, nrow = 2)