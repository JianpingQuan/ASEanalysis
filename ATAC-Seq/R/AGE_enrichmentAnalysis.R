##AGE gene enrichement analysis
dat <- read.csv("F:\\ATAC\\processFiles\\POE05112023\\MG_AllGenomePeakReadNumber.NoDistal.GroupByGene.csv", header = T, row.names = 1)
geneList <- rownames(dat)

AGE <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\D_AGEgeneSelection\\result\\AGE.list.csv", header = T, row.names = 1)
AGEList <- AGE$all_gene_unique

#cat *All.csv|tr ',' '\t'|awk '$6<=0.05'|cut -f1|sort|uniq > allAGE_peakGene.nodistal.txt

PeakGene <- read.table("F:\\ATAC\\processFiles\\POE05112023\\allAGE_peakGene.nodistal.txt", header = F)
PeakGeneList <- PeakGene$V1


ag <- c()

for(i in 1:10000){
  tmp <- sample(geneList, size = 394, replace = FALSE)
  overlap <- intersect(PeakGeneList, tmp)
  ag <- c(ag, length(overlap))
}


##overlap between AGE and PeakGene #51
length(intersect(AGEList, PeakGeneList))


#enrichement score #1.345378
(score_ag = 42/mean(ag))

#pvalue #0.02889711
(pvalue_ag = (sum(ag>=42)+1)/(10000+1))

age <- data.frame(Value=ag, Type=rep("AGE", length(ag)))


poe_age <- rbind(poe,age)

library(ggplot2)
library(ggsci)
library(ggbeeswarm)

p <- ggplot(poe_age, aes(x=Type, y=Value))+geom_boxplot(aes(fill = Type), alpha=0.8)+scale_fill_lancet()+theme_bw()+ guides(color = "none")+ylab("Overlap gene number") + xlab("Type") + labs(title="overlap number between randomly gene set of\nexpression and POE/AGE chromatin accessibility gene")+ theme(axis.title = element_text(size = 16)) + theme(axis.text = element_text(size = 14)) + geom_hline(yintercept = 42, linetype = "dashed", color = "blue") + geom_hline(yintercept = 25, linetype = "dashed", color = "red")+theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))


pdf("Enrichment_Analysis_RNA_ATAC.pdf", height=4, width=5)
p
dev.off()