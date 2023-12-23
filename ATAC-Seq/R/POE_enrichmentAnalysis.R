##POE gene enrichement analysis
dat <- read.csv("F:\\ATAC\\processFiles\\POE05112023\\MG_AllGenomePeakReadNumber.NoDistal.GroupByGene.csv", header = T, row.names = 1)
geneList <- rownames(dat)

POE <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result\\POE.list2.csv", header = T, row.names = 1)
POEList <- POE$all_gene_unique

PeakGene <- read.table("F:\\ATAC\\processFiles\\POE05112023\\allPOE_peakGene.nodistal.txt", header = F)
PeakGeneList <- PeakGene$V1


po <- c()

for(i in 1:10000){
  tmp <- sample(geneList, size = 179, replace = FALSE)
  overlap <- intersect(PeakGeneList, tmp)
  po <- c(po, length(overlap))
}

##overlap between POElist and Peakgenelist #25
length(intersect(POEList, PeakGeneList))

#enrichement score #2.100611
(score_po = 25/mean(po))

#pvalue #0.00039996
(pvalue_po = (sum(po>=25)+1)/(10000+1))

poe <- data.frame(Value=po, Type=rep("POE", length(po)))
