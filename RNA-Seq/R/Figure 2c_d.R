
#-----------------------------Identification of tissue-specific and developmental stage-specific expressed genes---------------------------------------
#
#1. It was performed in separate periods to screen out tissue-specific expressed genes
setwd("F:\\ATAC\\processFiles\\RNA-tissue_specific0711")
library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)

#Reads count based gene expression matrix
dataExp <- read.csv("F:\\RNA\\DEG_Tissues_Stages\\GeneExpressionAllSamples.csv", header = T, row.names = 1) %>% select(., !contains(c("DL70_F_Li3","LD168_F_Br1","LD168_F_Li1","LD168_F_Mu1","LD168_F_Br2","LD168_F_Li2","LD168_F_Mu2","DL168_F_Br1","DL168_F_Li1","DL168_F_Mu1")))

ref_tissue <- c("Br","Li","Mu","Pl")
ref_stage <- c("40","70","115","168")

for(j in 1:4){
	dataFilter <- select(dataExp, contains(ref_stage[j]))
	
	if(j >=3){
		n = 3
	}else{
		n=4
	}
	for(i in 1:n){
	##set reference tissue
	sex <- str_sub(colnames(dataFilter),-5,-5) %>% as.factor()
	tissue <- str_sub(colnames(dataFilter),-3,-2) %>% as.factor()
	tissue <- relevel(tissue, ref=ref_tissue[i])
	print(levels(tissue))
	
	design <- model.matrix(~sex+tissue)
	
	y <- DGEList(counts=dataFilter,genes=rownames(dataFilter))
	isexpr <- rowSums(cpm(y)>=3) >= 5
	y <- y[isexpr, , keep.lib.sizes=FALSE]
	y <- calcNormFactors(y,method = "TMM")
	y <- estimateDisp(y, design, robust=TRUE)
	fit <- glmFit(y, design)
	if(n==4){
		lrt_tissue <- glmLRT(fit, coef=3:5)
	}else{
		lrt_tissue <- glmLRT(fit, coef=3:4)
	}
	
	group_tissue <- lrt_tissue$table

	FDR_tissue <- p.adjust(group_tissue$PValue, method="fdr")
	
	group_tissue <- cbind(group_tissue, FDR_tissue)

	##filter out the gene that was compared to reference tissue, the logFC <=-4 in other three kinds of tissues
	if(n==4){
		colnames(group_tissue) <- paste0("V",seq(1,7))
			a <- filter(group_tissue, V7<0.05) %>% filter(V1<=-4,V2<=-4,V3<=-4)
	}else{
		colnames(group_tissue) <- paste0("V",seq(1,6))
		a <- filter(group_tissue, V6<0.05) %>% filter(V1<=-4,V2<=-4)
	}
	
	write.csv(a, paste0("Tis_spe.genes_basedOn","Ref",ref_tissue[i],"_Ref",ref_stage[j],".csv"), quote = F, row.names = T)
	}
}

##Statistics of tissue-specific expressed genes
path = c("F:\\ATAC\\processFiles\\RNA-tissue_specific0711\\")
tissue <- c("Br", "Li", "Mu", "Pl")
for(i in 1:4){
	pattern <- paste0("genes_basedOnRef",tissue[i],"_")
	files <- list.files(path=path, pattern = pattern )
	data <- lapply(files, read.csv, header=T)
	if(i<=3){
		gene_data_frame <- as.data.frame(table(c(data[[1]]$X, data[[2]]$X, data[[3]]$X, data[[4]]$X)))
		gene_filter <- filter(gene_data_frame, Freq==4)
	}else{
		gene_data_frame <- as.data.frame(table(c(data[[1]]$X, data[[2]]$X)))
		gene_filter <- filter(gene_data_frame, Freq==2)
	}
	assign(paste0(tissue[i],"_specific"), gene_filter)
}

Br_specific$Tissue <- rep("Br", nrow(Br_specific))
Li_specific$Tissue <- rep("Li", nrow(Li_specific))
Mu_specific$Tissue <- rep("Mu", nrow(Mu_specific))
Pl_specific$Tissue <- rep("Pl", nrow(Pl_specific))

Tissue_specific_gene <- rbind(Br_specific, Li_specific, Mu_specific, Pl_specific)
write.csv(Tissue_specific_gene,"RNA-Tissue_specific.genes.csv", quote=F, row.names=F)


#2. Intersections of tissue-specific genes at different times (using TBtools to plot venn diagram)


#3. The developmental stage-specific expressed genes were screened in different tissues
setwd("F:\\ATAC\\processFiles\\RNA-tissue_specific0711")
library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)

dataExp <- read.csv("F:\\RNA\\DEG_Tissues_Stages\\GeneExpressionAllSamples.csv", header = T, row.names = 1)

ref_tissue <- c("Br","Li","Mu","Pl")
ref_stage <- c("40","70","115","168")

for(j in 1:4){
	dataFilter <- select(dataExp, contains(ref_tissue[j]))
	tmp <- unlist(str_split(colnames(dataFilter),"_"))
	
	if(j >3){
		n = 2
	}else{
		n=4
	}
	for(i in 1:n){
	##set the reference stage
	sex <- str_sub(colnames(dataFilter),-5,-5) %>% as.factor()
	stage <- str_sub(tmp[seq(1,length(tmp),3)],3) %>% as.factor()
	stage <- relevel(stage, ref=ref_stage[i])
	
	print(levels(stage))

	design <- model.matrix(~sex+stage)
	
	y <- DGEList(counts=dataFilter,genes=rownames(dataFilter))
	isexpr <- rowSums(cpm(y)>=3) >= 5
	y <- y[isexpr, , keep.lib.sizes=FALSE]
	y <- calcNormFactors(y,method = "TMM")
	y <- estimateDisp(y, design, robust=TRUE)
	fit <- glmFit(y, design)
	
	if(n==4){
		lrt_stage <- glmLRT(fit, coef=3:5)
	}else{
		lrt_stage <- glmLRT(fit, coef=3)
	}

	group_stage <- lrt_stage$table

	FDR_stage <- p.adjust(group_stage$PValue, method="fdr")
	
	group_stage <- cbind(group_stage, FDR_stage)

	#filter out the gene that was compared to reference stage, the logFC <=-4 in other three stage
	if(n==4){
		colnames(group_stage) <- paste0("V",seq(1,7))
		a <- filter(group_stage, V7<0.05) %>% filter(V1<=-4,V2<=-4,V3<=-4)
	}else{
		colnames(group_stage) <- paste0("V",seq(1,5))
		a <- filter(group_stage, V5<0.05) %>% filter(V1<=-4)
	}
	
	write.csv(a, paste0("Sta_spe.genes_basedOn","Ref",ref_stage[i],"_Ref",ref_tissue[j],".csv"), quote = F, row.names = T)
	}
}


#3. The tissue-specifci and developmental stage-specific expressed genes were visualized using IGV software
# The input file were bigwig files that were generated in pipeline.