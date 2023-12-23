##Difference analysis of chromoatin accessibility 
setwd("F:\\Manuscript\\ATAC\\data\\Tissue_stage_specific")
library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)

dataCount <- read.table("F:\\Manuscript\\ATAC\\data\\AllPeakList.merged.TPM.txt", header = T, row.names=1)

datFilter <- select(dataCount, !contains(c("DL70_F_Li3","LD168_F_Br1","LD168_F_Li1","LD168_F_Mu1","LD168_F_Br2","LD168_F_Li2","LD168_F_Mu2","DL168_F_Br1","DL168_F_Li1","DL168_F_Mu1")))

ref_tissue <- c("Br","Li","Mu","Pl")
ref_stage <- c("40","70","115","168")


for(j in 1:4){
	dataFilter <- select(datFilter, contains(ref_stage[j]))
	
	if(j >=3){
		n = 3
	}else{
		n=4
	}
	for(i in 1:n){

	sex <- str_sub(colnames(dataFilter),-5,-5) %>% as.factor()
	tissue <- str_sub(colnames(dataFilter),-3,-2) %>% as.factor()
	tissue <- relevel(tissue, ref=ref_tissue[i])
	
	design <- model.matrix(~sex+tissue)
	
	y <- DGEList(counts=dataFilter,genes=rownames(dataFilter))

	#isexpr <- rowSums(cpm(y)>=3) >= 5
	#y <- y[isexpr, , keep.lib.sizes=FALSE]
	#y <- calcNormFactors(y,method = "TMM")
	
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

	##Compared with reference tissue, genes with logFC <=-4 in other 3 tissues were screened
	if(n==4){
		colnames(group_tissue) <- paste0("V",seq(1,7))
			a <- filter(group_tissue, V7<0.05) %>% filter(V1<=-1,V2<=-1,V3<=-1)
	}else{
		colnames(group_tissue) <- paste0("V",seq(1,6))
		a <- filter(group_tissue, V6<0.05) %>% filter(V1<=-1,V2<=-1)
	}
	
	write.csv(a, paste0("Tis_spe.peaks_basedOnRef",ref_tissue[i],"_Ref",ref_stage[j],".csv"), quote = F, row.names = T)
	}
}


##
setwd("F:\\Manuscript\\ATAC\\data\\Tissue_stage_specific")
library(dplyr)
path = c("F:\\Manuscript\\ATAC\\data\\Tissue_stage_specific")
tissue <- c("Br", "Li", "Mu", "Pl")

for(i in 1:4){
	pattern <- paste0("peaks_basedOnRef",tissue[i],"_")
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
write.csv(Tissue_specific_gene,"ATAC-Tissue_specific.peaks.csv", quote=F, row.names=F)

