setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\03_POE基因热图展示")
library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)
library(magrittr)
library(pheatmap)
### generate the bias score matrix grouping by stage
data_filter <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\A_GeneExpressionMatrix\\AutosomeGeneCount_filtered.csv",header = T, row.names = 1)

ASE <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\01_POE基因数量统计\\b_POE基因列表\\ASE_N-1_filter_final.PO2.csv", header = T, row.names = 1)

Gene_ref <- read.table("RefGeneList.txt", header = F, row.names = 1)

for(st in c("40","70","115","168")){

ase_stage <- select(ASE, contains(st))

## take the position of ASE and sort
gene_list <- c()
for(i in 1:ncol(ase_stage)){
  gene_list <- append(gene_list, ase_stage[,i])
}
gene_uniq <- sort(na.omit(unique(gene_list)))

## Extract the sample's gene expression data 
gene_uniq_expr <- data_filter %>% filter(rownames(data_filter) %in% gene_uniq) %>% select(contains(st))

##calculated the expression ratio mapping to paternal reference 
gene_uniq_expr_pat <- gene_uniq_expr %>% select(contains("Pat"))
gene_uniq_expr_mat <- gene_uniq_expr %>% select(contains("Mat"))
gene_uniq_expr_total <- gene_uniq_expr_pat + gene_uniq_expr_mat
gene_uniq_pat_bias_score <- gene_uniq_expr_pat/gene_uniq_expr_total

##change the column name 
colnames(gene_uniq_pat_bias_score) <- str_sub(colnames(gene_uniq_expr_pat), 1, -5)

##select the stage
stage <- select(gene_uniq_pat_bias_score, contains(st))

##remove the gene not expressed in all samples of this stage
stage_DL <- stage %>% select(contains("DL")) %>% filter(!if_all(.fns = is.na))
stage_LD <- stage %>% select(contains("LD")) %>% filter(!if_all(.fns = is.na))

gene_final <- intersect(rownames(stage_DL),rownames(stage_LD))
stage_final <- filter(stage, rownames(stage) %in% gene_final)

###
gene_uniq_pos <- filter(Gene_ref, rownames(Gene_ref) %in% gene_final) %>% arrange(rownames(.))
colnames(gene_uniq_pos) <- c("Chrom", "Start", "End")
stage_final_pos <- cbind(stage_final, gene_uniq_pos)
gene_uniq_pos_sort <- arrange(stage_final_pos, Chrom, Start)
gene_uniq_pos_sort_final <- gene_uniq_pos_sort[,1:(ncol(gene_uniq_pos_sort)-3)]

if(st=="40"){
	ord <- c(1,2,3,10,11,12,21,22,23,33,34,35,4,5,6,13,14,15,24,25,26,36,37,38,7,8,9,16,17,18,27,28,29,39,40,41,19,20,30,31,32,42,43,44)
}else if(st=="70"){
	ord <- c(1,2,3,12,13,14,23,24,25,34,35,36,4,5,6,15,16,26,27,28,37,38,39,7,8,9,17,18,19,29,30,31,40,41,42,10,11,20,21,22,32,33,43,44,45)
}else if(st=="115"){
	ord <- c(1,2,3,10,11,12,19,20,21,28,29,30,4,5,6,13,14,15,22,23,24,31,32,33,7,8,9,16,17,18,25,26,27,34,35,36)
}else if(st=="168"){
	ord <- c(1,4,7,10,13,16,19,22,25,2,5,8,11,14,17,20,23,26,3,6,9,12,15,18,21,24,27)
}

gene_uniq_pos_sort_final2 <- gene_uniq_pos_sort_final[,ord]

write.csv(gene_uniq_pos_sort_final2,paste0("Stage",st,"_heatmap.csv"), quote = F, row.names = T)
}


###Plot the heatmap 
output=c("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\03_POE基因热图展示\\figure")

dat <- read.csv("Stage70_heatmap.csv", header=T, row.names=1)

cc1 = grDevices::colorRampPalette(c("#7C1322","#F2F2F2","#143C7B") %>% ggplot2::alpha(0.8))
####
p1 <- pheatmap::pheatmap(dat,
                         cluster_rows = F,
                         cluster_cols=F,
                         #gaps_col = c(12,24,36),
						 #gaps_col = c(12,23,35),
						 gaps_col = c(12,24),
						 #gaps_col = c(9,18),
                         cellwidth=10,
                         cellheight=10,
                         color = cc1(30),
                         na_col="gray",
                         show_rownames = T,
                         show_colnames = T)
