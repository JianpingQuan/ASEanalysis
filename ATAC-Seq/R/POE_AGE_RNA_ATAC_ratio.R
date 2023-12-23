
stage <- c("40","70","115","168")
tissue <- c("Br", "Li", "Mu", "Pl")
POE <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result\\POE.list.csv", header = T, row.names = 1)

library(customLayout)
library(dplyr)

##Set canvas layout
lay <- lay_new(matrix(c(1,2,3,4,5,6,7,8,9,10,11,0,12,13,14,0), ncol=4, byrow = T), widths = c(1,1,1,1), heights = c(1,1,1,1))
lay_show(lay)

for(i in 1:4){
  ##Determine whether the satges are D1 and D168
  if(i == 3||i==4){
    n=3
  }else{
    n=4
  }
  
  for(j in 1:n){
    RNA <- read.csv(paste0("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result\\",tissue[j],"_",stage[i],"_bias.csv"),header = T, row.names=1)
    RNA_POE <- rownames(RNA)
    
    peak_gene <- read.csv(paste0("C:\\Users\\qjp\\Desktop\\POE_AGE\\",stage[i],tissue[j],"DiffPeaks_peakGene_NoDistal_POE1101.csv"), header = T, row.names=1)
    
    #peak_gene_filter <- filter(peak_gene, PValue<=0.05)
    
    com <- intersect(RNA_POE, rownames(peak_gene))
    
    RNA_POE_com <- filter(RNA, rownames(RNA) %in% com) %>% arrange(., rownames(.))
    
    Peak_POE_com <- filter(peak_gene, rownames(peak_gene) %in% com) %>% arrange(., rownames(.))
    
    ATAC=Peak_POE_com$logFC
    RNA=RNA_POE_com$logFC
    
	df <- data.frame(rna=RNA,atac=ATAC)
	
	df_same <- filter(df, (rna<0&atac<0)|(rna>0&atac>0))
	
	df <- data.frame(rna=RNA_POE_com$logFC, atac=Peak_POE_com$logFC)
	
	df_same <- filter(df, (rna<0&atac<0)|(rna>0&atac>0))
	print(paste0(tissue[j],"_",stage[i]))
	print(nrow(df))
	print(nrow(df_same))
	
    plot(2^ATAC/(1+2^ATAC), 2^RNA/(1+2^RNA), xlim=c(0,1), xlab=paste0("ATAC ",tissue[j],stage[i]), ylab = paste0("RNA ",tissue[j],stage[i]))    
    segments(x0=0.5,x1=0.5,y0=0,y1=1, col = "blue")
    segments(y0=0.5,y1=0.5,x0=0,x1=1, col = "blue")
  }
}


###AGE

AGE <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\D_AGEgeneSelection\\result\\AGE.list.csv", header = T, row.names = 1)

library(customLayout)
library(dplyr)

lay <- lay_new(matrix(c(1,2,3,4,5,6,7,8,9,10,11,0,12,13,14,0), ncol=4, byrow = T), widths = c(1,1,1,1), heights = c(1,1,1,1))
lay_show(lay)

for(i in 1:4){

  if(i == 3||i==4){
    n=3
  }else{
    n=4
  }
  
  for(j in 1:n){
    RNA <- read.csv(paste0("F:\\RNA_20221101\\02_ModifiedGenome\\D_AGEgeneSelection\\result\\",tissue[j],"_",stage[i],"_bias.csv"),header = T, row.names=1)
    RNA_AGE <- rownames(RNA)
    
    peak_gene <- read.csv(paste0("C:\\Users\\qjp\\Desktop\\POE_AGE\\",stage[i],tissue[j],"DiffPeaks_peakGene_Nofilter_AG1101.csv"), header = T, row.names=1)
    
    #peak_gene_filter <- filter(peak_gene, PValue<=0.05)
    
    com <- intersect(RNA_AGE, rownames(peak_gene))
    
    RNA_AGE_com <- filter(RNA, rownames(RNA) %in% com) %>% arrange(., rownames(.))
    
    Peak_AGE_com <- filter(peak_gene, rownames(peak_gene) %in% com) %>% arrange(., rownames(.))
    
    ATAC=Peak_AGE_com$logFC
    RNA=RNA_AGE_com$logFC
    
	df <- data.frame(rna=RNA,atac=ATAC)
	
	df_same <- filter(df, (rna<0&atac<0)|(rna>0&atac>0))
	
	df <- data.frame(rna=RNA_AGE_com$logFC, atac=Peak_AGE_com$logFC)
	
	df_same <- filter(df, (rna<0&atac<0)|(rna>0&atac>0))
	print(paste0(tissue[j],"_",stage[i]))
	print(nrow(df))
	print(nrow(df_same))
	
    plot(2^ATAC/(1+2^ATAC), 2^RNA/(1+2^RNA), xlim=c(0,1), xlab=paste0("ATAC ",tissue[j],stage[i]), ylab = paste0("RNA ",tissue[j],stage[i]))    
    segments(x0=0.5,x1=0.5,y0=0,y1=1, col = "blue")
    segments(y0=0.5,y1=0.5,x0=0,x1=1, col = "blue")
  }
}