##POE
stage <- c("40","70","115","168")
tissue <- c("Br", "Li", "Mu", "Pl")

library(customLayout)
library(dplyr)


lay <- lay_new(matrix(c(1,2,3,4,5,6,7,8,9,10,11,0,12,13,14,0), ncol=4, byrow = T), widths = c(1,1,1,1), heights = c(1,1,1,1))
lay_show(lay)

a <- c()
for(i in 1:4){

  if(i == 3||i==4){
    n=3
  }else{
    n=4
  }
  
  for(j in 1:n){
    RNA <- read.table(paste0("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result\\",tissue[j],"_",stage[i],"_bias_final.txt"),header = F)
    RNA_POE <- RNA$V1
    RNA_data <- read.csv(paste0("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result\\",tissue[j],"_",stage[i],"_AllGenesPO2.csv"),header = T, row.names=1)
	
    peak_gene <- read.csv(paste0("F:\\ATAC\\processFiles\\POE05112023\\ATAC_POE\\method2\\FilterDistal\\",stage[i],tissue[j],"DiffPeaks_peakGene_POE.method2_nodistalAll.csv"), header = T, row.names = 1)
    com <- intersect(RNA_POE, rownames(peak_gene))
    
    RNA_POE_com <- filter(RNA_data, rownames(RNA_data) %in% com)  %>% arrange(.,rownames(.))
    Peak_POE_com <- filter(peak_gene, rownames(peak_gene) %in% com) %>% arrange(.,rownames(.))
    
    ATAC_fc=Peak_POE_com$logFC
    RNA_fc=RNA_POE_com$logFC
	
	RNA_ATAC_combine <- data.frame(ATAC=2^ATAC_fc/(1+2^ATAC_fc), RNA=2^RNA_fc/(1+2^RNA_fc))
	rownames(RNA_ATAC_combine) <- rownames(RNA_POE_com)
	#write.csv(RNA_ATAC_combine, paste0(stage[i],tissue[j],"_RNA_ATAC_biasScorePOE.csv"), quote = F, row.names = TRUE)
	
	#Filter the number of points in the first and third quadrants
	RNA_ATAC_1_3 <- filter(RNA_ATAC_combine, (RNA<0.5&ATAC<0.5)|(RNA>0.5&ATAC>0.5))
	#write.csv(RNA_ATAC_1_3, paste0(stage[i],tissue[j],"_RNA_ATAC_biasScore_contPOE.csv"), quote = F, row.names = TRUE)
	print(paste0(stage[i],tissue[j]))
	print(nrow(RNA_ATAC_1_3)/nrow(RNA_ATAC_combine))
	a <- c(a,nrow(RNA_ATAC_1_3)/nrow(RNA_ATAC_combine))
	
	
	x <- 2^ATAC_fc/(1+2^ATAC_fc)
	y <- 2^RNA_fc/(1+2^RNA_fc)
	
	#Spearman correlation coefficient was calculated
	cor_res <- cor.test(x, y, method = "spearman")
	r_spearman <- round(cor_res$estimate, 2)
	colors <- ifelse((x >= 0.5 & y >= 0.5) | (x < 0.5 & y < 0.5), "#ED0000FF", "gray")
	
    plot(x, y, pch = 16, cex= 1, col = colors, xlim=c(0,1), xlab=paste0("ATAC ",tissue[j],stage[i]), ylab = paste0("RNA ",tissue[j],stage[i]))
    segments(x0=0.5,x1=0.5,y0=0,y1=1, col = "blue",lty = "dashed")
    segments(y0=0.5,y1=0.5,x0=0,x1=1, col = "blue",lty = "dashed")
	text(x = 0.2, y = 0.9, labels = paste0("r = ", r_spearman))
  }
}



##AGE
stage <- c("40","70","115","168")
tissue <- c("Br", "Li", "Mu", "Pl")

library(customLayout)
library(dplyr)


lay <- lay_new(matrix(c(1,2,3,4,5,6,7,8,9,10,11,0,12,13,14,0), ncol=4, byrow = T), widths = c(1,1,1,1), heights = c(1,1,1,1))
lay_show(lay)

b <- c()
for(i in 1:4){

  if(i == 3||i==4){
    n=3
  }else{
    n=4
  }
  
  for(j in 1:n){
    RNA <- read.table(paste0("F:\\RNA_20221101\\02_ModifiedGenome\\D_AGEgeneSelection\\result\\",tissue[j],"_",stage[i],".txt"),header = F)
    RNA_AGE <- RNA$V1
	
    RNA_data <- read.csv(paste0("F:\\RNA_20221101\\02_ModifiedGenome\\D_AGEgeneSelection\\result\\",tissue[j],"_",stage[i],"_AllGenesAG2.csv"),header = T, row.names=1)
	
    peak_gene <- read.csv(paste0("F:\\ATAC\\processFiles\\POE05112023\\ATAC_AGE\\method2\\",stage[i],tissue[j],"DiffPeaks_peakGene_AGE.method2_nodistalAll.csv"), header = T, row.names = 1)
    com <- intersect(RNA_AGE, rownames(peak_gene))
    
    RNA_AGE_com <- filter(RNA_data, rownames(RNA_data) %in% com)  %>% arrange(.,rownames(.))
    Peak_AGE_com <- filter(peak_gene, rownames(peak_gene) %in% com) %>% arrange(.,rownames(.))
    
    ATAC_fc=Peak_AGE_com$logFC
    RNA_fc=RNA_AGE_com$logFC
	
	RNA_ATAC_combine <- data.frame(ATAC=2^ATAC_fc/(1+2^ATAC_fc), RNA=2^RNA_fc/(1+2^RNA_fc))
	rownames(RNA_ATAC_combine) <- rownames(RNA_AGE_com)
	#write.csv(RNA_ATAC_combine, paste0(stage[i],tissue[j],"_RNA_ATAC_biasScoreAGE.csv"), quote = F, row.names = TRUE)
	
	RNA_ATAC_1_3 <- filter(RNA_ATAC_combine, (RNA<0.5&ATAC<0.5)|(RNA>0.5&ATAC>0.5))
	#write.csv(RNA_ATAC_1_3, paste0(stage[i],tissue[j],"_RNA_ATAC_biasScore_contAGE.csv"), quote = F, row.names = TRUE)
	print(paste0(stage[i],tissue[j]))
	print(nrow(RNA_ATAC_1_3)/nrow(RNA_ATAC_combine))
	b <- c(b,nrow(RNA_ATAC_1_3)/nrow(RNA_ATAC_combine))
	
	
    x <- 2^ATAC_fc/(1+2^ATAC_fc)
	y <- 2^RNA_fc/(1+2^RNA_fc)
	
	cor_res <- cor.test(x, y, method = "spearman")
	r_spearman <- round(cor_res$estimate, 2)
	colors <- ifelse((x >= 0.5 & y >= 0.5) | (x < 0.5 & y < 0.5), "#ED0000FF", "gray")
	
    plot(x, y, pch = 16, cex= 1, col = colors, xlim=c(0,1), xlab=paste0("ATAC ",tissue[j],stage[i]), ylab = paste0("RNA ",tissue[j],stage[i]))
    segments(x0=0.5,x1=0.5,y0=0,y1=1, col = "blue",lty = "dashed")
    segments(y0=0.5,y1=0.5,x0=0,x1=1, col = "blue",lty = "dashed")
	text(x = 0.2, y = 0.9, labels = paste0("r = ", r_spearman))
  }
}
