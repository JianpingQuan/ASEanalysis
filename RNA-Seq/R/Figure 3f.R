library(dplyr)
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\02_POE基因位置分布")

input=c("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\02_POE基因位置分布\\data\\")

##pig chromosome length
chr_bed <- read.table(paste0(input,"pig_karyotype.txt"), header = 1)

##gene position information in reference
gene_list <- read.table(paste0(input,"RefGeneList.txt"), header = F, row.names = 1)

##AGE gene identified
AGE_gene <- read.csv(paste0(input,"AGE.list.csv"), header = T, row.names = 1)

AGE_gene_pos <- filter(gene_list, rownames(gene_list) %in% AGE_gene$all_gene_unique)
colnames(AGE_gene_pos) <- c("Chrom","Start", "End")

gene_chr <- unique(AGE_gene_pos$Chrom)
chr_bed_t <- filter(chr_bed, chr_bed$Chr %in% gene_chr) %>% arrange(desc(Chr)) %>% t()
nchr <- ncol(chr_bed_t)

colnames(chr_bed_t) <- chr_bed_t[1,]
chr_len <- chr_bed_t[3,]
chr_len_df <- as.data.frame(chr_len)

gene_position_changed <- c()
gene_chrom <- c()

AGE_gene_pos[,1] <- as.numeric(AGE_gene_pos[,1])
AGE_gene_pos_sort <- arrange(AGE_gene_pos,Chrom, Start)


for(i in 1:nrow(AGE_gene_pos_sort)){
  if(i==1){
    gene_chrom[i] = AGE_gene_pos_sort[i,1]
    gene_position_changed[i] = AGE_gene_pos_sort[i,2]
  } else if(i>1){
    
    if(AGE_gene_pos_sort[i,1] == AGE_gene_pos_sort[i-1,1]){
    
      if(AGE_gene_pos_sort[i,2] >= gene_position_changed[i-1]+1500000){
      
        gene_chrom[i] = AGE_gene_pos_sort[i,1]
        gene_position_changed[i] = AGE_gene_pos_sort[i,2]
   
      } else if(AGE_gene_pos_sort[i,2] < gene_position_changed[i-1]+1500000){
        gene_chrom[i] = AGE_gene_pos_sort[i,1]
        gene_position_changed[i] <- gene_position_changed[i-1]+1500000
      }
     
    } else if(AGE_gene_pos_sort[i,1] != AGE_gene_pos_sort[i-1,1]){
      gene_chrom[i] = AGE_gene_pos_sort[i,1]
      gene_position_changed[i] = AGE_gene_pos_sort[i,2]
    }
  }
}


AGE_gene_pos_sort_changed <- data.frame(Chrom=gene_chrom, Start= gene_position_changed)
rownames(AGE_gene_pos_sort_changed) <- rownames(AGE_gene_pos_sort)

AGE_gene_pos_sort_combine <- cbind(AGE_gene_pos_sort, AGE_gene_pos_sort_changed)
colnames(AGE_gene_pos_sort_combine) <- c("Chrom_ori","Start_ori", "End", "Chrom_alt", "Start_alt")

AGE_gene_pos_sort_final <- AGE_gene_pos_sort_combine

barplot(chr_len, horiz = TRUE, border = "white", col = "#FCF8CD", width = 2, space = rep(2.5,nchr),las=1, xlim = c(0,3.0e+08), xaxt="n")
axis(1,c(0.0e+00,5.0e+07,1.0e+08,1.5e+08,2.0e+08,2.5e+08,3.0e+08),labels=c("0","50Mb","100Mb","150Mb","200Mb","250Mb","300Mb"))

points(x=rep(1000000,nchr),y=seq(6,7*nchr,by=7),pch=20, col = "#7C1322", bg="#7C1322",cex = 3)

for(i in 1:nrow(AGE_gene_pos_sort_final)){
  ord <- which(as.numeric(rownames(chr_len_df))==AGE_gene_pos_sort_final$Chrom_ori[i])
  y_0 <- 5.2+7*(ord-1)
  y_1<- 6.8+7*(ord-1)
  x_0 <- AGE_gene_pos_sort_final$Start_ori[i]
  x_1 <- AGE_gene_pos_sort_final$Start_ori[i]
  segments(x0=x_0,y0=y_0,x1=x_1,y1=y_1, col = "#2c7bb6",lwd = 0.5)
  
  x_point <- AGE_gene_pos_sort_final$Start_alt[i]
  
  points(x=x_point,y=y_1+1,pch=20, col = "#2c7bb6", bg=NULL,cex = 1)
  
  segments(x0=x_0,y0=y_1+0.3,x1=x_point,y1=y_1+0.7, col = "#2c7bb6",lwd = 0.5)
}

