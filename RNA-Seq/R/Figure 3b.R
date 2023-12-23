library(dplyr)
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\02_POE基因位置分布")

input=c("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\02_POE基因位置分布\\data\\")

##pig chromosome length
chr_bed <- read.table(paste0(input,"pig_karyotype.txt"), header = 1)

##gene position information in reference
gene_list <- read.table(paste0(input,"RefGeneList.txt"), header = T, row.names = 1)

##POE gene identified
stage_gene <- read.csv(paste0(input,"AllSampleGroupByTissue_heatmap.csv"), header = T, row.names = 1)

##imprinting gene database
im_db <- read.csv(paste0(input,"Im.csv"),header = T, row.names = 1)

stage_gene_pos <- filter(gene_list, rownames(gene_list) %in% rownames(stage_gene))
colnames(stage_gene_pos) <- c("Chrom","Start", "End")

im_db_pos <- filter(gene_list, rownames(gene_list) %in% im_db$Ens_id)
colnames(im_db_pos) <- c("Chrom","Start", "End")

gene_chr <- unique(stage_gene_pos$Chrom)
chr_bed_t <- filter(chr_bed, chr_bed$Chr %in% gene_chr) %>% arrange(desc(Chr)) %>% t()

nchr <- ncol(chr_bed_t)
im_db_filter <- filter(im_db_pos, Chrom %in% gene_chr)

colnames(chr_bed_t) <- chr_bed_t[1,]
chr_len <- chr_bed_t[3,]
chr_len_df <- as.data.frame(chr_len)

gene_position_changed <- c()
gene_chrom <- c()

stage_gene_pos[,1] <- as.numeric(stage_gene_pos[,1])
stage_gene_pos_sort <- arrange(stage_gene_pos,Chrom, Start)

for(i in 1:nrow(stage_gene_pos_sort)){
  if(i==1){
    gene_chrom[i] = stage_gene_pos_sort[i,1]
    gene_position_changed[i] = stage_gene_pos_sort[i,2]
  } else if(i>1){
    ##Determine if the chromosomes are the same
    if(stage_gene_pos_sort[i,1] == stage_gene_pos_sort[i-1,1]){
      ##In the case of the same chromosome, the position of the gene before and after is determined to be less than 1.5M
      if(stage_gene_pos_sort[i,2] >= gene_position_changed[i-1]+1500000){
        ##The gene distance before and after the same chromosome is greater than 1.5M
        gene_chrom[i] = stage_gene_pos_sort[i,1]
        gene_position_changed[i] = stage_gene_pos_sort[i,2]
        ##The gene distance before and after the same chromosome is smaller than 1.5M
      } else if(stage_gene_pos_sort[i,2] < gene_position_changed[i-1]+1500000){
        gene_chrom[i] = stage_gene_pos_sort[i,1]
        gene_position_changed[i] <- gene_position_changed[i-1]+1500000
      }
      ##Retention of gene location information of different chromosomes
    } else if(stage_gene_pos_sort[i,1] != stage_gene_pos_sort[i-1,1]){
      gene_chrom[i] = stage_gene_pos_sort[i,1]
      gene_position_changed[i] = stage_gene_pos_sort[i,2]
    }
  }
}

stage_gene_pos_sort_changed <- data.frame(Chrom=gene_chrom, Start= gene_position_changed)
rownames(stage_gene_pos_sort_changed) <- rownames(stage_gene_pos_sort)

stage_gene_pos_sort_combine <- cbind(stage_gene_pos_sort, stage_gene_pos_sort_changed)
colnames(stage_gene_pos_sort_combine) <- c("Chrom_ori","Start_ori", "End", "Chrom_alt", "Start_alt")


#########################Modify the imprinted gene coordinates in the imprinted gene database#####################

im_db_filter[,1] <- as.numeric(im_db_filter[,1])
im_db_filter_sort <- arrange(im_db_filter, Chrom, Start)
im_position_changed <- c()
im_chrom <- c()

for(j in 1:nrow(im_db_filter_sort)){
  if(j == 1){
    im_chrom[j] = im_db_filter_sort[j,1]
    im_position_changed[j] = im_db_filter_sort[j,2]
  } else if(j > 1){
    ##Determine if the chromosomes are the same
    if(im_db_filter_sort[j,1] == im_db_filter_sort[j-1,1]){
      ##In the case of the same chromosome, the position of the gene before and after is determined to be less than 1.5M
      if(im_db_filter_sort[j,2] >= im_position_changed[j-1]+1500000){
        ##The gene distance before and after the same chromosome is greater than 1.5M
        im_chrom[j] = im_db_filter_sort[j,1]
        im_position_changed[j] = im_db_filter_sort[j,2]
        ##The gene distance before and after the same chromosome is smaller than 1.5M
      } else if(im_db_filter_sort[j,2] < im_position_changed[j-1]+1500000){
        im_chrom[j] = im_db_filter_sort[j,1]
        im_position_changed[j] = im_position_changed[j-1]+1500000
      }
      ##Retention of gene location information of different chromosomes
    } else if(im_db_filter_sort[j,1] != im_db_filter_sort[j-1,1]){
      im_chrom[j] = im_db_filter_sort[j,1]
      im_position_changed[j] = im_db_filter_sort[j,2]
    }
  }
}

im_db_filter_sort_changed <- data.frame(Chrom=im_chrom, Start= im_position_changed)
rownames(im_db_filter_sort) <- rownames(im_db_filter_sort)

im_db_filter_sort_combine <- cbind(im_db_filter_sort, im_db_filter_sort_changed)
colnames(im_db_filter_sort_combine) <- c("Chrom_ori","Start_ori", "End", "Chrom_alt", "Start_alt")


#####The overlapping genes of POE gene and imprinted gene pool were screened
overlap <- intersect(rownames(stage_gene_pos_sort_combine),rownames(im_db_filter_sort_combine))
overlap_pos <- filter(stage_gene_pos_sort_combine,rownames(stage_gene_pos_sort_combine) %in% overlap)

##The overlapping genes in POE gene and imprinted gene bank were removed respectively
stage_gene_pos_sort_final <- filter(stage_gene_pos_sort_combine, !(rownames(stage_gene_pos_sort_combine) %in% overlap))
im_db_filter_sort_final <- filter(im_db_filter_sort_combine, !(rownames(im_db_filter_sort_combine) %in% overlap))


barplot(chr_len, horiz = TRUE, border = "white", col = "#FCF8CD", width = 2, space = rep(2.5,nchr),las=1, xlim = c(0,3.0e+08), xaxt="n")
axis(1,c(0.0e+00,5.0e+07,1.0e+08,1.5e+08,2.0e+08,2.5e+08,3.0e+08),labels=c("0","50Mb","100Mb","150Mb","200Mb","250Mb","300Mb"))

##plot point
points(x=rep(1000000,nchr),y=seq(6,7*nchr,by=7),pch=20, col = "#7C1322", bg="#7C1322",cex = 3)

###Plot POE gene
for(i in 1:nrow(stage_gene_pos_sort_final)){
  ord <- which(as.numeric(rownames(chr_len_df))==stage_gene_pos_sort_final$Chrom_ori[i])
  y_0 <- 5.2+7*(ord-1)
  y_1<- 6.8+7*(ord-1)
  x_0 <- stage_gene_pos_sort_final$Start_ori[i]
  x_1 <- stage_gene_pos_sort_final$Start_ori[i]
  segments(x0=x_0,y0=y_0,x1=x_1,y1=y_1, col = "#2c7bb6",lwd = 0.5)
  
  x_point <- stage_gene_pos_sort_final$Start_alt[i]
  
  points(x=x_point,y=y_1+1,pch=20, col = "#2c7bb6", bg=NULL,cex = 1)
  
  segments(x0=x_0,y0=y_1+0.3,x1=x_point,y1=y_1+0.7, col = "#2c7bb6",lwd = 0.5)
}

###Draw imprinted genes database genes
for(i in 1:nrow(im_db_filter_sort_final)){
  ord <- which(as.numeric(rownames(chr_len_df))==im_db_filter_sort_final$Chrom_ori[i])
  y_0 <- 5.2+7*(ord-1)
  y_1<- 6.8+7*(ord-1)
  x_0 <- im_db_filter_sort_final$Start_ori[i]
  x_1 <- im_db_filter_sort_final$Start_ori[i]
  segments(x0=x_0,y0=y_0,x1=x_1,y1=y_1, col = "#fdae61",lwd = 0.5)
  
  x_point <- im_db_filter_sort_final$Start_alt[i]
  
  points(x=x_point,y=y_1+1,pch=20, col = "#fdae61", bg=NULL,cex = 1)
  
  segments(x0=x_0,y0=y_1+0.3,x1=x_point,y1=y_1+0.7, col = "#fdae61",lwd = 0.5)
}


##Drawing overlap gene
for(i in 1:nrow(overlap_pos)){
  ord <- which(as.numeric(rownames(chr_len_df))==overlap_pos$Chrom_ori[i])
  y_0 <- 5.2+7*(ord-1)
  y_1<- 6.8+7*(ord-1)
  x_0 <- overlap_pos$Start_ori[i]
  x_1 <- overlap_pos$Start_ori[i]
  segments(x0=x_0,y0=y_0,x1=x_1,y1=y_1, col = "#d7191c",lwd = 0.5)
  
  x_point <- overlap_pos$Start_alt[i]
  
  points(x=x_point,y=y_1+1,pch=20, col = "#d7191c", bg=NULL,cex = 1)
  
  segments(x0=x_0,y0=y_1+0.3,x1=x_point,y1=y_1+0.7, col = "#d7191c",lwd = 0.5)
}


##width=15,height=10