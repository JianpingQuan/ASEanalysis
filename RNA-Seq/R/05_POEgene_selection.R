##POE gene identification
Input="F:\\RNA_20221101\\02_ModifiedGenome\\A_GeneExpressionMatrix\\"

library(dplyr)
library(edgeR)
library(ggplot2)

data_filter <- read.csv(paste0(Input,"AutosomeGeneCount_filtered2.csv"), header = T, row.names = 1)

###------------------------------------Function--------------------------------------------------
ASE_finder <- function(tissue,stage,fc,component,writeout){

	group <- select(data_filter,contains(tissue)) %>% select(contains(stage))
    
	y <- DGEList(counts=group,genes=rownames(group))
    
	isexpr <- rowSums(cpm(y)>1) >= 6
    
	y <- y[isexpr, , keep.lib.sizes=FALSE]
    
	y <- calcNormFactors(y,method = "TMM")
			
	po  <- rep(0,ncol(group))
	ag  <- rep(0,ncol(group))
	mg  <- rep(0,ncol(group))
	sex <- rep(0,ncol(group))
    
	po_ord  <- grep("Mat",colnames(group))
	ag_ord  <- grep("_D" ,colnames(group))
	mg_ord  <- grep("DL" ,colnames(group))
	sex_ord <- grep("_F_",colnames(group))
    
    po[po_ord]   <- 1 #set Maternal Genome to 1
    ag[ag_ord]   <- 1 #set Duroc allele to 1
    mg[mg_ord]   <- 1 #set the offspring fo Duroc mother to 1
	sex[sex_ord] <- 1 #set female to 1
    
    PO  <- po
    AG  <- ag
    MG  <- mg
	Sex <- sex
	
	data.frame(Sample=colnames(y), AG, PO, MG, Sex)
	design <- model.matrix(~PO+AG+MG+Sex)
	y <- estimateDisp(y, design, robust=TRUE)
	y$common.dispersion
	fit <- glmFit(y, design)
	
	##for PO component
	lrt_PO <- glmLRT(fit, coef=2)
	group_PO <- lrt_PO$table
	FDR <- p.adjust(group_PO$PValue, method="fdr")
	group_PO <- cbind(group_PO, FDR)
	group_PO_final <- filter(group_PO,FDR<0.05)%>%filter(abs(logFC)>=fc)
    
	if(component=="PO"&&writeout){
	write.csv(group_PO, paste0(tissue,"_",stage,"_AllGenes","PO2",".csv"), quote = F,row.names = T)
	write.csv(group_PO_final, paste0(tissue,"_",stage,"_",fc,"PO_removed2",".csv"), quote = F,row.names = T)
	}
	
    ##for AG component
	lrt_AG <- glmLRT(fit, coef=3)
	group_AG <- lrt_AG$table
	FDR <- p.adjust(group_AG$PValue, method="fdr")
	group_AG <- cbind(group_AG, FDR)
	group_AG_final <- filter(group_AG,FDR<0.05)%>%filter(abs(logFC)>=fc)
	
	if(component=="AG"&&writeout){
	write.csv(group_AG, paste0(tissue,"_",stage,"_AllGenes","AG2",".csv"), quote = F,row.names = T)
	write.csv(group_AG_final, paste0(tissue,"_",stage,"_",fc,"AG_removed2",".csv"), quote = F,row.names = T)
	}
	
    ##For MG component
	lrt_MG <- glmLRT(fit, coef=4)
	group_MG <- lrt_MG$table
	FDR <- p.adjust(group_MG$PValue, method="fdr")
	group_MG <- cbind(group_MG, FDR)
	group_MG_final <- filter(group_MG,FDR<0.05)%>%filter(abs(logFC)>=fc)
    if(component=="MG"&&writeout){
	write.csv(group_MG, paste0(tissue,"_",stage,"_AllGenes","MG2",".csv"), quote = F,row.names = T)
	write.csv(group_MG_final, paste0(tissue,"_",stage,"_",fc,"MG_removed2",".csv"), quote = F,row.names = T)
	}
	
    ##For sex component
	lrt_Sex <- glmLRT(fit, coef=5)
	group_Sex <- lrt_Sex$table
	FDR <- p.adjust(group_Sex$PValue, method="fdr")
	group_Sex <- cbind(group_Sex, FDR)
	group_Sex_final <- filter(group_Sex,FDR<0.05)%>%filter(abs(logFC)>=fc)
	
    if(component=="Sex"&&writeout){
	write.csv(group_Sex, paste0(tissue,"_",stage,"_AllGenes","Sex2",".csv"), quote = F,row.names = T)
	write.csv(group_Sex_final, paste0(tissue,"_",stage,"_",fc,"Sex_removed2",".csv"), quote = F,row.names = T)
	}
	
	
	##Component selection
	if(component=="PO"){
	###filtering the PO effect ASE
	duroc <- group %>% dplyr::select(starts_with("DL")) %>% dplyr::filter(rownames(group) %in% rownames(group_PO_final))
	lulai <- group %>% dplyr::select(starts_with("LD")) %>% dplyr::filter(rownames(group) %in% rownames(group_PO_final))
  
	gene_pre <- c()
	lulai_pat <-  lulai %>% select(contains("LD")) %>% select(contains("Pat"))
	lulai_mat <-  lulai %>% select(contains("LD")) %>% select(contains("Mat"))
  
	duroc_pat <- duroc %>% select(contains("DL")) %>% select(contains("Pat"))
	duroc_mat <- duroc %>% select(contains("DL")) %>% select(contains("Mat"))
  
	lulai_ration <- round((rowMeans(lulai_pat))/(rowMeans(lulai_mat)),4)
	duroc_ration <- round((rowMeans(duroc_pat))/(rowMeans(duroc_mat)),4)
	
	##
	for(k in 1:length(lulai_ration)){
		if((lulai_ration[k]>1&&duroc_ration[k]>1)||(duroc_ration[k]<1&&lulai_ration[k]<1)){
			gene_pre <- append(gene_pre,rownames(duroc_pat)[k])
		}
	}
  
	group_2 <- group %>% dplyr::filter(rownames(group) %in% gene_pre)
  
	pat_po <- group_2 %>% dplyr::select(contains("Pat"))
	mat_po <- group_2 %>% dplyr::select(contains("Mat"))
  
	gene_final <- c()
  
	for(m in 1:nrow(group_2)){
		if(sum((pat_po[m,]-mat_po[m,])>0)>= ncol(pat_po)-1||sum((pat_po[m,]-mat_po[m,])<0)>= ncol(pat_po)-1){
			gene_final <- append(gene_final,rownames(group_2)[m])
		}
	}
	return(assign(paste0(tissue,stage,"PO"),gene_final))
	
	}else if(component=="AG"){
		duroc <- group %>% dplyr::select(starts_with("DL")) %>% dplyr::filter(rownames(group) %in% rownames(group_AG_final))
		lulai <- group %>% dplyr::select(starts_with("LD")) %>% dplyr::filter(rownames(group) %in% rownames(group_AG_final))
  
		gene_pre <- c()
		lulai_pat <-  lulai %>% select(contains("LD")) %>% select(contains("Pat"))
		lulai_mat <-  lulai %>% select(contains("LD")) %>% select(contains("Mat"))
  
		duroc_pat <- duroc %>% select(contains("DL")) %>% select(contains("Pat"))
		duroc_mat <- duroc %>% select(contains("DL")) %>% select(contains("Mat"))
		
		Duroc_LD <- round((rowMeans(lulai_pat))/(rowMeans(lulai_mat)),4)
		Duroc_DL <- round((rowMeans(duroc_mat))/(rowMeans(duroc_pat)),4)
		
		for(k in 1:length(Duroc_LD)){
		if((Duroc_LD[k]>1&&Duroc_DL[k]>1)||(Duroc_LD[k]<1&&Duroc_DL[k]<1)){
			gene_pre <- append(gene_pre,rownames(duroc_pat)[k])
			}
		}
		group_2 <- group %>% dplyr::filter(rownames(group) %in% gene_pre)
  
		pat_ag <- group_2 %>% dplyr::select(contains("Pat"))
		mat_ag <- group_2 %>% dplyr::select(contains("Mat"))
		
  		duroc <- cbind(select(pat_ag,contains("LD")),select(mat_ag,contains("DL")))
		lulai <- cbind(select(mat_ag,contains("LD")),select(pat_ag,contains("DL")))
		
		gene_final <- c()
  
		for(m in 1:nrow(group_2)){
			if(sum((duroc[m,]-lulai[m,])>0) >= ncol(pat_ag)-1||sum((duroc[m,]-lulai[m,])<0) >= ncol(pat_ag)-1){
				gene_final <- append(gene_final,rownames(group_2)[m])
			}
		}
		return(assign(paste0(tissue,stage,"AG"),gene_final))
	}else if(component=="Sex"){
		duroc <- group %>% dplyr::select(starts_with("DL")) %>% dplyr::filter(rownames(group) %in% rownames(group_Sex_final))
		lulai <- group %>% dplyr::select(starts_with("LD")) %>% dplyr::filter(rownames(group) %in% rownames(group_Sex_final))
  
		gene_final <- c()
		lulai_F <-  lulai %>% select(contains("LD")) %>% select(contains("F"))
		lulai_M <-  lulai %>% select(contains("LD")) %>% select(contains("M"))
  
		duroc_F <- duroc %>% select(contains("DL")) %>% select(contains("F"))
		duroc_M <- duroc %>% select(contains("DL")) %>% select(contains("M"))
		
		lulai_sex <- round((rowMeans(lulai_F))/(rowMeans(lulai_M)),4)
		duroc_sex <- round((rowMeans(duroc_F))/(rowMeans(duroc_M)),4)
		
		for(k in 1:length(duroc_F)){
		if((lulai_sex[k]>1&&duroc_sex[k]>1)||(lulai_sex[k]<1&&duroc_sex[k]<1)){
			gene_final <- append(gene_final,rownames(duroc_sex)[k])
			}
		}
		return(assign(paste0(tissue,stage,"Sex"),gene_final))
		
	}else if(component=="MG"){
		gene_final <- rownames(group_MG_final)}
		return(assign(paste0(tissue,stage,"MG"),gene_final))
}



##------------------------------------------Function using--------------------------------------
##----------Generate N-1FC2_filter.PO.csv--------------
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result")
stage <- c("40","70","115","168")
tissue <- c("Br", "Li", "Mu")

for(i in 1:4){
  for(j in 1:3){
    st <- stage[i]
    ti <- tissue[j]
    assign(paste0(ti,st), ASE_finder(tissue=ti,stage=st,fc=2,component="PO",writeout=TRUE))
  }
}

Pl40 <- ASE_finder(tissue="Pl",stage="40",fc=2,component="PO",writeout=TRUE)
Pl70 <- ASE_finder(tissue="Pl",stage="70",fc=2,component="PO",writeout=TRUE)

sq <- seq(max(c(length(Br40),length(Br70),length(Br115),length(Br168),length(Li40),length(Li70),length(Li115),length(Li168),length(Mu40),length(Mu70),length(Mu115),length(Mu168),length(Pl40),length(Pl70))))

all <- data.frame(Br_40=Br40[sq],Br_70=Br70[sq],Br_115=Br115[sq],Br_168=Br168[sq],Li_40=Li40[sq],Li_70=Li70[sq],Li_115=Li115[sq],Li_168=Li168[sq],Mu_40=Mu40[sq],Mu_70=Mu70[sq],Mu_115=Mu115[sq],Mu_168=Mu168[sq],Pl_40=Pl40[sq],Pl_70=Pl70[sq])

write.csv(all, "ASE_N-1FC2_filter.PO2.csv", quote=F)



##------------------------------------------add gene with loose standard-------------------------------------
#--------Generate FC2_filter_final.PO.csv-----------
###use a looser threshold value to find all possible genes
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result")
library(dplyr)
tissue <- c("Br", "Li", "Mu","Pl")
stage <- c(40,70,115,168)

dat <- read.csv("ASE_N-1FC2_filter.PO2.csv", header = T, row.names = 1)

needremove <- read.csv("poe_age_both.csv",header=T,row.names=1) %>% filter(.,type=="AGE")

dat <- apply(dat, 2, function(x) ifelse(x %in% needremove$gene, NA, x))

gene <- c()
for(i in 1:ncol(dat)){
  gene <- append(gene, dat[,i])
}
all_gene_unique <- unique(na.omit(gene))
a <- as.data.frame(all_gene_unique)

write.csv(a, "POE.list2.csv", quote = F)


##Brain,Liver,Muscle
for(i in 1:3){
  for(j in 1:4){
    Group <- read.csv(paste0(tissue[i], "_", stage[j], "_AllGenesPO.csv"), header = T, row.names = 1)
	
    col_name <- paste0(tissue[i], "_", stage[j])
    index <- col_name == colnames(dat)
    int_gene <- filter(Group, rownames(Group) %in% all_gene_unique)
    int_gene_filter <- filter(int_gene, FDR < 0.01)
    gene_final <- union(dat[,index],rownames(int_gene_filter))
    assign(paste0(tissue[i],stage[j]), gene_final)
  }
}

##placenta
for(i in 4:4){
  for(j in 1:2){
    Group <- read.csv(paste0(tissue[i], "_", stage[j], "_AllGenesPO.csv"), header = T, row.names = 1)
    col_name <- paste0(tissue[i], "_", stage[j])
    index <- col_name == colnames(dat)
    int_gene <- filter(Group, rownames(Group) %in% all_gene_unique)
    int_gene_filter <- filter(int_gene, FDR < 0.01)
    gene_final <- union(dat[,index],rownames(int_gene_filter))
    assign(paste0(tissue[i],stage[j]), gene_final)
  }
}

sq <- seq(max(c(length(Br40),length(Br70),length(Br115),length(Br168),length(Li40),length(Li70),length(Li115),length(Li168),length(Mu40),length(Mu70),length(Mu115),length(Mu168),length(Pl40),length(Pl70))))

all_final <- data.frame(Br_40=Br40[sq],Br_70=Br70[sq],Br_115=Br115[sq],Br_168=Br168[sq],Li_40=Li40[sq],Li_70=Li70[sq],Li_115=Li115[sq],Li_168=Li168[sq],Mu_40=Mu40[sq],Mu_70=Mu70[sq],Mu_115=Mu115[sq],Mu_168=Mu168[sq],Pl_40=Pl40[sq],Pl_70=Pl70[sq])

write.csv(all_final, "ASE_N-1_filter_final.PO2.csv",quote=F)




##------------Output POE genes expreesion origin-----------
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result")
library(dplyr)
dat <- read.csv("ASE_N-1_filter_final.PO2.csv", header = T, row.names = 1)
tissue <- c("Br","Li","Mu","Pl")
stage <- c("40","70","115","168")

for(i in 1:4){
  if(i==4){
    n=2
  }else{
    n=4
  }
  for(j in 1:n){
    raw <- read.csv(paste0(tissue[i], "_", stage[j], "_AllGenesPO2.csv"), header = T, row.names = 1)
    col_name <- paste0(tissue[i], "_", stage[j])
    index <- col_name == colnames(dat)
    gene_unique <- na.omit(dat[,index])
    gene_FC <- filter(raw, rownames(raw) %in% gene_unique)
    bias <- c()
    for(k in 1:nrow(gene_FC)){
      if(gene_FC[k,1]< 0){
        bias <- append(bias,"Pat")
      }else if(gene_FC[k,1] >0){
        bias <- append(bias,"Mat")
      }
    }
    gene_FC$bias <- bias
    gene_FC$name <- rownames(gene_FC)
    assign(paste0(tissue[i], "_", stage[j]), gene_FC)
    write.csv(gene_FC, paste0(tissue[i], "_", stage[j],"_bias2.csv"), quote = F)
  }
}

