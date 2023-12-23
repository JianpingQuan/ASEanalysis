setwd("F:\\ATAC\\processFiles\\RNA-tissue_specific0711")
library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)

#Reads count based gene expression matrix
dataExp <- read.csv("F:\\RNA\\DEG_Tissues_Stages\\GeneExpressionAllSamples.csv", header = T, row.names = 1)

ref_tissue <- c("Br","Li","Mu","Pl")
ref_stage <- c("40","70","115","168")

for(j in 1:4){
  dataFilter <- select(dataExp, contains(ref_tissue[j]))
  tmp <- unlist(str_split(colnames(dataFilter),"_"))
  
  if(j >3){
    n = 1
  }else{
    n=3
  }
  for(i in 1:n){
  
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
    
    if(i==1){
      lrt_stage <- glmLRT(fit, coef=5)
      group_stage <- lrt_stage$table
      FDR_stage <- p.adjust(group_stage$PValue, method="fdr")
      group_stage <- cbind(group_stage, FDR_stage)
      write.csv(group_stage, paste0("Sta_F40_F70","Ref",ref_stage[i],"_Ref",ref_tissue[j],".csv"), quote = F, row.names = T)
    }else if (i==2){
      lrt_stage <- glmLRT(fit, coef=3)
      group_stage <- lrt_stage$table
      FDR_stage <- p.adjust(group_stage$PValue, method="fdr")
      group_stage <- cbind(group_stage, FDR_stage)
      write.csv(group_stage, paste0("Sta_F70_D1","Ref",ref_stage[i],"_Ref",ref_tissue[j],".csv"), quote = F, row.names = T)
    }else if(i==3){
      lrt_stage <- glmLRT(fit, coef=3)
      group_stage <- lrt_stage$table
      FDR_stage <- p.adjust(group_stage$PValue, method="fdr")
      group_stage <- cbind(group_stage, FDR_stage)
      write.csv(group_stage, paste0("Sta_D1_D168","Ref",ref_stage[i],"_Ref",ref_tissue[j],".csv"), quote = F, row.names = T)
    }
    
  }
}


##for placenta
j = 4

dataFilter <- select(dataExp, contains(ref_tissue[j]))
tmp <- unlist(str_split(colnames(dataFilter),"_"))

for(i in 1:1){
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
  
  if(i==1){
    lrt_stage <- glmLRT(fit, coef=3)
    group_stage <- lrt_stage$table
    FDR_stage <- p.adjust(group_stage$PValue, method="fdr")
    group_stage <- cbind(group_stage, FDR_stage)
    write.csv(group_stage, paste0("Sta_F40_F70","Ref",ref_stage[i],"_Ref",ref_tissue[j],".csv"), quote = F, row.names = T)
  }
}


##Brain Liver Muscle

F40_70 <- read.csv("Sta_F40_F70Ref40_RefMu.csv", header = T, row.names = 1)
  
F70_1 <- read.csv("Sta_F70_D1Ref70_RefMu.csv", header = T, row.names = 1)

D1_168 <- read.csv("Sta_D1_D168Ref115_RefMu.csv", header = T, row.names = 1)


F40_70_nosig <- F40_70 %>% filter(., FDR_stage>0.05)
F40_70_up <- F40_70 %>% filter(., FDR_stage<=0.05&logFC > 0)
F40_70_down <- F40_70 %>% filter(., FDR_stage<=0.05&logFC < 0)


F70_1_nosig <- F70_1 %>% filter(., FDR_stage>0.05)
F70_1_up <- F70_1 %>% filter(., FDR_stage<=0.05&logFC > 0)
F70_1_down <- F70_1 %>% filter(., FDR_stage<=0.05&logFC < 0)

D1_168_nosig <- D1_168 %>% filter(., FDR_stage>0.05)
D1_168_up <- D1_168 %>% filter(., FDR_stage<=0.05&logFC > 0)
D1_168_down <- D1_168 %>% filter(., FDR_stage<=0.05&logFC < 0)


#
(nnn <- intersect(intersect(rownames(F40_70_nosig),rownames(F70_1_nosig)),rownames(D1_168_nosig)) %>% length())
(nnu <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_nosig)), rownames(D1_168_up)) %>% length())
(nnd <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_nosig)), rownames(D1_168_down)) %>% length())

(nun <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_up)), rownames(D1_168_nosig)) %>% length())
(nuu <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_up)), rownames(D1_168_up)) %>% length())
(nud <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_up)), rownames(D1_168_down)) %>% length())

(ndn <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_down)), rownames(D1_168_nosig)) %>% length())
(ndu <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_down)), rownames(D1_168_up)) %>% length())
(ndd <- intersect(intersect(rownames(F40_70_nosig), rownames(F70_1_down)), rownames(D1_168_down)) %>% length())

#
(unn <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_nosig)), rownames(D1_168_nosig)) %>% length())
(unu <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_nosig)), rownames(D1_168_up)) %>% length())
(und <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_nosig)), rownames(D1_168_down)) %>% length())

(uun <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_up)), rownames(D1_168_nosig)) %>% length())
(uuu <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_up)), rownames(D1_168_up)) %>% length())
(uud <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_up)), rownames(D1_168_down)) %>% length())

(udn <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_down)), rownames(D1_168_nosig)) %>% length())
(udu <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_down)), rownames(D1_168_up)) %>% length())
(udd <- intersect(intersect(rownames(F40_70_up), rownames(F70_1_down)), rownames(D1_168_down)) %>% length())

#
(dnn <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_nosig)), rownames(D1_168_nosig)) %>% length())
(dnu <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_nosig)), rownames(D1_168_up)) %>% length())
(dnd <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_nosig)), rownames(D1_168_down)) %>% length())

(dun <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_up)), rownames(D1_168_nosig)) %>% length())
(duu <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_up)), rownames(D1_168_up)) %>% length())
(dud <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_up)), rownames(D1_168_down)) %>% length())

(ddn <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_down)), rownames(D1_168_nosig)) %>% length())
(ddu <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_down)), rownames(D1_168_up)) %>% length())
(ddd <- intersect(intersect(rownames(F40_70_down), rownames(F70_1_down)), rownames(D1_168_down)) %>% length())


##Placenta
F40_70 <- read.csv("Sta_F40_F70Ref40_RefPl.csv", header = T, row.names = 1)
F40_70_nosig <- F40_70 %>% filter(., FDR_stage>0.05)
F40_70_up <- F40_70 %>% filter(., FDR_stage<=0.05&logFC > 0)
F40_70_down <- F40_70 %>% filter(., FDR_stage<=0.05&logFC < 0)


dim(F40_70_nosig)
dim(F40_70_up)
dim(F40_70_down)
