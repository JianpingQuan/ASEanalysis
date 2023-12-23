
##POE dominate expression identification
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result")
library(dplyr)

POE <- read.csv("ASE_N-1_filter_final.PO2.csv", header = T, row.names = 1)

Br_40 <- POE$Br_40[!is.na(POE$Br_40)]
Br_70 <- POE$Br_70[!is.na(POE$Br_70)]
Br_115 <- POE$Br_115[!is.na(POE$Br_115)]
Br_168 <- POE$Br_168[!is.na(POE$Br_168)]

Li_40 <- POE$Li_40[!is.na(POE$Li_40)]
Li_70 <- POE$Li_70[!is.na(POE$Li_70)]
Li_115 <- POE$Li_115[!is.na(POE$Li_115)]
Li_168 <- POE$Li_168[!is.na(POE$Li_168)]

Mu_40 <- POE$Mu_40[!is.na(POE$Mu_40)]
Mu_70 <- POE$Mu_70[!is.na(POE$Mu_70)]
Mu_115 <- POE$Mu_115[!is.na(POE$Mu_115)]
Mu_168 <- POE$Mu_168[!is.na(POE$Mu_168)]

Pl_40 <- POE$Pl_40[!is.na(POE$Pl_40)]
Pl_70 <- POE$Pl_70[!is.na(POE$Pl_70)]


Br40_filter <- read.csv("Br_40_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_40)
Br70_filter <- read.csv("Br_70_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_70)
Br115_filter <- read.csv("Br_115_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_115)
Br168_filter <- read.csv("Br_168_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_168)

Li40_filter <- read.csv("Li_40_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_40)
Li70_filter <- read.csv("Li_70_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_70)
Li115_filter <- read.csv("Li_115_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_115)
Li168_filter <- read.csv("Li_168_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_168)

Mu40_filter <- read.csv("Mu_40_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_40)
Mu70_filter <- read.csv("Mu_70_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_70)
Mu115_filter <- read.csv("Mu_115_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_115)
Mu168_filter <- read.csv("Mu_168_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_168)

Pl40_filter <- read.csv("Pl_40_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Pl_40)
Pl70_filter <- read.csv("Pl_70_AllGenesPO2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Pl_70)



Br40_filter$Group <- rep("F40Brain", nrow(Br40_filter))
Br70_filter$Group <- rep("F70Brain", nrow(Br70_filter))
Br115_filter$Group <- rep("D1Brain", nrow(Br115_filter))
Br168_filter$Group <- rep("D168Brain", nrow(Br168_filter))

Li40_filter$Group <- rep("F40Liver", nrow(Li40_filter))
Li70_filter$Group <- rep("F70Liver", nrow(Li70_filter))
Li115_filter$Group <- rep("D1Liver", nrow(Li115_filter))
Li168_filter$Group <- rep("D168Liver", nrow(Li168_filter))

Mu40_filter$Group <- rep("F40Muscle", nrow(Mu40_filter))
Mu70_filter$Group <- rep("F70Muscle", nrow(Mu70_filter))
Mu115_filter$Group <- rep("D1Muscle", nrow(Mu115_filter))
Mu168_filter$Group <- rep("D168Muscle", nrow(Mu168_filter))

Pl40_filter$Group <- rep("F40Placenta", nrow(Pl40_filter))
Pl70_filter$Group <- rep("F70Placenta", nrow(Pl70_filter))

##--------------------------------------------------------
Br40_filter$Gene <- rownames(Br40_filter)
Br70_filter$Gene <- rownames(Br70_filter)
Br115_filter$Gene <- rownames(Br115_filter)
Br168_filter$Gene <- rownames(Br168_filter)

Li40_filter$Gene <- rownames(Li40_filter)
Li70_filter$Gene <- rownames(Li70_filter)
Li115_filter$Gene <- rownames(Li115_filter)
Li168_filter$Gene <- rownames(Li168_filter)

Mu40_filter$Gene <- rownames(Mu40_filter)
Mu70_filter$Gene <- rownames(Mu70_filter)
Mu115_filter$Gene <- rownames(Mu115_filter)
Mu168_filter$Gene <- rownames(Mu168_filter)

Pl40_filter$Gene <- rownames(Pl40_filter)
Pl70_filter$Gene <- rownames(Pl70_filter)


POE_rowbind <- rbind(Br40_filter,Br70_filter,Br115_filter,Br168_filter, Li40_filter,Li70_filter,Li115_filter,Li168_filter, Mu40_filter,Mu70_filter, Mu115_filter, Mu168_filter, Pl40_filter, Pl70_filter)

Dir <- c()

for(i in 1:nrow(POE_rowbind)){
  
  if(POE_rowbind[i,1]>0){
    
    Dir[i] <- c("Maternal")
  }else{
    Dir[i] <- c("Paternal")
  }
}


POE_rowbind$DominantAllele <- Dir

rownames(POE_rowbind) <- seq(1,nrow(POE_rowbind))


##
Pat <- POE_rowbind %>% filter(.,DominantAllele == "Paternal")
Mat <- POE_rowbind %>% filter(.,DominantAllele == "Maternal")

a <- intersect(Pat$Gene, Mat$Gene)

Reverse <- POE_rowbind %>% filter(.,Gene %in% a) %>% arrange(., Gene)


write.csv(Pat, "Paternal_dominate_expression_POE.csv", quote = F, row.names = F)
write.csv(Mat, "Maternal_dominate_expression_POE.csv", quote = F, row.names = F)
write.csv(Reverse, "Reverse_dominate_expression_POE.csv", quote = F, row.names = F)

write.csv(POE_rowbind, "POE_dominate_expression.csv", quote = F, row.names = F)


###########################################################################################
##AGE diminate expression identification
setwd("F:\\RNA_20221101\\02_ModifiedGenome\\D_AGEgeneSelection\\result")
library(dplyr)

AGE <- read.csv("ASE_N-1_filter_final.AG2.csv", header = T, row.names = 1)

Br_40 <- AGE$Br_40[!is.na(AGE$Br_40)]
Br_70 <- AGE$Br_70[!is.na(AGE$Br_70)]
Br_115 <- AGE$Br_115[!is.na(AGE$Br_115)]
Br_168 <- AGE$Br_168[!is.na(AGE$Br_168)]

Li_40 <- AGE$Li_40[!is.na(AGE$Li_40)]
Li_70 <- AGE$Li_70[!is.na(AGE$Li_70)]
Li_115 <- AGE$Li_115[!is.na(AGE$Li_115)]
Li_168 <- AGE$Li_168[!is.na(AGE$Li_168)]

Mu_40 <- AGE$Mu_40[!is.na(AGE$Mu_40)]
Mu_70 <- AGE$Mu_70[!is.na(AGE$Mu_70)]
Mu_115 <- AGE$Mu_115[!is.na(AGE$Mu_115)]
Mu_168 <- AGE$Mu_168[!is.na(AGE$Mu_168)]

Pl_40 <- AGE$Pl_40[!is.na(AGE$Pl_40)]
Pl_70 <- AGE$Pl_70[!is.na(AGE$Pl_70)]


Br40_filter <- read.csv("Br_40_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_40)
Br70_filter <- read.csv("Br_70_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_70)
Br115_filter <- read.csv("Br_115_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_115)
Br168_filter <- read.csv("Br_168_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Br_168)

Li40_filter <- read.csv("Li_40_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_40)
Li70_filter <- read.csv("Li_70_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_70)
Li115_filter <- read.csv("Li_115_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_115)
Li168_filter <- read.csv("Li_168_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Li_168)

Mu40_filter <- read.csv("Mu_40_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_40)
Mu70_filter <- read.csv("Mu_70_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_70)
Mu115_filter <- read.csv("Mu_115_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_115)
Mu168_filter <- read.csv("Mu_168_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Mu_168)

Pl40_filter <- read.csv("Pl_40_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Pl_40)
Pl70_filter <- read.csv("Pl_70_AllGenesAG2.csv", header = T, row.names = 1) %>% filter(., rownames(.) %in% Pl_70)



Br40_filter$Group <- rep("F40Brain", nrow(Br40_filter))
Br70_filter$Group <- rep("F70Brain", nrow(Br70_filter))
Br115_filter$Group <- rep("D1Brain", nrow(Br115_filter))
Br168_filter$Group <- rep("D168Brain", nrow(Br168_filter))

Li40_filter$Group <- rep("F40Liver", nrow(Li40_filter))
Li70_filter$Group <- rep("F70Liver", nrow(Li70_filter))
Li115_filter$Group <- rep("D1Liver", nrow(Li115_filter))
Li168_filter$Group <- rep("D168Liver", nrow(Li168_filter))

Mu40_filter$Group <- rep("F40Muscle", nrow(Mu40_filter))
Mu70_filter$Group <- rep("F70Muscle", nrow(Mu70_filter))
Mu115_filter$Group <- rep("D1Muscle", nrow(Mu115_filter))
Mu168_filter$Group <- rep("D168Muscle", nrow(Mu168_filter))

Pl40_filter$Group <- rep("F40Placenta", nrow(Pl40_filter))
Pl70_filter$Group <- rep("F70Placenta", nrow(Pl70_filter))

##--------------------------------------------------------
Br40_filter$Gene <- rownames(Br40_filter)
Br70_filter$Gene <- rownames(Br70_filter)
Br115_filter$Gene <- rownames(Br115_filter)
Br168_filter$Gene <- rownames(Br168_filter)

Li40_filter$Gene <- rownames(Li40_filter)
Li70_filter$Gene <- rownames(Li70_filter)
Li115_filter$Gene <- rownames(Li115_filter)
Li168_filter$Gene <- rownames(Li168_filter)

Mu40_filter$Gene <- rownames(Mu40_filter)
Mu70_filter$Gene <- rownames(Mu70_filter)
Mu115_filter$Gene <- rownames(Mu115_filter)
Mu168_filter$Gene <- rownames(Mu168_filter)

Pl40_filter$Gene <- rownames(Pl40_filter)
Pl70_filter$Gene <- rownames(Pl70_filter)


AGE_rowbind <- rbind(Br40_filter,Br70_filter,Br115_filter,Br168_filter, Li40_filter,Li70_filter,Li115_filter,Li168_filter, Mu40_filter,Mu70_filter, Mu115_filter, Mu168_filter, Pl40_filter, Pl70_filter)

Dir <- c()

for(i in 1:nrow(AGE_rowbind)){
  
  if(AGE_rowbind[i,1]>0){
    
    Dir[i] <- c("Duroc")
  }else{
    Dir[i] <- c("Lulai")
  }
}


AGE_rowbind$DominantAllele <- Dir

rownames(AGE_rowbind) <- seq(1,nrow(AGE_rowbind))


##
Duroc <- AGE_rowbind %>% filter(.,DominantAllele == "Duroc")
Lulai <- AGE_rowbind %>% filter(.,DominantAllele == "Lulai")

a <- intersect(Duroc$Gene, Lulai$Gene)

Reverse <- AGE_rowbind %>% filter(.,Gene %in% a) %>% arrange(., Gene)


write.csv(Duroc, "Duroc_dominate_expression_AGE.csv", quote = F, row.names = F)
write.csv(Lulai, "Lulai_dominate_expression_AGE.csv", quote = F, row.names = F)
write.csv(Reverse, "Reverse_dominate_expression_AGE.csv", quote = F, row.names = F)

write.csv(AGE_rowbind, "AGE_dominate_expression.csv", quote = F, row.names = F)



