setwd("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\03_POE基因热图展示\\data")
library(dplyr)
library(ggplot2)
library(basicTrendline)
library(customLayout)
library(magrittr)
library(pheatmap)

data_filter <- read.csv("AutosomeGeneCount_filtered.csv", header = T, row.names = 1)

POE <- read.table("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\05_已知POE基因reads数量探索\\ImprintingGeneList.txt")

data_filterPOE <- data_filter %>% filter(rownames(.) %in% POE$V1)

##computed Ratio
##Brain,Liver,Muscle

tissue <- c("Br","Li","Mu","Pl")
stage <- c("40","70","115","168")
ref <- c("Pat","Mat")

for(i in 1:3){
  for(j in 1:4){
    for(k in 1:2){
      group <- select(data_filterPOE,contains(tissue[i])) %>% select(contains(stage[j])) %>% select(contains(ref[k]))
      rowmean <- rowMeans(group)
      assign(paste0(tissue[i],stage[j],ref[k]),rowmean)
    }
  }
}

##Placenta
for(i in 4:4){
  for(j in 1:2){
    for(k in 1:2){
      group <- select(data_filterPOE,contains(tissue[i])) %>% select(contains(stage[j])) %>% select(contains(ref[k]))
      rowmean <- rowMeans(group)
      assign(paste0(tissue[i],stage[j],ref[k]),rowmean)
    }
  }
}


Mean <- cbind(Br40Mat,Br40Pat,Li40Mat,Li40Pat,Mu40Mat,Mu40Pat,Pl40Mat,Pl40Pat,Br70Mat,Br70Pat,Li70Mat,Li70Pat,Mu70Mat,Mu70Pat,Pl70Mat,Pl70Pat,Br115Mat,Br115Pat,Li115Mat,Li115Pat,Mu115Mat,Mu115Pat,Br168Mat,Br168Pat,Li168Mat,Li168Pat,Mu168Mat,Mu168Pat) %>% as.data.frame()

colnames(Mean) <- c("Br40Mat","Br40Pat","Li40Mat","Li40Pat","Mu40Mat","Mu40Pat","Pl40Mat","Pl40Pat","Br70Mat","Br70Pat","Li70Mat","Li70Pat","Mu70Mat","Mu70Pat","Pl70Mat","Pl70Pat","Br115Mat","Br115Pat","Li115Mat","Li115Pat","Mu115Mat","Mu115Pat","Br168Mat","Br168Pat","Li168Mat","Li168Pat","Mu168Mat","Mu168Pat")

for(i in seq(1,14,2)){
	for(j in 1:nrow(Mean)){
		if(Mean[j,i]<=2&&Mean[j,i+1]<=2){
			Mean[j,i] <- 2
			Mean[j,i+1] <- 2
		}
	}
}

write.csv(Mean, "GroupAbundanceMean.csv", quote = F)

Mean <- read.table("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\figures_tables\\05_已知POE基因reads数量探索\\Mean.txt", header = T, row.names = 1)

#MeanRatio <- cbind(Br40Mat/(Br40Pat+Br40Mat),Li40Mat/(Li40Pat+Li40Mat),Mu40Mat/(Mu40Pat+Mu40Mat),Pl40Mat/(Pl40Pat+Pl40Mat),Br70Mat/(Br70Pat+Br70Mat),Li70Mat/(Li70Pat+Li70Mat),Mu70Mat/(Mu70Pat+Mu70Mat),Pl70Mat/(Pl70Pat+Pl70Mat),Br115Mat/(Br115Pat+Br115Mat),Li115Mat/(Li115Pat+Li115Mat),Mu115Mat/(Mu115Pat+Mu115Mat),Br168Mat/(Br168Pat+Br168Mat),Li168Mat/(Li168Pat+Li168Mat),Mu168Mat/(Mu168Pat+Mu168Mat)) %>% as.data.frame()


#MeanRatio <- cbind(Br40Mat/Br40Pat,Li40Mat/Li40Pat,Mu40Mat/Mu40Pat,Pl40Mat/Pl40Pat,Br70Mat/Br70Pat,Li70Mat/Li70Pat,Mu70Mat/Mu70Pat,Pl70Mat/Pl70Pat,Br115Mat/Br115Pat,Li115Mat/Li115Pat,Mu115Mat/Mu115Pat,Br168Mat/Br168Pat,Li168Mat/Li168Pat,Mu168Mat/Mu168Pat) %>% as.data.frame()


MeanRatio <- cbind(Mean$Br40Mat/(Mean$Br40Pat+Mean$Br40Mat),Mean$Li40Mat/(Mean$Li40Pat+Mean$Li40Mat),Mean$Mu40Mat/(Mean$Mu40Pat+Mean$Mu40Mat),Mean$Pl40Mat/(Mean$Pl40Pat+Mean$Pl40Mat),Mean$Br70Mat/(Mean$Br70Pat+Mean$Br70Mat),Mean$Li70Mat/(Mean$Li70Pat+Mean$Li70Mat),Mean$Mu70Mat/(Mean$Mu70Pat+Mean$Mu70Mat),Mean$Pl70Mat/(Mean$Pl70Pat+Mean$Pl70Mat),Mean$Br115Mat/(Mean$Br115Pat+Mean$Br115Mat),Mean$Li115Mat/(Mean$Li115Pat+Mean$Li115Mat),Mean$Mu115Mat/(Mean$Mu115Pat+Mean$Mu115Mat),Mean$Br168Mat/(Mean$Br168Pat+Mean$Br168Mat),Mean$Li168Mat/(Mean$Li168Pat+Mean$Li168Mat),Mean$Mu168Mat/(Mean$Mu168Pat+Mean$Mu168Mat)) %>% as.data.frame()


colnames(MeanRatio) <- c("Br40","Li40","Mu40","Pl40","Br70","Li70","Mu70","Pl70","Br115","Li115","Mu115","Br168","Li168","Mu168")

rownames(MeanRatio) <- rownames(Mean)


cc1 = grDevices::colorRampPalette(c("#7C1322","#F8E5C6","#143C7B") %>% ggplot2::alpha(0.99))

p1 <- pheatmap::pheatmap(MeanRatio,
                         cluster_rows = F,
                         cluster_cols=F,
                         #gaps_col = c(6,12,18,24,30,36,41),
                         cellwidth=10,
                         cellheight=10,
                         color = cc1(30),
                         na_col="black",
                         show_rownames = T,
                         show_colnames = T)