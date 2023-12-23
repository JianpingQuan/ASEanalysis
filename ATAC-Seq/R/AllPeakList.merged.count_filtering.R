setwd("F:\\Manuscript\\ATAC\\data")
library(dplyr)
dat <- read.table("AllPeakList.merged.count.txt", header = T, row.names = 1)

datFilter <- select(dat, !contains(c("DS70_Li3","LS168_B_Br1","LS168_B_Li1","LS168_B_Mu1","LS168_C_Br2","LS168_C_Li2","LS168_C_Mu2","DS168_M_Br1","DS168_M_Li1","DS168_M_Mu1")))

#TPM (Transcripts Per Kilobase Million) 
counts2TPM <- function(count=count, efflength=efflen){   
  RPK <- count/(efflength/1000)     
  PMSC_rpk <- sum(RPK)/1e6           
  RPK/PMSC_rpk                       
}     

datTPM <- matrix(nrow = nrow(datFilter), ncol = 153)

for(i in 6:ncol(datFilter)){
  
  dat1 <- counts2TPM(count=datFilter[,i], efflength = datFilter$Length)
  
  datTPM[,i-5] <- dat1
  
}

datTPM <- as.data.frame(datTPM)
colnames(datTPM) <- colnames(datFilter[,6:158])
rownames(datTPM) <- rownames(datFilter)

##peaks filtering according reads count
datTPM$readcount <- rowMeans(datTPM)
datTPM_sort <- arrange(datTPM,desc(readcount)) %>% filter(readcount>0)

datTPM_sort_final <- filter(datTPM_sort, readcount>=200)

datTPM_filter <- filter(datTPM, !(rownames(datTPM) %in% rownames(datTPM_sort_final)))
datTPM_filter <- datTPM_filter[,-154] %>% select(.,!contains("LB70_Pl3"))

datFilter_filter <- filter(datFilter, rownames(datFilter) %in% rownames(datTPM_filter))

#TPM
write.table(datTPM_filter, "AllPeakList.merged.TPM2.txt",quote = F, sep = "\t",row.names = T)
#Read count
write.table(datFilter_filter, "AllPeakList.merged.ReadCount.txt",quote = F, sep = "\t",row.names = T)