InputPath="F:\\RNA_20221101\\02_ModifiedGenome\\A_GeneExpressionMatrix\\"
OutputPath="F:\\RNA_20221101\\02_ModifiedGenome\\B_SampleEvaluation\\a_MappingBias\\Figure\\"

library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggsci)

RawReadCount <-read.table(paste(InputPath,"AutosomeGeneCountFinal.matrix.txt",sep=""), header=F, row.names = 1) 
SampleParent <- read.table(paste(InputPath, "Sample_Parent.list.txt", sep=""), header = F, row.names = 1) 

colnames(RawReadCount) <- rownames(SampleParent)

RC.colSum <- colSums(RawReadCount)

dat <- data.frame(Sample = rownames(SampleParent), SumCount = RC.colSum)

dat$Parent <- as.factor(str_sub(dat$Sample,-3)) 

dat$Sample <- str_sub(dat$Sample,1,-5)

Mat <- filter(dat, Parent=="Mat")

Pat <- filter(dat, Parent=="Pat")

Ratio <- Mat$SumCount/Pat$SumCount

datPlot1 <- data.frame(Sample=Mat$Sample, Mat=Mat$SumCount, Pat=Pat$SumCount, Ratio=Ratio)

datPlot1$Tissue <- str_sub(datPlot1$Sample,-3,-2)

write.csv(datPlot1, paste0(OutputPath,"Pat_Mat_totalReadsCount.csv"), quote = F, row.names = F)

datPlot2 <- dat

#bar plot
p <- ggplot(datPlot2, aes(Sample,SumCount, fill= Parent))

p1 <- p + geom_bar(stat='identity',position='fill',width = 0.8) + scale_fill_manual(values=c("#7C1322","#143C7B")) + labs(title = "The count of mapped reads on autosome between two personalize genomes")+xlab(NULL) + ylab("Read Count Percentage")

p2 <- p1 + theme_bw() + theme(plot.title = element_text(size = 15, colour = "black",hjust=0.5), axis.text.x = element_text(angle = 90, vjust=0.3, size = 5), panel.grid.major = element_line(color = "black", size = 0.2), legend.position ="top")

p3 <-  p2 + scale_y_continuous(expand = c(0,0))

p4 <- p3 + geom_hline(aes(yintercept=0.60), color = "grey",lty=2) + geom_hline(aes(yintercept=0.40), color = "grey",lty=2)+geom_hline(aes(yintercept=0.50), color = "purple")

pdf(paste0(OutputPath,"Paternal_Maternal.stackbar.pdf"),width = 15, height = 3)
p4
dev.off()


