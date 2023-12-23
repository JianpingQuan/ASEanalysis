library(ggplot2)
library(dplyr)
library(stringr)

stage <- c("40","70","115","168")
tissue <- c("Br", "Li", "Mu", "Pl")

##set color
color <- c("#8DD2C6", "#FFFFB3", "#BDB9DA", "#FA8071", "#80B0D2", "#FCB461", "#B3DE69", "#FBCDE4", "#D9D9D9", "#BB80BC", "#CCEBC5", "#FFEC6E")

###Br Li,Mu
for(i in 1:4){
  for(j in 1:4){
    setwd(paste0("F:\\ATAC\\ATAC-seqDatafinal_plot\\Insert size\\",stage[j],"_",tissue[i]))
    
    files <- list.files(paste0("F:\\ATAC\\ATAC-seqDatafinal_plot\\Insert size\\",stage[j],"_", tissue[i]), pattern = "*.txt") 
    samples <- str_remove(files, ".InsertSizes.txt")
    
    data <- lapply(list.files(paste0("F:\\ATAC\\ATAC-seqDatafinal_plot\\Insert size\\",stage[j],"_", tissue[i]), pattern = "*.txt"), read.table, header = T)
    
    ##filtering
    data <- lapply(data, function(df) {df %>% filter(insert_size <= 1000)})
    names(data) <- samples
    
    p <- ggplot()
    
    ##add layer
    for(n in 1:length(samples)){
      p <- p + geom_line(data=data[[n]],aes(x=insert_size, y=All_Reads.fr_count+All_Reads.rf_count),color = color[n])
    }
    
    p1 <- p+ theme_bw() + scale_x_continuous(breaks = seq(0,1000,200)) + labs(x= "The length of insert size", y="Reads count")
    
    ###add sample name
    for(m in 1:length(samples)){
      p1 <- p1 + annotate("text", x=950, y=500000-30000*(m-1), label=samples[m], size=2)
    }
    
    ##Add line segments and set colors
    p2 <- p1+geom_segment(aes(x=800, y=500000-30000*(seq(1,length(samples),1)-1), xend=850, yend=500000-30000*(seq(1,length(samples),1)-1)),color = color[seq(1,length(samples),1)],size = 1)
    
    assign(paste0(tissue[i],stage[j]),p2)
  }
}

##


###save plot

setwd(paste0("F:\\ATAC\\ATAC-seqDatafinal_plot\\Insert size\\"))
pdf("InsertSize.pdf", width = 4.5, height = 5.5)
Br40
Br70
Br115
Br168
Li40
Li70
Li115
Li168
Mu40
Mu70
Mu115
Mu168
Pl40
Pl70
dev.off()
