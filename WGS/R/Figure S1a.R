####the diversity of variats in each chromosome 
setwd("C:/Users/qjp/Desktop/Manuscript/WGS/data")
chr_len <- read.table("pig_karyotype.txt",header=T)
p <- c(paste0('chr',1:19))
for(i in 1:nrow(chr_len)){
  len=chr_len$End[i]
  s=seq(1,len,100000)
  e=seq(100000,len,100000)
  chr=rep(chr_len$Chr[i],length(s))
  if(length(s)>length(e)){
    e[length(s)]=len }
  assign(p[i],data.frame(Chr=chr,start=s,end=e))
}

chr_all <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19)

write.csv(chr_all,"chr_all.csv",quote=F,row.names = F)


##make file used for building ideogram
pig_karyotype_98<-read.table("Sus_scrofa.Sscrofa11.1.98.gff3", stringsAsFactors = F, header = F, comment.char = "#", sep = "\t", quote = "")
pig_karyotype_98 <- pig_karyotype_98[pig_karyotype_98$V3=="chromosome",]
pig_karyotype_98 <- pig_karyotype_98[,c(1,4,5)]
colnames(pig_karyotype_98) <- c("Chr","Start","End")
pig_karyotype_98 <- pig_karyotype_98[c(order(as.numeric(pig_karyotype_98$Chr[1:18]))),]
ideogram(karyotype = pig_karyotype_98)
convertSVG("chromosome.svg", device = "png")
write.table(pig_karyotype_98, "pig_karyotype2.txt",sep = "\t",row.names = F)

#density file format ：Chr Start End Value
#the value snp density file is the SNP number in each windows (End - Start)  


####chromosome summary plot
library(RIdeogram)
pear_karyotype <- read.table("pig_karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
duroc_SNP_density <- read.csv("Duroc_density.100K.csv", header = T, stringsAsFactors = F)
duroc_SNP_density2 <- data.frame(Chr=duroc_SNP_density$Chr, Start=duroc_SNP_density$Start, End=duroc_SNP_density$End, Value=log2(duroc_SNP_density$Value))
lulai_SNP_density <- read.csv("Lulai_density.100K.csv", header = T, stringsAsFactors = F)
lulai_SNP_density2 <- data.frame(Chr=lulai_SNP_density$Chr, Start=lulai_SNP_density$Start, End=lulai_SNP_density$End, Value=log2(lulai_SNP_density$Value))


#SNP_markers <- read.table("SNP_markers.txt", sep = "\t", header = T, stringsAsFactors = F)

#ideogram(karyotype = pear_karyotype, overlaid = SNP_density2,label = SNP_markers, label_type = "marker", colorset1 = c("darkblue", "white", "red"),Lx = 80, Ly = 25)

#ideogram(karyotype = pear_karyotype, overlaid = SNP_density2,colorset1 = c("#F2AF66","white", "darkblue"),Lx = 80, Ly = 25)

ideogram(karyotype = pear_karyotype, overlaid = duroc_SNP_density2,label = lulai_SNP_density2,label_type = "heatmap", colorset1 = c("white", "darkblue"), colorset2 = c("white","#458A45"))

convertSVG("chromosome.svg", device = "png", dpi = 600)


