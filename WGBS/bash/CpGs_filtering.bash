#1. found the SNPs that need to be removed from subsequent.
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF
# extract individual data from merged VCF file
bcftools view -s DB2,LS2,LB2,DS2 allChr.snp.excessHet.hardFilter.vcf.gz > F0_Stage70.snp.excessHet.hardFilter.vcf

# ectract genotype information from VCF file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' F0_Stage70.snp.excessHet.hardFilter.vcf > F0_Stage70.snp.genotype.txt

# because we only consider the CpG site in autosome, so we also only analysis the SNPs in autosome.
# we according to the genotype information of four F0 individual to filter SNPs that need to be remove in subsequent analysis. 
# need removed SNPs standard: DB2/LS2 have one or two not 0/0 or 0|0, or DS2/LB2 have one or two not 0/0 or 0|0.
# 
awk '$1 ~/[^1-9]/' F0_Stage70.snp.genotype.txt|awk '((($5!="0/0")&&($5!="0|0")&&($5!="./.")&&($5!=".|."))||(($6!="0/0")&&($6!="0|0")&&($6!="./.")&&($6!=".|."))||(($7!="0/0")&&($7!="0|0")&&($7!="./.")&&($7!=".|."))||(($8!="0/0")&&($8!="0|0")&&($8!="./.")&&($8!=".|.")))' > NeedRemoveSNP.txt


#2. remove the CpG sites not in autosome and convert to dss format
#covert CpG site in CpG_report.txt.gz file to dss format, and only contain CpGs in autosome. 
#CpG_report.txt.gz contain all CpGs in genome, so all CpG_report.txt.gz files contain the same lines.

cd /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt
zcat "$sample"_OriRef.deduplicated.sortbyName.CpG_report.txt.gz|awk 'FS=OFS="\t" {print $1,$2,$4+$5,$4}'|awk '$1 ~/^[0-9]/'|sort -k1,1n -k2,2n > /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt/"$sample".dss.txt

# the format of dss.txt was: chromosome,start,N,X
# N was coverage number of each CpGs
# X was methylated number of each CpGs


#3. remove the CpG site with low coverage (coverage >= 10)
#module purge
#module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2
#R

library(stringr)
library(dplyr)

## For each tissue
for(i in 1:4){
 
 ##get the total coverage number of each CpG sites
 get_nC <- function(x) {
 dat_tmp <- read.table(x, header = F, stringsAsFactors = F)
 nC <- dat_tmp[,3]
 return(nC)
 }

 ##get the position of each CpG
 position <- read.table("DB70M_B1.dss.txt", header = F, stringsAsFactors = F)
 CpG_pos <- paste(position$V1, position$V2, sep="_") %>% as.data.frame()
 colnames(CpG_pos) <- c("chr_start")
 
 Tis <- c("B", "L", "M", "P")
 Tissue <- c("Brain", "Liver", "Muscle", "Placenta")

 pattern <- paste0("*",Tis[i],"[123].dss.txt")
 files <- list.files(pattern = pattern)
 data <- lapply(files, function(x) get_nC(x))
 tis <- as.data.frame(data)
 colnames(tis) <- str_sub(files, 1, 8)

 ##filtering by converage(read number >10)
 tis_if <- as.data.frame(tis[,1:ncol(tis)]>=10)
 tis$num <- rowSums(tis_if)
 tis_combine <- cbind(CpG_pos, tis)
 
 ##every group should have at least 8 sample >10 at each CpG site
 tis_filter <- filter(tis_combine, num>=8)
 dataSave <- tis_filter[,1]

 ##leaved CpG site
 outputName <- paste0(Tissue[i], ".filtered.CpG.txt") 
 write.table(dataSave, outputName, quote = F, sep="\t", row.names=F, col.names=T)
 
 ## remove all menmory
 rm(list=ls())
}


#4.  convert NeedRemoveSNP.txt to bed file and remove NeedRemoveSNP from unique CpGs got of above
# A. convert NeedRemoveSNP.txt to bed file
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF
awk '{print $1"\t"$2-1"\t"$2}' NeedRemoveSNP.txt|sort -k1,1 -k2,2n > NeedRemoveSNP.bed

# B. Get the unique CpG site of all Tissues combind CpG sites
cd /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt
cat <(awk 'NR>1' Brain.filtered.CpG.txt) <(awk 'NR>1' Liver.filtered.CpG.txt) <(awk 'NR>1' Muscle.filtered.CpG.txt) <(awk 'NR>1' Placenta.filtered.CpG.txt)|cut -f1|sort|uniq|tr '_' '\t'|sort -k1,1 -k2,2n|awk '{print $1"\t"$2-1"\t"$2}' > All.CovNumFilterd.MergedUniqueCpG.txt

# C. Remove the SNPs in NeedRemoveSNP.bed from All.CovNumFilterd.UniqueCpG.txt
bedtools intersect -a /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt/All.CovNumFilterd.MergedUniqueCpG.txt -b /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/NeedRemoveSNP.bed -v|awk '{print $1"_"$2}' > /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt/Autosome.MergedUniqueCpG.SNPfilter.bed


#5. substract the CpGs sites in Autosome.MergedUniqueCpG.SNPfilter.bed from each sample.dss.txt file
cd /mnt/home/quanjian/quan/BS/runfile
while read line
do
 sample=$(echo $line |cut -d " " -f1)
 sbatch --export=sample="$sample"  ~/quan/BS/script/Sample.CpG.removeSNP.CovFilter.batch
done < ~/quan/BS/runfile/Samples.list.txt


#6. combine the extracted CpGs of each sample, used for subsequent analysis
library(dplyr)
library(stringr)
get_Meth <- function(x) {
 dat_tmp <- read.table(x, header = F, stringsAsFactors = F)
 Meth <- dat_tmp[,4]/dat_tmp[,3]*100
 return(Meth)
}

get_nC <- function(x) {
 dat_tmp <- read.table(x, header = F, stringsAsFactors = F)
 nC <- dat_tmp[,3]
 return(nC)
}

get_nX <- function(x) {
 dat_tmp <- read.table(x, header = F, stringsAsFactors = F)
 nX <- dat_tmp[,4]
 return(nX)
}

position <- read.table("DB70M_B1.CpG.removeSNP.CovFilter.txt", header = F, stringsAsFactors = F)

CpG_pos <- paste(position$V1, position$V2, sep="_") %>% as.data.frame()

colnames(CpG_pos) <- c("start")

files <- list.files(pattern = "*CpG.removeSNP.CovFilter.txt")

dat_level <- lapply(files, function(x) get_Meth(x))
dat_nC <- lapply(files, function(x) get_nC(x))
dat_nX <- lapply(files, function(x) get_nX(x))

tis_level <- as.data.frame(dat_level)
tis_nC <- as.data.frame(dat_nC)
tis_nX <- as.data.frame(dat_nX)

colnames(tis_level) <- str_sub(files, 1, 8)
colnames(tis_nC) <- str_sub(files, 1, 8)
colnames(tis_nX) <- str_sub(files, 1, 8)

tis_combine_level <- cbind(CpG_pos, tis_level)
tis_combine_nC <- cbind(CpG_pos, tis_nC)
tis_combine_nX <- cbind(CpG_pos, tis_nX)

write.table(tis_combine_level, "AllSamples.CpG.removeSNP.CovFilter.methLevel.txt", quote = F, sep="\t", row.names=F, col.names=T)
write.table(tis_combine_nC, "AllSamples.CpG.removeSNP.CovFilter.nC.txt", quote = F, sep="\t", row.names=F, col.names=T)
write.table(tis_combine_nX, "AllSamples.CpG.removeSNP.CovFilter.nX.txt", quote = F, sep="\t", row.names=F, col.names=T)


#7. convert the file to bed file
paste <(awk 'NR>1' AllSamples.CpG.removeSNP.CovFilter.methLevel.txt|cut -f1|tr '_' '\t'|awk '{print $1"\t"$2-1"\t"$2}') <(awk 'NR>1' AllSamples.CpG.removeSNP.CovFilter.methLevel.txt|cut -f2-48) > AllSamples.CpG.removeSNP.CovFilter.methLevel.txt2

paste <(awk 'NR>1' AllSamples.CpG.removeSNP.CovFilter.nC.txt|cut -f1|tr '_' '\t'|awk '{print $1"\t"$2-1"\t"$2}') <(awk 'NR>1' AllSamples.CpG.removeSNP.CovFilter.nC.txt|cut -f2-48) > AllSamples.CpG.removeSNP.CovFilter.nC.txt2

paste <(awk 'NR>1' AllSamples.CpG.removeSNP.CovFilter.nX.txt|cut -f1|tr '_' '\t'|awk '{print $1"\t"$2-1"\t"$2}') <(awk 'NR>1' AllSamples.CpG.removeSNP.CovFilter.nX.txt|cut -f2-48) > AllSamples.CpG.removeSNP.CovFilter.nX.txt2

mv AllSamples.CpG.removeSNP.CovFilter.methLevel.txt2 AllSamples.CpG.removeSNP.CovFilter.methLevel.txt
mv AllSamples.CpG.removeSNP.CovFilter.nC.txt2 AllSamples.CpG.removeSNP.CovFilter.nC.txt
mv AllSamples.CpG.removeSNP.CovFilter.nX.txt2 AllSamples.CpG.removeSNP.CovFilter.nX.txt