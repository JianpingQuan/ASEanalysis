#1. methylation level of promoter, utr'5, exon, intron, utr'3 and genebody

#1.1 get the position of each genome feature

#1) promoter (TSS forward 2000bp and TSS downstream 200bp)
# get gene bed information (should keep strand direction information)
cat Sus_scrofa.Sscrofa11.1.98.gtf|awk -v OFS="\t" '$3=="gene"{print $1,$4,$5,$4,$5,$7,$10}'|sed 's/"//g'|sed 's/;//g'|awk '$1 ~/^[1-9]/' > gene.bed

# get the TSS position of each gene according to the direction in strands
awk -v OFS="\t" '{if($6=="+"){print $1,$2,$2+1,$4,$5,$6,$7}else{print $1,$3-1,$3,$4,$5,$6,$7}}' gene.bed > gene.TSS.bed

# get the promoter region according to TSS position
bedtools slop -i gene.TSS.bed -g pig_chromosome_length.txt -l 2000 -r 199 -s|awk -v OFS="\t" '{print $1,$2,$3,$6,$7}' > gene.promoter.bed

#2) utr'5
paste <(awk '$1 !~ /^#/' Sus_scrofa.Sscrofa11.1.98.gff3|awk '$3=="five_prime_UTR"'|cut -f1,4,5) <(awk '$1 !~ /^#/' Sus_scrofa.Sscrofa11.1.98.gff3|awk '$3=="five_prime_UTR"'|cut -f7,9|sed 's/Parent=transcript://g')|sort -k1,1 -k2,2n > 5_utr.bed

#3) utr'3
paste <(awk '$1 !~ /^#/' Sus_scrofa.Sscrofa11.1.98.gff3|awk '$3=="three_prime_UTR"'|cut -f1,4,5) <(awk '$1 !~ /^#/' Sus_scrofa.Sscrofa11.1.98.gff3|awk '$3=="three_prime_UTR"'|cut -f7,9|sed 's/Parent=transcript://g')|sort -k1,1 -k2,2n > 3_utr.bed

#4) 1st_exon
# get 1st exon of each transcript
grep 'exon_number "1"' Sus_scrofa.Sscrofa11.1.98.gtf|awk -v OFS="\t" '$3=="exon" {print $1,$4,$5,$7,$10}'|sed 's/"//g'|sed 's/;//g'|awk '$1 ~/^[1-9]/' > exon.1st.bed

#keep the frist of exon_number "1" of the first transcript
awk '!gene[$5]++' exon.1st.bed > exon.1st.bed2
mv exon.1st.bed2 exon.1st.bed

#5) intron
# get chromSize bed file
cat Sus_scrofa.Sscrofa11.1.98.gff3 | awk -v OFS="\t" '$1 ~ /^#/ {print $2,"0",$4}'|awk '$1~ /^[1-9XY]/'|awk 'NR>1'|awk '$1<=18||$1=="X"||$1=="Y"'|sort -k1,1 -k2,2n > chromSizes.bed

# get whole exons position
awk '$1 !~ /^#/' Sus_scrofa.Sscrofa11.1.98.gff3|awk -v OFS="\t" '$3=="exon" {print $1,$4,$5}'|awk '$1 ~/^[1-9XY]/'|sort -k1,1 -k2,2n > exon.bed

# get gene body postion
cat Sus_scrofa.Sscrofa11.1.98.gtf | awk -v OFS="\t" '$3=="gene" {print $1,$4,$5,$7,$10}'|awk '$1 ~/^[1-9]/'|sed 's/"//g'|sed 's/;//g' > geneBody.bed

# get the intergenic position
bedtools subtract -a chromSizes.bed -b gene.bed > intergenic.bed

bedtools subtract -a gene.bed -b exon.bed > introns_intergenic.bed

bedtools subtract -a introns_intergenic.bed -b intergenic.bed|awk -v OFS="\t" '{print $1,$2,$3,$6,$7}' > intron.bed


#1.2 get the methylation level of each feature
Input1=`echo /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt`
Input2=`echo ~/resources`
OutPut=`echo /mnt/home/quanjian/featuresMethLevel`

##total coverage number and methylated number of each feature
#1) gene body 
bedtools intersect -a $Input2/geneBody.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nC.txt -wa -wb|cut -f1-3,5,9-55> $OutPut/AllSample.geneBody.methy.nCG.txt
bedtools intersect -a $Input2/geneBody.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nX.txt -wa -wb|cut -f1-3,5,9-55 > $OutPut/AllSample.geneBody.methy.nX.txt

#2) promoter
bedtools intersect -a $Input2/gene.promoter.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nC.txt -wa -wb|cut -f1-3,5,9-55> $OutPut/AllSample.genePromoter.methy.nCG.txt
bedtools intersect -a $Input2/gene.promoter.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nX.txt -wa -wb|cut -f1-3,5,9-55 > $OutPut/AllSample.genePromoter.methy.nX.txt

#3) utr'5
bedtools intersect -a $Input2/5_utr.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nC.txt -wa -wb|cut -f1-3,5,9-55> $OutPut/AllSample.geneUTR5.methy.nCG.txt
bedtools intersect -a $Input2/5_utr.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nX.txt -wa -wb|cut -f1-3,5,9-55 > $OutPut/AllSample.geneUTR5.methy.nX.txt

#4) utr'3
bedtools intersect -a $Input2/3_utr.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nC.txt -wa -wb|cut -f1-3,5,9-55> $OutPut/AllSample.geneUTR3.methy.nCG.txt
bedtools intersect -a $Input2/3_utr.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nX.txt -wa -wb|cut -f1-3,5,9-55 > $OutPut/AllSample.geneUTR3.methy.nX.txt

#5) 1st exon
bedtools intersect -a $Input2/exon.1st.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nC.txt -wa -wb|cut -f1-3,5,9-55> $OutPut/AllSample.gene1stExon.methy.nCG.txt
bedtools intersect -a $Input2/exon.1st.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nX.txt -wa -wb|cut -f1-3,5,9-55 > $OutPut/AllSample.gene1stExon.methy.nX.txt

#6) intron
bedtools intersect -a $Input2/intron.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nC.txt -wa -wb|cut -f1-3,5,9-55> $OutPut/AllSample.geneIntron.methy.nCG.txt
bedtools intersect -a $Input2/intron.bed -b $Input1/AllSamples.CpG.removeSNP.CovFilter.nX.txt -wa -wb|cut -f1-3,5,9-55 > $OutPut/AllSample.geneIntron.methy.nX.txt


## sum CpG site coverage number included in same gene through bedtools groupby module
#1) gene body
bedtools groupby -i AllSample.geneBody.methy.nCG.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneBody.methy.nCG.groupby.txt
bedtools groupby -i AllSample.geneBody.methy.nX.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneBody.methy.nX.groupby.txt

#2) promoter
bedtools groupby -i AllSample.genePromoter.methy.nCG.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.genePromoter.methy.nCG.groupby.txt
bedtools groupby -i AllSample.genePromoter.methy.nX.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.genePromoter.methy.nX.groupby.txt

#3) utr'5
bedtools groupby -i AllSample.geneUTR5.methy.nCG.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneUTR5.methy.nCG.groupby.txt
bedtools groupby -i AllSample.geneUTR5.methy.nX.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneUTR5.methy.nX.groupby.txt

#4) utr'3
bedtools groupby -i AllSample.geneUTR3.methy.nCG.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneUTR3.methy.nCG.groupby.txt
bedtools groupby -i AllSample.geneUTR3.methy.nX.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneUTR3.methy.nX.groupby.txt

#5) 1st exon
bedtools groupby -i AllSample.gene1stExon.methy.nCG.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.gene1stExon.methy.nCG.groupby.txt
bedtools groupby -i AllSample.gene1stExon.methy.nX.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.gene1stExon.methy.nX.groupby.txt

#6) intron
bedtools groupby -i AllSample.geneIntron.methy.nCG.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneIntron.methy.nCG.groupby.txt
bedtools groupby -i AllSample.geneIntron.methy.nX.txt -g 1,2,3,4 -c 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51 -o sum > AllSample.geneIntron.methy.nX.groupby.txt



## cd /mnt/home/quanjian/featuresMethLevel
#module purge
#module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2
#R
#
library(dplyr)
feature <- c("geneBody", "genePromoter", "geneUTR5", "geneUTR3", "gene1stExon", "geneIntron")

for(i in 1:6){
 nCG_feature <- paste0("AllSample.",feature[i],".methy.nCG.groupby.txt")
 nX_feature <- paste0("AllSample.",feature[i],".methy.nX.groupby.txt")

 nCG <- read.table(nCG_feature, header = F)
 nX  <- read.table(nX_feature, header = F)
 
 nCG <- nCG[,4:51]
 nX <- nX[,4:51]
 
 ## groupby by gene name, beacuse the utr5, utr3 and intron had multiple line with same gene
 nCG_groupby <- aggregate(. ~ V4, nCG, sum)
 nX_groupby <- aggregate(. ~ V4, nX, sum)

 MethyLevel <- nX_groupby[,2:48]/nCG_groupby[,2:48]
 rownames(MethyLevel) <- nCG_groupby[,1]
 
 output_feature <- paste0("AllSample.",feature[i], ".methyLevel.txt")
 
 write.table(MethyLevel, output_feature, quote = F, row.names = T,col.names = F,sep = '\t')
 }
