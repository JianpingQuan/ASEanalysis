
#----------------------------------------------------------------
#              0.1 personal bismark database prepration           #
#----------------------------------------------------------------
##building index
for sample in DB2 DS2 LB2 LS2
do
	sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/BS/log/BuildIndex/"$sample".fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/BS/log/BuildIndex/"$sample".fq2gvcf.sbatch.err /mnt/home/quanjian/quan/BS/script/bismark_index_alignment/BismarkGenomePreparation.batch
done


#-----------------------------------------------------------------
#                1.Alignment,Deduplication                       #
#-----------------------------------------------------------------
#1)using Bismark to mapped reads based on modified genome
cd /mnt/scratch/quanjian/BS/Alignment
while read line
do
sample=$(echo $line|cut -f2|cut -d '/' -f7|cut -d '_' -f1,2)
r1=$(echo $line | cut -d " " -f2)
r2=$(echo $line | cut -d " " -f3)
ref=$(echo $line | cut -d " " -f1)

cd /mnt/scratch/quanjian/BS/Alignment/"$ref"
sbatch --export=sample="$sample",r1="$r1",r2="$r2",ref="$ref" --output=/mnt/home/quanjian/quan/BS/log/Alignment/fq2gvcf/"$sample"_Ref"$ref".fq2gvcf2.sbatch.log --error=/mnt/home/quanjian/quan/BS/log/Alignment/fq2gvcf/"$sample"_Ref"$ref".fq2gvcf2.sbatch.err /mnt/ufs18/home-081/quanjian/quan/BS/script/bismark_index_alignment/BismarkAlignment.batch

done < /mnt/home/quanjian/quan/BS/runfile/BSCutadaptTrimmedData_LS2.txt


#2)BismarkDeduplication.batch
while read line
do
sample=$(echo $line|cut -f2|cut -d '/' -f7|cut -d '_' -f1,2)
ref=$(echo $line | cut -d " " -f1)

sbatch --export=sample="$sample",ref="$ref" --output=/mnt/home/quanjian/quan/BS/log/Deduplication/"$sample"_Ref"$ref".sbatch.log --error=/mnt/home/quanjian/quan/BS/log/Deduplication/"$sample"_Ref"$ref".sbatch.err /mnt/home/quanjian/quan/BS/script/bismark_index_alignment/BismarkDeduplication.batch

done < /mnt/home/quanjian/quan/BS/runfile/BSCutadaptTrimmedData_LS2.txt


#3)check the integrity of bam file
cd /mnt/scratch/quanjian/BS/Deduplication
samtools quickcheck -v *deduplicated.bam > bad_bams.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'


#---------------------------------------------------------------------------------
#                            2. Assign reads to parental genome                  #
#---------------------------------------------------------------------------------
#0)methylation count step 
#
paste <(awk '$2 ~ /S2$/' BSDeduplicated.txt|awk '{print substr($0,1,8)}') <(awk '$2 ~ /S2$/' BSDeduplicated.txt|cut -f2) <(awk '$2 ~ /B2$/' BSDeduplicated.txt|cut -f2)|tr ' ' '\t' > Samples.methylCounts.txt


#1)Methylation site statistics of paternal and maternal parents (p1 is maternal, p2 is paternal), which includes three methylation types: CHH,CHG, and CpG
#
while read line
do
sample=$(echo $line|cut -d ' ' -f1)
ref1=$(echo $line|cut -d ' ' -f2)
ref2=$(echo $line|cut -d ' ' -f3)

p1bed=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome/"$ref1"/"$ref1"_lift.bed
p2bed=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome/"$ref2"/"$ref2"_lift.bed
dir=/mnt/scratch/quanjian/BS/methylCounts
env=/mnt/research/qgg/resource/ase/resource.env
p1bam=/mnt/scratch/quanjian/BS/Deduplication/"$sample"_Ref"$ref1".deduplicated.bam
p2bam=/mnt/scratch/quanjian/BS/Deduplication/"$sample"_Ref"$ref2".deduplicated.bam

sbatch --mem=256G --export=p1bed="$p1bed",p2bed="$p2bed",dir="$dir"/"$sample",env="$env",p1bam="$p1bam",p2bam="$p2bam" --out="$dir"/"$sample"/methyl.out --error="$dir"/"$sample"/methyl.err /mnt/research/qgg/software/scripts/methylCounts.sbatch
done <  /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt



#2)The row containing Z is extracted, and the bed file with one base before and after the C site is constructed
cd /mnt/scratch/quanjian/BS/methylCounts
while read line
do

sample=$(echo $line|cut -d ' ' -f1)
sbatch --export=file=$sample /mnt/home/quanjian/quan/BS/script/Filtering_CpG.batch
	
done < /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt



#3)Separate the C-base sites of different chromosomes, then merge by chromosome, then sort|uniq, and finally merge the different chromosomes together
#1. split Chr
cd /mnt/scratch/quanjian/BS/methylCounts
while read line
do

sample=$(echo $line|cut -d ' ' -f1)
sbatch --export=file=$sample /mnt/home/quanjian/quan/BS/script/Chr_CpG_Sep.batch
	
done < /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt


#2. Combine all chromosome and sort them
cd /mnt/scratch/quanjian/BS/methylCounts
seq 1 18|while read line
do
chr=$(echo $line)
sbatch --export=chr=$chr /mnt/home/quanjian/quan/BS/script/Chr_CpG_Union.batch
done


#3)Final CpG list combining different chromosomes
cat /mnt/scratch/quanjian/BS/methylCounts/meth.count.chr.*.merge |sort > /mnt/scratch/quanjian/BS/methylCounts/meth.count.chr.all.merge


#4)join meth.count.final for all samples with the union file
while read line
do
sample=$(echo $line|cut -d ' ' -f1)
sbatch --export=sample="$sample" /mnt/ufs18/home-081/quanjian/quan/BS/script/meth.count.join.batch
done < /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt


#5)Preliminary screening CpG was combined for all samples
cd /mnt/scratch/quanjian/BS/methylCounts
paste <(cut -d ' ' -f1 ./DB70M_B1/meth.count.final.join) <(paste *_B*/meth.count.final.join|awk '{ for (i = 2; i <= NF; i+=3) { printf $i " " $(i+1) " " } printf "\n" }')|tr ' ' '\t' >  /mnt/scratch/quanjian/BS/methylCounts/meth.count.Brain.txt
paste <(cut -d ' ' -f1 ./DB70M_B1/meth.count.final.join) <(paste *_L*/meth.count.final.join|awk '{ for (i = 2; i <= NF; i+=3) { printf $i " " $(i+1) " " } printf "\n" }')|tr ' ' '\t' >  /mnt/scratch/quanjian/BS/methylCounts/meth.count.Liver.txt
paste <(cut -d ' ' -f1 ./DB70M_B1/meth.count.final.join) <(paste *_M*/meth.count.final.join|awk '{ for (i = 2; i <= NF; i+=3) { printf $i " " $(i+1) " " } printf "\n" }')|tr ' ' '\t' >  /mnt/scratch/quanjian/BS/methylCounts/meth.count.Muscle.txt
paste <(cut -d ' ' -f1 ./DB70M_B1/meth.count.final.join) <(paste *_P*/meth.count.final.join|awk '{ for (i = 2; i <= NF; i+=3) { printf $i " " $(i+1) " " } printf "\n" }')|tr ' ' '\t' >  /mnt/scratch/quanjian/BS/methylCounts/meth.count.Placenta.txt


#Brain Liver Muscle Placenta
for tissue in Placenta Brain Liver Muscle
do

cd /mnt/scratch/quanjian/BS/methylCounts

#cat /mnt/scratch/quanjian/BS/methylCounts/meth.count."$tissue".txt|sed 's/Z://g'|sed 's/U://g'|sed 's/H://g'|sed 's/X://g' > OnlyCpG."$tissue".txt

paste <(cut -d ' ' -f1 ./DB70M_P1/meth.count.final.join) <(awk '{for (i = 2; i <= NF; i+=2) { printf $i " " } printf "\n" }'   /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".txt) > OnlyCpG."$tissue".p1.txt
paste <(cut -d ' ' -f1 ./DB70M_P1/meth.count.final.join) <(awk '{for (i = 3; i <= NF; i+=2) { printf $i " " } printf "\n" }' /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".txt) > OnlyCpG."$tissue".p2.txt
paste OnlyCpG."$tissue".p1.txt <(cut -f2- OnlyCpG."$tissue".p2.txt) > OnlyCpG."$tissue".p1_p2.txt
rm OnlyCpG."$tissue".p1.txt OnlyCpG."$tissue".p2.txt
done

cat OnlyCpG."$tissue".p1_p2.txt|tr ':' '\t'|awk '!($2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13==0)'|awk '!($14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25==0)'|awk '!($26+$27+$28+$29+$30+$31+$32+$33+$34+$35+$36+$37==0)'|awk '!($38+$39+$40+$41+$42+$43+$44+$45+$46+$47+$48+$49==0)' > OnlyCpG."$tissue".p1_p2.filter.txt

#Merge methylation and unmethylation sites for p1 and p2
awk '{print $1,$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24,$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25,$26+$28+$30+$32+$34+$36+$38+$40+$42+$44+$46+$48,$27+$29+$31+$33+$35+$37+$39+$41+$43+$45+$47+$49}' OnlyCpG."$tissue".p1_p2.filter.txt > OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth.txt
sed 's/_/\t/g' OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth.txt|sed 's/ /\t/g'|sort -k1,1 -k2,2n > OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth2.txt
rm OnlyCpG."$tissue".p1.txt OnlyCpG."$tissue".p2.txt OnlyCpG."$tissue".p1_p2.txt

done



for tissue in Placenta Brain Liver Muscle
do

#2.1)groupby according to the methLevel
cd ~/resources/genomeFeature
bedtools intersect -a gene.promoter.bed -b /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth2.txt -wa -wb|awk 'OFS="\t" {print $1,$2,$3,$4,$5,$9/($9+$10),$11/($11+$12)}'|bedtools groupby -i - -g 1,2,3,4,5 -c 6,7 -o mean > OnlyCpG."$tissue".genePromoter.methyLevel.mean.txt

#2.2)groupby according to the read count
bedtools intersect -a gene.promoter.bed -b /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth2.txt -wa -wb|cut -f1-5,9-12 |bedtools groupby -i - -g 1,2,3,4,5 -c 6,7,8,9 -o sum | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6/($6+$7),$8/($8+$9)}'> OnlyCpG."$tissue".genePromoter.readCount.sum.methLevel.txt

done



##6)Re-screening was performed based on the reads count number
#-----------------------------------------------------------------------------
## Coverage number greater than 5 in at least 4 samples

for tissue in Brain Liver Muscle Placenta
do
cd /mnt/scratch/quanjian/BS/methylCounts

awk 'OFS="\t"{print $1,$2+$3,$4+$5,$6+$7,$8+$9,$10+$11,$12+$13,$14+$15,$16+$17,$18+$19,$20+$21,$22+$23,$24+$25,$26+$27,$28+$29,$30+$31,$32+$33,$34+$35,$36+$37,$38+$39,$40+$41,$42+$43,$44+$45,$46+$47,$48+$49}' OnlyCpG."$tissue".p1_p2.filter.txt|awk '{
    count1=0; count2=0; count3=0; count4=0;
    for (i=2; i<=7; i++) {
        if ($i >= 5) {
            count1++;
            if (count1 >= 4) break;
        }
    }
    for (i=8; i<=13; i++) {
        if ($i >= 5) {
            count2++;
            if (count2 >= 4) break;
        }
    }
    for (i=14; i<=19; i++) {
        if ($i >= 5) {
            count3++;
            if (count3 >= 4) break;
        }
    }
    for (i=20; i<=25; i++) {
        if ($i >= 5) {
            count4++;
            if (count4 >= 4) break;
        }
    }
    if (count1 >= 4 && count2 >= 4 && count3 >= 4 && count4 >= 4) {
        print;
    }
}' > OnlyCpG."$tissue".p1_p2.filter.N_2.5.txt


##
join OnlyCpG."$tissue".p1_p2.filter.N_2.5.txt OnlyCpG."$tissue".p1_p2.filter.txt|cut -d ' ' -f1,26-|tr ' ' '\t' > OnlyCpG."$tissue".p1_p2.filter.N_2.5.final.txt


##Merge methylation and unmethylation sites for p1 and p2
awk '{print $1,$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24,$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25,$26+$28+$30+$32+$34+$36+$38+$40+$42+$44+$46+$48,$27+$29+$31+$33+$35+$37+$39+$41+$43+$45+$47+$49}' OnlyCpG."$tissue".p1_p2.filter.N_2.5.final.txt > OnlyCpG."$tissue".p1_p2.filter.N_2.5.Meth.UnMeth.txt

awk '{print $1,$2+$4+$6+$8+$10+$12+$38+$40+$42+$44+$46+$48,$3+$5+$7+$9+$11+$13+$39+$41+$43+$45+$47+$49,$14+$16+$18+$20+$22+$24+$26+$28+$30+$32+$34+$36,$15+$17+$19+$21+$23+$25+$27+$29+$31+$33+$35+$37}' OnlyCpG."$tissue".p1_p2.filter.N_2.5.final.txt > OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth.txt

done

for tissue in Brain Liver Muscle Placenta
do
sed 's/_/\t/g' OnlyCpG."$tissue".p1_p2.filter.N_2.5.Meth.UnMeth.txt|sed 's/ /\t/g'|sort -k1,1 -k2,2n > OnlyCpG."$tissue".p1_p2.filter.N_2.5.Meth.UnMeth2.txt

sed 's/_/\t/g' OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth.txt|sed 's/ /\t/g'|sort -k1,1 -k2,2n > OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth2.txt

done

for tissue in Brain Liver Muscle Placenta
do
#2.1)groupby according to the methLevel
bedtools intersect -a ~/resources/genomeFeature/gene.promoter3000.bed -b /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".p1_p2.filter.N_2.5.Meth.UnMeth2.txt -wa -wb|awk 'OFS="\t" {print $1,$2,$3,$4,$5,$9/($9+$10),$11/($11+$12)}'|bedtools groupby -i - -g 1,2,3,4,5 -c 6,7 -o mean > OnlyCpG."$tissue".p1_p2.N_2.5.genePromoter.methyLevel.mean3000.txt
bedtools intersect -a ~/resources/genomeFeature/gene.promoter3000.bed -b /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth2.txt -wa -wb|awk 'OFS="\t" {print $1,$2,$3,$4,$5,$9/($9+$10),$11/($11+$12)}'|bedtools groupby -i - -g 1,2,3,4,5 -c 6,7 -o mean > OnlyCpG."$tissue".D_L.N_2.5.genePromoter.methyLevel.mean3000.txt


#2.2)groupby according to the read count
bedtools intersect -a ~/resources/genomeFeature/gene.promoter3000.bed -b /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".p1_p2.filter.N_2.5.Meth.UnMeth2.txt -wa -wb|cut -f1-5,9-12 |bedtools groupby -i - -g 1,2,3,4,5 -c 6,7,8,9 -o sum | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6/($6+$7),$8/($8+$9)}'> OnlyCpG."$tissue".p1_p2.N_2.5.genePromoter.readCount.sum.methLevel3000.txt
bedtools intersect -a ~/resources/genomeFeature/gene.promoter3000.bed -b /mnt/scratch/quanjian/BS/methylCounts/OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth2.txt -wa -wb|cut -f1-5,9-12 |bedtools groupby -i - -g 1,2,3,4,5 -c 6,7,8,9 -o sum | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6/($6+$7),$8/($8+$9)}'> OnlyCpG."$tissue".D_L.N_2.5.genePromoter.readCount.sum.methLevel3000.txt

done

####how many AGE genes overlap with CpG
##33
join <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Brain.txt) <(cut -f5 OnlyCpG.Brain.D_L.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l

##32
join <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Liver.txt) <(cut -f5 OnlyCpG.Liver.D_L.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l

##42
join <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Muscle.txt) <(cut -f5 OnlyCpG.Muscle.D_L.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l

##38
join <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Placenta.txt) <(cut -f5 OnlyCpG.Placenta.D_L.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l


##how many POE genes overlap with CpG
##8
join <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Brain.txt) <(cut -f5 OnlyCpG.Brain.p1_p2.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l

##7
join <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Liver.txt) <(cut -f5 OnlyCpG.Liver.p1_p2.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l

##10
join <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Muscle.txt) <(cut -f5 OnlyCpG.Muscle.p1_p2.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l

##11
join <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Placenta.txt) <(cut -f5 OnlyCpG.Placenta.p1_p2.N_2.5.genePromoter.methyLevel.mean3000.txt|sort)| wc -l



###genome-wide methylation level between duroc and lulai
for tissue in Brain Liver Muscle Placenta
do
awk '{sum +=$6}END{print sum/NR}' OnlyCpG."$tissue".D_L.N_2.5.genePromoter.methyLevel.mean3000.txt
awk '{sum +=$7}END{print sum/NR}' OnlyCpG."$tissue".D_L.N_2.5.genePromoter.methyLevel.mean3000.txt
done

for tissue in Brain Liver Muscle Placenta
do
join -1 1 -2 5 <(sort /mnt/home/quanjian/resources/genomeFeature/AGE."$tissue".txt) <(sort -k5,5 OnlyCpG."$tissue".D_L.N_2.5.genePromoter.methyLevel.mean3000.txt) > "$tissue".AGE.overlap.genePromoter.methyLevel.mean3000.txt
done


##
join -o '0,1.2' /mnt/home/quanjian/resources/genomeFeature/AGE.Brain.Direction.txt Brain.AGE.overlap.genePromoter.methyLevel.mean3000.txt
join -o '0,1.2' /mnt/home/quanjian/resources/genomeFeature/AGE.Liver.Direction.txt Liver.AGE.overlap.genePromoter.methyLevel.mean3000.txt
join -o '0,1.2' /mnt/home/quanjian/resources/genomeFeature/AGE.Muscle.Direction.txt Muscle.AGE.overlap.genePromoter.methyLevel.mean3000.txt
join -o '0,1.2' /mnt/home/quanjian/resources/genomeFeature/AGE.Placenta.Direction.txt Placenta.AGE.overlap.genePromoter.methyLevel.mean3000.txt


##7)Calculate p value
sbatch --export=file="OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt",meta="BrainSample.txt" --out=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Brain.p1_p2.filter.out --error=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Brain.p1_p2.filter.err /mnt/ufs18/home-081/quanjian/quan/BS/script/Autosome.Csite.pvalue.batch

sbatch --export=file="OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt",meta="LiverSample.txt" --out=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Liver.p1_p2.filter.out --error=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Liver.p1_p2.filter.err /mnt/ufs18/home-081/quanjian/quan/BS/script/Autosome.Csite.pvalue.batch

sbatch --export=file="OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt",meta="MuscleSample.txt" --out=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Muscle.p1_p2.filter.out --error=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Muscle.p1_p2.filter.err /mnt/ufs18/home-081/quanjian/quan/BS/script/Autosome.Csite.pvalue.batch

sbatch --export=file="OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt",meta="PlacentaSample.txt" --out=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Placenta.p1_p2.filter.out --error=/mnt/home/quanjian/quan/BS/log/DifferentCsite.pvalue/OnlyCpG.Placenta.p1_p2.filter.err /mnt/ufs18/home-081/quanjian/quan/BS/script/Autosome.Csite.pvalue.batch


##
awk '$2<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l
awk '$2<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l
awk '$2<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l
awk '$2<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l

awk '$3<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l
awk '$3<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l
awk '$3<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l
awk '$3<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|uniq|wc -l


##----------------------------- poeMRs in the IGF2R gene range-------------------------------------------------------------------------------------------------
awk '$2<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$1==1&&$2>7367617&&$3<7477552'|awk '$4>1'
awk '$2<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$1==1&&$2>7367617&&$3<7477552'|awk '$4>1'
awk '$2<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$1==1&&$2>7367617&&$3<7477552'|awk '$4>1'
awk '$2<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$1==1&&$2>7367617&&$3<7477552'|awk '$4>1'


##Statistics about the number of peoMRs
awk '$2<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l
awk '$2<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l
awk '$2<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l
awk '$2<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l


#The number of ageMRs is collected
awk '$3<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l
awk '$3<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l
awk '$3<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l
awk '$3<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1'|wc -l


##The number of genes that overlap with poeMRs
paste gene.promoter3000.bed geneBody.bed|awk 'BEGIN {FS="\t"; OFS="\t"} {if ($4 == "+") {print $1, $2, $8, $4, $5} else if ($4 == "-") {print $1, $7, $3, $4, $5}}' > geneBody.ExtendPromoter.bed

bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|wc -l
bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|wc -l
bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|wc -l
bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|wc -l

##The number of POE genes that overlap with poeMRs
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Brain.txt)
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Liver.txt)
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Muscle.txt)
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$2<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/POE.Placenta.txt)


##The number of AGE genes that overlap with ageMRs
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$3<=0.05' OnlyCpG.Brain.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Brain.txt)
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$3<=0.05' OnlyCpG.Liver.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Liver.txt)
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$3<=0.05' OnlyCpG.Muscle.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Muscle.txt)
join <(bedtools intersect -a /mnt/home/quanjian/resources/genomeFeature/geneBody.ExtendPromoter.bed -b <(awk '$3<=0.05' OnlyCpG.Placenta.p1_p2.filter.N_2.5.final.txt.pvalue.txt|cut -f1|tr '_' '\t'|sort -k1,1 -k2,2n|bedtools merge -d 1500 -c 3 -o count|awk '$4>1') -wa -wb|cut -f5|sort|uniq|sort) <(sort /mnt/home/quanjian/resources/genomeFeature/AGE.Placenta.txt)




#-----------------------------------------------------------------------------------
#          3. Compare read number mapped to paternal and maternal genome           #
#-----------------------------------------------------------------------------------
#1)Paternal and maternal methylation levels of whole genome Promoter, Utr5, 1stExon, geneBody, intron and utr3 were measured
#cd /mnt/home/quanjian/resources
#gene.promoter.bed
#exon.1st.bed
#geneBody.bed
#geneid maternal_meth maternal_unmeth paternal_meth maternal_meth

##genome feature methylation level count
#
while read line
do	
	sample=$(echo $line|cut -d ' ' -f1)
	dir=/mnt/scratch/quanjian/BS/methylCounts
	
	sbatch --export=sample="$sample" --out="$dir"/"$sample"/methyl.genomeFeature.out --error="$dir"/"$sample"/methyl.genomeFeature.err /mnt/ufs18/home-081/quanjian/quan/BS/script/genomeFeature.meth.count.batch
	
done < ~/quan/BS/runfile/Samples.methylCounts.txt


##2)The methylation levels of POE gene and AGE gene in gene body, 1st exon and gene promoter were obtained
#1. POE gene list
#
cd /mnt/home/quanjian/resources
dos2unix POE.finalList.txt
cat ~/quan/BS/runfile/Samples.methylCounts.txt|cut -f1|sort|uniq|while read id
do
	grep -f /mnt/ufs18/home-081/quanjian/resources/POE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/gene1stExon.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/POE.gene1stExon.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/POE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/geneBody.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/POE.geneBody.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/POE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/genePromoter.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/POE.genePromoter.bed
done 


#2. AGE gene list
dos2unix AGE.finalList.txt

cat ~/quan/BS/runfile/Samples.methylCounts.txt|cut -f1|sort|uniq|while read id
do
	grep -f /mnt/ufs18/home-081/quanjian/resources/AGE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/gene1stExon.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/AGE.gene1stExon.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/AGE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/geneBody.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/AGE.geneBody.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/AGE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/genePromoter.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/AGE.genePromoter.bed
done 



###3.Statistical POE and AGE gene list
##POE
#
for file in `realpath */POE.gene1stExon.bed`
do 
	cat $file  >> allPOE.gene1stExon.list
done && cut -f1 allPOE.gene1stExon.list|sort|uniq > allUniqPOE.gene1stExon.list.txt && rm allPOE.gene1stExon.list
#
for file in `realpath */POE.geneBody.bed`
do 
	cat $file  >> allPOE.geneBody.list
done && cut -f1 allPOE.geneBody.list|sort|uniq > allUniqPOE.geneBody.list.txt && rm allPOE.geneBody.list
#
for file in `realpath */POE.genePromoter.bed`
do 
	cat $file  >> allPOE.genePromoter.list
done && cut -f1 allPOE.genePromoter.list|sort|uniq > allUniqPOE.genePromoter.list.txt && rm allPOE.genePromoter.list


##AGE
#
for file in `realpath */AGE.gene1stExon.bed`
do 
	cat $file  >> allAGE.gene1stExon.list
done && cut -f1 allAGE.gene1stExon.list|sort|uniq > allUniqAGE.gene1stExon.list.txt && rm allAGE.gene1stExon.list
#
for file in `realpath */AGE.geneBody.bed`
do 
	cat $file  >> allAGE.geneBody.list
done && cut -f1 allAGE.geneBody.list|sort|uniq > allUniqAGE.geneBody.list.txt && rm allAGE.geneBody.list
#
for file in `realpath */AGE.genePromoter.bed`
do 
	cat $file  >> allAGE.genePromoter.list
done && cut -f1 allAGE.genePromoter.list|sort|uniq > allUniqAGE.genePromoter.list.txt && rm allAGE.genePromoter.list


##Construct the union list methylation matrix of POE and AGE
realpath */AGE.genePromoter.bed|cut -d '/' -f8|while read sample
do
	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqPOE.geneBody.list.txt) <(sort -k1,1 ./$sample/POE.geneBody.bed) |tr ' ' '\t' > ./$sample/POE.geneBody.Union.bed

	join -a 1 -e 'O' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqPOE.genePromoter.list.txt) <(sort -k1,1 ./$sample/POE.genePromoter.bed) |tr ' ' '\t' > ./$sample/POE.genePromoter.Union.bed

	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqPOE.gene1stExon.list.txt) <(sort -k1,1 ./$sample/POE.gene1stExon.bed) |tr ' ' '\t' > ./$sample/POE.gene1stExon.Union.bed
	
	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqAGE.geneBody.list.txt) <(sort -k1,1 ./$sample/AGE.geneBody.bed) |tr ' ' '\t' > ./$sample/AGE.geneBody.Union.bed

	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqAGE.genePromoter.list.txt) <(sort -k1,1 ./$sample/AGE.genePromoter.bed) |tr ' ' '\t' > ./$sample/AGE.genePromoter.Union.bed

	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqAGE.gene1stExon.list.txt) <(sort -k1,1 ./$sample/AGE.gene1stExon.bed) |tr ' ' '\t' > ./$sample/AGE.gene1stExon.Union.bed
done


##Methylation level compute
awk '{print $1,$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24,$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25,$26+$28+$30+$32+$34+$36+$38+$40+$42+$44+$46+$48,$27+$29+$31+$33+$35+$37+$39+$41+$43+$45+$47+$49}' OnlyCpG.Brain.p1_p2.filter.txt | tr '_' '\t' |sort -k1,1 -k2,2n  > OnlyCpG.Brain.p1_p2.filter.Meth.UnMeth.txt

awk '{print $1,$2+$4+$6+$8+$10+$12+$38+$40+$42+$44+$46+$48,$3+$5+$7+$9+$11+$13+$39+$41+$43+$45+$47+$49,$14+$16+$18+$20+$22+$24+$26+$28+$30+$32+$34+$36,$15+$17+$19+$21+$23+$25+$27+$29+$31+$33+$35+$37}' OnlyCpG.Brain.p1_p2.filter.txt| tr '_' '\t' |sort -k1,1 -k2,2n > OnlyCpG.Brain.D_L.filter.Meth.UnMeth.txt


awk '{print $1,$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24,$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25,$26+$28+$30+$32+$34+$36+$38+$40+$42+$44+$46+$48,$27+$29+$31+$33+$35+$37+$39+$41+$43+$45+$47+$49}' OnlyCpG.Liver.p1_p2.filter.txt | tr '_' '\t' |sort -k1,1 -k2,2n  > OnlyCpG.Liver.p1_p2.filter.Meth.UnMeth.txt

awk '{print $1,$2+$4+$6+$8+$10+$12+$38+$40+$42+$44+$46+$48,$3+$5+$7+$9+$11+$13+$39+$41+$43+$45+$47+$49,$14+$16+$18+$20+$22+$24+$26+$28+$30+$32+$34+$36,$15+$17+$19+$21+$23+$25+$27+$29+$31+$33+$35+$37}' OnlyCpG.Liver.p1_p2.filter.txt| tr '_' '\t' |sort -k1,1 -k2,2n > OnlyCpG.Liver.D_L.filter.Meth.UnMeth.txt


awk '{print $1,$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24,$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25,$26+$28+$30+$32+$34+$36+$38+$40+$42+$44+$46+$48,$27+$29+$31+$33+$35+$37+$39+$41+$43+$45+$47+$49}' OnlyCpG.Muscle.p1_p2.filter.txt | tr '_' '\t' |sort -k1,1 -k2,2n  > OnlyCpG.Muscle.p1_p2.filter.Meth.UnMeth.txt

awk '{print $1,$2+$4+$6+$8+$10+$12+$38+$40+$42+$44+$46+$48,$3+$5+$7+$9+$11+$13+$39+$41+$43+$45+$47+$49,$14+$16+$18+$20+$22+$24+$26+$28+$30+$32+$34+$36,$15+$17+$19+$21+$23+$25+$27+$29+$31+$33+$35+$37}' OnlyCpG.Muscle.p1_p2.filter.txt| tr '_' '\t' |sort -k1,1 -k2,2n > OnlyCpG.Muscle.D_L.filter.Meth.UnMeth.txt


awk '{print $1,$2+$4+$6+$8+$10+$12+$14+$16+$18+$20+$22+$24,$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25,$26+$28+$30+$32+$34+$36+$38+$40+$42+$44+$46+$48,$27+$29+$31+$33+$35+$37+$39+$41+$43+$45+$47+$49}' OnlyCpG.Placenta.p1_p2.filter.txt | tr '_' '\t' |sort -k1,1 -k2,2n  > OnlyCpG.Placenta.p1_p2.filter.Meth.UnMeth.txt

awk '{print $1,$2+$4+$6+$8+$10+$12+$38+$40+$42+$44+$46+$48,$3+$5+$7+$9+$11+$13+$39+$41+$43+$45+$47+$49,$14+$16+$18+$20+$22+$24+$26+$28+$30+$32+$34+$36,$15+$17+$19+$21+$23+$25+$27+$29+$31+$33+$35+$37}' OnlyCpG.Placenta.p1_p2.filter.txt| tr '_' '\t' |sort -k1,1 -k2,2n > OnlyCpG.Placenta.D_L.filter.Meth.UnMeth.txt


##IGV visualization
for tissue in Brain Liver Muscle Placenta
do
awk 'OFS="\t" {print $1,$2,$3,$4/($4+$5)}' OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth.txt > OnlyCpG."$tissue".p1.methyLevel.bedgraph
awk 'OFS="\t" {print $1,$2,$3,$6/($6+$7)}' OnlyCpG."$tissue".p1_p2.filter.Meth.UnMeth.txt > OnlyCpG."$tissue".p2.methyLevel.bedgraph
done

for tissue in Brain Liver Muscle Placenta
do
awk 'OFS="\t" {print $1,$2,$3,$4/($4+$5)}' OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth2.txt > OnlyCpG."$tissue".Duroc.methyLevel.bedgraph
awk 'OFS="\t" {print $1,$2,$3,$6/($6+$7)}' OnlyCpG."$tissue".D_L.filter.N_2.5.Meth.UnMeth2.txt > OnlyCpG."$tissue".Lulai.methyLevel.bedgraph
done



