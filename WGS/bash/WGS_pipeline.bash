#1)Add read group information #
#=============================#
for sample in DB1 DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
  file1="$sample"_R1.fq.gz
  file2="$sample"_R2.fq.gz
  mkdir $sample
  sbatch --time=24:00:00 --export=dir=$sample,file1=$file1,file2=$file2 --output=$sample/fq.split.out --error=$sample/fq.split.err /mnt/research/qgg/resource/swim/splitFastQ.sbatch
done



##2)fastq file to gvcf###
#=======================#
for sample in DB1 DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
	bash /mnt/research/qgg/resource/swim/fq2gvcf.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/scratch/quanjian/WGS/orginial_rawdata/rawdata/ --tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp --out /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf --sample $sample > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf/$sample.out 2>&1 &
done



#3)creat hc.list##
#================#
# hc1 = chr 1-4
# hc2 = chr 5-9
# hc3 = chr 10-14
# hc4 = chr 15-18, X (ploidy = 2)
# hc5 = chr X (plody = 1)
find ./ -name "hc1.*.g.vcf.gz" | sed 's/^\.\///' | awk '{print "/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf/"$1}' > hc1.auto.file.list

find ./ -name "hc2.*.g.vcf.gz" | sed 's/^\.\///' | awk '{print "/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf/"$1}' > hc2.auto.file.list

find ./ -name "hc3.*.g.vcf.gz" | sed 's/^\.\///' | awk '{print "/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf/"$1}' > hc3.auto.file.list

find ./ -name "hc4.*.g.vcf.gz" | sed 's/^\.\///' | awk '{print "/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf/"$1}' > hc4.auto.file.list

find ./ -name "hc*.*.g.vcf.gz" | awk '$1 ~ /hc4/ || $1 ~ /hc5/' | awk '($1 ~ /B/ && $1 ~ /hc5/) || ($1 ~ /S/ && $1 ~ /hc4/)' | sed 's/^\.\///' | sort | awk '{print "/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/gvcf/"$1}' > hc5.x.file.list


#4.1)combine gvcf file of autosome#
#=================================#
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
	if (($chr>=1)) && (($chr<=4));then
		list=$(echo hc1.auto.file.list|awk '{print $0}')
		bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF --tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp --list $list --int $chr --pre $chr > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF/log/$chr.out 2>&1 &
	elif (($chr>=5))&&(($chr<=9));then
		list=$(echo hc2.auto.file.list|awk '{print $0}')
		bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF --tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp --list $list --int $chr --pre $chr > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF/log/$chr.out 2>&1 &
	elif (($chr>=10))&&(($chr<=14));then
		list=$(echo hc3.auto.file.list|awk '{print $0}')
		bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF --tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp --list $list --int $chr --pre $chr > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF/log/$chr.out 2>&1 &
	elif (($chr>=15))&&(($chr<=18));then
		list=$(echo hc4.auto.file.list|awk '{print $0}')
		bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash --resource /mnt/research/qgg/resource/swim --dir /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF --tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp --list $list --int $chr --pre $chr > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF/log/$chr.out 2>&1 &
	fi
done

#4.2) combine gvcf file of X chromosome#
#======================================#			
chr=X
bash /mnt/research/qgg/resource/swim/combineGVCFswait.bash \
--resource /mnt/research/qgg/resource/swim \
--dir /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF \
--tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp \
--list hc5.x.file.list \
--int $chr \
--pre $chr > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF/log/$chr.out 2>&1 &



##5)genotypeGVCF#
#===============#
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/combineGVCF
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X
do
	bash /mnt/research/qgg/resource/swim/genotypeGVCFswait.bash \
	--resource /mnt/research/qgg/resource/swim \
	--dir /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF \
	--tmp /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp \
	--input "$chr".g.vcf.gz \
	--pre "$chr" > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/"$chr".out 2>&1 &
done



#6)Merge all Chromosome #
#=======================#
for chr in $(seq 1 18) X; do realpath $chr.vcf.gz; done > allChr.vcf.list

sbatch --export=mem=64G,\
env=/mnt/research/qgg/resource/swim/resource.env,\
input=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/allChr.vcf.list,\
output=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/allChr.vcf.gz,\
tmp=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp,\
out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/,prefix=allChr \
--output=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.mergeVCFs.out \
--error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.mergeVCFs.err \
/mnt/research/qgg/resource/swim/sbatch/mergeVCFs.sbatch



#7)split VCF into Indel and SNP#
#==============================#
sbatch --export=mem=64G,\
env=/mnt/research/qgg/resource/swim/resource.env,\
vcf=allChr.vcf.gz,output1=allChr.snp.vcf.gz,\
output2=allChr.indel.vcf.gz,\
tmp=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp,\
out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/ \
--output=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.splitVCF.out \
--error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.splitVCF.err \
/mnt/research/qgg/resource/swim/sbatch/splitVCF.sbatch


##8)excessHet filtering##
#=======================#
for variant in indel snp
do
	sbatch --export=mem=32G,\
	env=/mnt/research/qgg/resource/swim/resource.env,\
	input=allChr.$variant.vcf.gz,\
	tmp=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp,\
	output=allChr.$variant.excessHet.vcf.gz,\
	dir=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF,\
	prefix=allChr.$variant \
	--output=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.$variant.excessHet.out \
	--error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.$variant.excessHet.err \
	/mnt/research/qgg/resource/swim/sbatch/excessHet.sbatch
done


#9)Indel filtering #
#==================#
sbatch --export=mem=64G,\
env=/mnt/research/qgg/resource/swim/resource.env,\
vcf=allChr.indel.excessHet.vcf.gz,\
output=allChr.indel.excessHet.hardFilter.vcf.gz,\
tmp=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp,\
out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF \
--error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.indel.excessHet.hardFilter.err \
/mnt/research/qgg/resource/swim/sbatch/hardFilterIndel.sbatch


#10)SNP filtering###
#==================#
sbatch --export=mem=64G,\
env=/mnt/research/qgg/resource/swim/resource.env,\
vcf=allChr.snp.excessHet.vcf.gz,\
output=allChr.snp.excessHet.hardFilter.vcf.gz,\
tmp=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/tmp,\
out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF \
--error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/log/allChr.snp.excessHet.hardFilter.err \
/mnt/home/quanjian/quan/WGS/scripts/hardFilterSNP.sbatch


#11)combine the hardfiltered SNP and Indel,select the PASS SNP and Indel vairants#
#=========================================#
bcftools concat -a allChr.snp.excessHet.hardFilter.vcf.gz allChr.indel.excessHet.hardFilter.vcf.gz > allChr.variant.hardFilter.vcf.gz
bcftools view -f 'PASS,.' allChr.variant.hardFilter.vcf.gz -o allChr.variant.hardFilter.finalPASS.vcf.gz



#12) extract the hard filter genotype information fron VCF file##
#===============================================================#
bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n' /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/allChr.variant.hardFilter.finalPASS.vcf.gz > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome/allChr.variant.hardFilter.FinalGenotype.txt


#13)prepared bed file based on vcf file (autosome and X chromosome)##
#===================================================================#
#This step produces BED file with animals whose variants are output for subsequent steps. The output contains information on reference sequence, the variant sequence.
#upper case = homo with certainty controled by some of the flags; lower case = het or homo with high uncertainty.
#--vcf vcf file name
#--geno non-missing rate in the groups
#--freq frequency of the major allele to call the consensus is a homo
#--group1 animals in first group, this must include all animals belonging to the group in the VCF file. The list can be generated by `grep \#CHROM 1_re_8.vcf | cut -f 10- | tr '\t' '\n' | grep -v semen | tr '\n' ','`, in this case, a `yorkshire` is added at the end so a consensus sequence is called
#--group2 only the `semen` sample in this example

perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/prepSeqBED.pl --vcf /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF/allChr.variant.hardFilter.finalPASS.vcf.gz --group1 DB1,DB2,DB3,DB4,DS1,DS2,DS3,DS4 --group2 LB1,LB2,LB3,LB4,LS1,LS2,LS3,LS4  --freq 0.99 --geno 0.874 > /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome/allChr.variant.FinalPass.bed


#14)extract the variants of match paired parental cross with Paternal_Maternal order
#remove sites that were same with reference in both parental genome
#only keep sites that different with reference at one least parent
#====================================================================================#
path="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome"
file="allChr.variant.FinalPass.bed"

awk '{print $1,$2,$3,$4,$5,$6,$18}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB1_LS1_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$7,$19}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB2_LS2_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$8,$20}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB3_LS3_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$9,$21}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/DB4_LS4_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$14,$10}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB1_DS1_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$15,$11}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB2_DS2_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$16,$12}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB3_DS3_AutoX.bed

awk '{print $1,$2,$3,$4,$5,$17,$13}' "$path"/"$file"|awk '!(toupper($6) == toupper($4) && toupper($7) == toupper($4))' > "$path"/LB4_DS4_AutoX.bed


##15)split the paired cross variants to each parent##
#===================================================#
#-----------------------------For Maternal genome--------------------------
#"chr,start,end,ref,Mat,id"
path="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome"
for file in `echo [LD]B*_AutoX.bed`
do
	tail -n+2 "$file"| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$7,$5}' > "$path"/`echo "$file" | cut -d "_" -f 2 `_prelift.bed
done

#----------------------------For Paternal genome---------------------------
#"chr,start,end,ref,Pat,id"
path="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome"
for file in `echo [LD]B*_AutoX.bed`
do
	tail -n+2 "$file"| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' > "$path"/`echo "$file" | cut -d "_" -f 1 `_prelift.bed
done


#16)lift gtf and produce new reference based on above bed file
#=============================================================#
#According to parental variant sites to generate their own gtf and reference.
dir="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome"
for sample in DB1 LS1 DB2 LS2 DB3 LS3 DB4 LS4 LB1 DS1 LB2 DS2 LB3 DS3 LB4 DS4
do
	perl /mnt/home/quanjian/quan/WGS/scripts/aseUtils/liftGTF.pl \
	--gtf /mnt/ufs18/home-081/quanjian/resources/Sus_scrofa.Sscrofa11.1.98.gtf \
	--bed "$dir"/"$sample"_prelift.bed \
	--fasta /mnt/research/qgg/resource/sscrofa11.1/genome/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa \
	--gtfout "$dir"/"$sample"_lift.gtf \
	--varout "$dir"/"$sample"_lift.bed \
	--fastaout "$dir"/"$sample"_lift.fa
done



#----------------------------------------------------------------------------------------------------------------
#17)convert finalPass vcf file to plink file by vcftools#
#======================================================#
sbatch --export=mem=30G,env=/mnt/home/quanjian/quan/WGS/scripts/resource.env,gzvcf=allChr.variant.hardFilter.finalPASS.vcf.gz,dir=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/genotypeGVCF,out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/plink --error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/plink/log/allChr.variant.err /mnt/ufs18/home-081/quanjian/quan/WGS/scripts/gVcfToPlink.batch 


#18)creat relatedness matrix#
#===========================#
/mnt/research/qgg/software/plink-v1.90b6.18/plink --gzvcf allChr.variant.hardFilter.finalPASS.vcf.gz --maf 0.05 --make-bed --out snp --chr-set 18 no-xy
/mnt/research/qgg/software/gcta_1.93.2beta/gcta64 --make-grm --out snp.gcta --bfile snp --autosome-num 18
/mnt/research/qgg/software/gcta_1.93.2beta/gcta64 --grm snp.gcta --pca 10 --out snp.gcta


#19)Fst value between two breeds
#calculate Fst value of each slid windows
/mnt/home/quanjian/software/vcftools/bin/vcftools --gzvcf allChr.variant.hardFilter.finalPASS.vcf.gz --weir-fst-pop Duroc_population.txt --weir-fst-pop Lulai_population.txt  --out Duroc_Lulai_bin --fst-window-size 100000 --fst-window-step 100000



#20)get summary data of F0 sample
#===============================#

#raw reads number
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/map

ls *flagstat|while read id
do
awk 'NR==1' $id|cut -d ' ' -f1
done >> RawReadNum.txt

#properly paired reads number
ls *flagstat|while read id
do
awk 'NR==9' $id|cut -d ' ' -f1
done >> ProperlyPairedReadNum.txt

#duplicated reads number
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/markDuplicates
ls *alg.mark.metrics.txt|while read id
do
awk 'NR==8' $id|cut -f9
done >> PercentageDup.txt

#deduplicate reads average coverage (%) and average depth
for sample in DB1 DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
sbatch --export=sample=$sample /mnt/ufs18/home-081/quanjian/quan/WGS/scripts/CoverageDepth.batch
done

ls *.alg.mark.depth|while read id
do
cat $id|tr ',' '\t'|awk 'NR>1{print $0}'|awk '{sum+=$3}END{print sum/NR}'
done

#extract SNP and INDEL information of each sample from vcf file
for file in DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
bcftools view -s "$file" allChr.variant.hardFilter.finalPASS.vcf.gz|/mnt/home/quanjian/software/vcftools/bin/vcftools --gzvcf - --remove-indels --recode --recode-INFO-all --out "$file"_SNP_only

bcftools view -s "$file" allChr.variant.hardFilter.finalPASS.vcf.gz|/mnt/home/quanjian/software/vcftools/bin/vcftools --gzvcf - --keep-only-indels --recode --recode-INFO-all --out "$file"_INDEL_only
done


#SNP number
for file in DB1 DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n' "$file"_SNP_only.recode.vcf|cut -d ' ' -f6|grep -E '1|2|3|4'|wc -l
done

#INDEL number
for file in LS3 LS4
do
bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n'  "$file"_INDEL_only.recode.vcf|cut -d ' ' -f3,4,6|grep -E '1|2|3|4|5|6|7|8|9'|wc -l 
done


#informative variants number in exon
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome
for file in DB1 DB2 DB3 DB4 LB1 LB2 LB3 LB4
do
	awk '$6 ~ /[ATCG]/ && $7 ~ /[ATCG]/ && $6 != $7' */"$file"_*_AutoX.bed|tr ' ' '\t'|bedtools intersect -a - -b <(awk '$3 == "exon"' ~/resources/Sus_scrofa.Sscrofa11.1.98.gtf) -wa | uniq | cut -f 5 | sort | uniq | wc -l
done


#gene number can be probed by informatibe variants
for file in DB1 DB2 DB3 DB4 LB1 LB2 LB3 LB4
do
	awk '$6 ~ /[ATCG]/ && $7 ~ /[ATCG]/ && $6 != $7' */"$file"_*_AutoX.bed|tr ' ' '\t'| bedtools intersect -a - -b <(awk '$3 == "exon"' ~/resources/Sus_scrofa.Sscrofa11.1.98.gtf) -wa -wb |cut -f 5,16 | awk -F "\"" '{print $1"\t"$2}' | awk '{print $3"\t"$1}' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o count_distinct|wc -l 
done


##Total gene in Autosomes and X chromosome of reference genome
awk '$3=="gene"' ~/resources/Sus_scrofa.Sscrofa11.1.98.gtf|awk '$1<=18||$1=="X"'|wc -l

##Total gene that expressed in autosome and X chromosome
dat <- read.csv("F:\\RNA\\DEG_Tissues_Stages\\GeneExpressionAllSamples.csv", header = T, row.names = 1)
library(dplyr)
dat2 <- filter(dat, rowSums(dat)>=6)
dim(dat2)

##informative number(in exon)
for file in DB1 DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do
	bedtools intersect -a $file/"$file"_informative.bed -b <(sort -k1,1 -k2,2n $file/"$file".lift.exon.bed) -wa|wc -l
done


##gene
for file in DB1 DB2 DB3 DB4 DS2 DS3 DS4 LB1 LB2 LB3 LB4 LS1 LS2 LS3 LS4
do 
	bedtools intersect -a $file/"$file"_informative.bed -b <(sort -k1,1 -k2,2n $file/"$file".lift.exon.bed) -wa -wb|cut -f11|sort|uniq|wc -l
done
