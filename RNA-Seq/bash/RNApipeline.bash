#-------------------------------------------------------------------------------
#                                A_Orginal Genome                              #
#-------------------------------------------------------------------------------

#0. combine the fastq files from two pairs to one pair,so that we can use the same alignment pipeline
#generate the reference index
sbatch Hisat2ExtractSpliceSite.batch
sbatch Hisat2ExtractExon.batch
sbatch Hisat2Build.batch


#1. Use the Hisat2 and samtools to map reads and statistic reads number and remove duplication
while read line
do
sample=$(echo $line | cut -d " " -f 1)
fq1=$(echo $line | cut -d " " -f 2)
fq2=$(echo $line | cut -d " " -f 3)

out=/mnt/scratch/quanjian/RNA_seq/oriAlignment
env=/mnt/ufs18/home-081/quanjian/resources/resourceHome.env

sbatch --export=sample="$sample",fq1="$fq1",fq2="$fq2",out="$out",env="$env" /mnt/home/quanjian/quan/RNA/script/oriHisat2_AlignAndReDup_OnePair.batch
done < /mnt/home/quanjian/quan/RNA/ori.alignment.OnePairSampleList.txt


#2. Get gene expression of each sample
while read line
do
sample=$(echo $line | cut -d " " -f 1)

sbatch --export=sample="$sample" /mnt/ufs18/home-081/quanjian/quan/RNA/script/featureCount.batch
done < /mnt/home/quanjian/quan/RNA/ori.alignment.OnePairSampleList.txt


#3. variants calling using RNA-Seq data
#3.1 baseRecalibrator
cd /mnt/home/quanjian/quan/RNA/
while read line
do
sample=$(echo $line | cut -d " " -f 1)

sbatch --export=sample="$sample" --output=log/gatk/baseRecalibrator/"$sample".BQSR.sbatch.log --error=log/gatk/baseRecalibrator/"$sample".BQSR.sbatch.err /mnt/home/quanjian/quan/RNA/script/baseRecalibrator.batch
done < Sample.bam.txt

#3.2 applyBQSR
cd /mnt/home/quanjian/quan/RNA/
while read line
do
sample=$(echo $line | cut -d " " -f 1)

sbatch --export=sample="$sample" --output=log/gatk/applyBQSR/"$sample".apply.sbatch.log --error=log/gatk/applyBQSR/"$sample".apply.sbatch.err /mnt/home/quanjian/quan/RNA/script/applyBQSR.batch
done < Sample.bam.txt

#3.3 haplotype caller
cd /mnt/home/quanjian/quan/RNA/
myarr=("1:1-12500000" "1:12500001-25000000" "1:25000001-37500000" "1:37500001-50000000" "1:50000001-62500000" "1:62500001-75000000" "1:75000001-87500000" "1:87500001-100000000" "1:100000001-112500000" "1:112500001-125000000" "1:125000001-137500000" "1:137500001-150000000" "1:150000001-162500000" "1:162500001-175000000" "1:175000001-187500000" "1:187500001-200000000" "1:200000001-212500000" "1:212500001-225000000" "1:225000001-237500000" "1:237500001-250000000" "1:250000001-262500000" "1:262500001-274330532" "2:1-12500000" "2:12500001-25000000" "2:25000001-37500000" "2:37500001-50000000" "2:50000001-62500000" "2:62500001-75000000" "2:75000001-87500000" "2:87500001-100000000" "2:100000001-112500000" "2:112500001-125000000" "2:125000001-137500000" "2:137500001-150000000" "2:150000001-151935994" "3:1-12500000" "3:12500001-25000000" "3:25000001-37500000" "3:37500001-50000000" "3:50000001-62500000" "3:62500001-75000000" "3:75000001-87500000" "3:87500001-100000000" "3:100000001-112500000" "3:112500001-125000000" "3:125000001-132848913" "4:1-12500000" "4:12500001-25000000" "4:25000001-37500000" "4:37500001-50000000" "4:50000001-62500000" "4:62500001-75000000" "4:75000001-87500000" "4:87500001-100000000" "4:100000001-112500000" "4:112500001-125000000" "4:125000001-130910915" "5:1-12500000" "5:12500001-25000000" "5:25000001-37500000" "5:37500001-50000000" "5:50000001-62500000" "5:62500001-75000000" "5:75000001-87500000" "5:87500001-100000000" "5:100000001-104526007" "6:1-12500000" "6:12500001-25000000" "6:25000001-37500000" "6:37500001-50000000" "6:50000001-62500000" "6:62500001-75000000" "6:75000001-87500000" "6:87500001-100000000" "6:100000001-112500000" "6:112500001-125000000" "6:125000001-137500000" "6:137500001-150000000" "6:150000001-162500000" "6:162500001-170843587" "7:1-12500000" "7:12500001-25000000" "7:25000001-37500000" "7:37500001-50000000" "7:50000001-62500000" "7:62500001-75000000" "7:75000001-87500000" "7:87500001-100000000" "7:100000001-112500000" "7:112500001-121844099" "8:1-12500000" "8:12500001-25000000" "8:25000001-37500000" "8:37500001-50000000" "8:50000001-62500000" "8:62500001-75000000" "8:75000001-87500000" "8:87500001-100000000" "8:100000001-112500000" "8:112500001-125000000" "8:125000001-137500000" "8:137500001-138966237" "9:1-12500000" "9:12500001-25000000" "9:25000001-37500000" "9:37500001-50000000" "9:50000001-62500000" "9:62500001-75000000" "9:75000001-87500000" "9:87500001-100000000" "9:100000001-112500000" "9:112500001-125000000" "9:125000001-137500000" "9:137500001-139512083" "10:1-12500000" "10:12500001-25000000" "10:25000001-37500000" "10:37500001-50000000" "10:50000001-62500000" "10:62500001-69359453" "11:1-12500000" "11:12500001-25000000" "11:25000001-37500000" "11:37500001-50000000" "11:50000001-62500000" "11:62500001-75000000" "11:75000001-79169978" "12:1-12500000" "12:12500001-25000000" "12:25000001-37500000" "12:37500001-50000000" "12:50000001-61602749" "13:1-12500000" "13:12500001-25000000" "13:25000001-37500000" "13:37500001-50000000" "13:50000001-62500000" "13:62500001-75000000" "13:75000001-87500000" "13:87500001-100000000" "13:100000001-112500000" "13:112500001-125000000" "13:125000001-137500000" "13:137500001-150000000" "13:150000001-162500000" "13:162500001-175000000" "13:175000001-187500000" "13:187500001-200000000" "13:200000001-208334590" "14:1-12500000" "14:12500001-25000000" "14:25000001-37500000" "14:37500001-50000000" "14:50000001-62500000" "14:62500001-75000000" "14:75000001-87500000" "14:87500001-100000000" "14:100000001-112500000" "14:112500001-125000000" "14:125000001-137500000" "14:137500001-141755446" "15:1-12500000" "15:12500001-25000000" "15:25000001-37500000" "15:37500001-50000000" "15:50000001-62500000" "15:62500001-75000000" "15:75000001-87500000" "15:87500001-100000000" "15:100000001-112500000" "15:112500001-125000000" "15:125000001-137500000" "15:137500001-140412725" "16:1-12500000" "16:12500001-25000000" "16:25000001-37500000" "16:37500001-50000000" "16:50000001-62500000" "16:62500001-75000000" "16:75000001-79944280" "17:1-12500000" "17:12500001-25000000" "17:25000001-37500000" "17:37500001-50000000" "17:50000001-62500000" "17:62500001-63494081" "18:1-12500000" "18:12500001-25000000" "18:25000001-37500000" "18:37500001-50000000" "18:50000001-55982971")
for loc in ${myarr[@]}
myarr=("1:1-12500000" "1:12500001-25000000" "1:25000001-37500000" "1:37500001-50000000" "1:50000001-62500000" "1:62500001-75000000" "1:75000001-87500000" "1:87500001-100000000" "1:100000001-112500000" "1:112500001-125000000" "1:125000001-137500000" "1:137500001-150000000" "1:150000001-162500000" "1:162500001-175000000" "1:175000001-187500000" "1:187500001-200000000" "1:200000001-212500000" "1:212500001-225000000" "1:225000001-237500000" "1:237500001-250000000" "1:250000001-262500000" "1:262500001-274330532" "2:1-12500000" "2:12500001-25000000" "2:25000001-37500000" "2:37500001-50000000" "2:50000001-62500000" "2:62500001-75000000" "2:75000001-87500000" "2:87500001-100000000" "2:100000001-112500000" "2:112500001-125000000" "2:125000001-137500000" "2:137500001-150000000" "2:150000001-151935994" "3:1-12500000" "3:12500001-25000000" "3:25000001-37500000" "3:37500001-50000000" "3:50000001-62500000" "3:62500001-75000000" "3:75000001-87500000" "3:87500001-100000000" "3:100000001-112500000" "3:112500001-125000000" "3:125000001-132848913" "4:1-12500000" "4:12500001-25000000" "4:25000001-37500000" "4:37500001-50000000" "4:50000001-62500000" "4:62500001-75000000" "4:75000001-87500000" "4:87500001-100000000" "4:100000001-112500000" "4:112500001-125000000" "4:125000001-130910915" "5:1-12500000" "5:12500001-25000000" "5:25000001-37500000" "5:37500001-50000000" "5:50000001-62500000" "5:62500001-75000000" "5:75000001-87500000" "5:87500001-100000000" "5:100000001-104526007" "6:1-12500000" "6:12500001-25000000" "6:25000001-37500000" "6:37500001-50000000" "6:50000001-62500000" "6:62500001-75000000" "6:75000001-87500000" "6:87500001-100000000" "6:100000001-112500000" "6:112500001-125000000" "6:125000001-137500000" "6:137500001-150000000" "6:150000001-162500000" "6:162500001-170843587" "7:1-12500000" "7:12500001-25000000" "7:25000001-37500000" "7:37500001-50000000" "7:50000001-62500000" "7:62500001-75000000" "7:75000001-87500000" "7:87500001-100000000" "7:100000001-112500000" "7:112500001-121844099" "8:1-12500000" "8:12500001-25000000" "8:25000001-37500000" "8:37500001-50000000" "8:50000001-62500000" "8:62500001-75000000" "8:75000001-87500000" "8:87500001-100000000" "8:100000001-112500000" "8:112500001-125000000" "8:125000001-137500000" "8:137500001-138966237" "9:1-12500000" "9:12500001-25000000" "9:25000001-37500000" "9:37500001-50000000" "9:50000001-62500000" "9:62500001-75000000" "9:75000001-87500000" "9:87500001-100000000" "9:100000001-112500000" "9:112500001-125000000" "9:125000001-137500000" "9:137500001-139512083" "10:1-12500000" "10:12500001-25000000" "10:25000001-37500000" "10:37500001-50000000" "10:50000001-62500000" "10:62500001-69359453" "11:1-12500000" "11:12500001-25000000" "11:25000001-37500000" "11:37500001-50000000" "11:50000001-62500000" "11:62500001-75000000" "11:75000001-79169978" "12:1-12500000" "12:12500001-25000000" "12:25000001-37500000" "12:37500001-50000000" "12:50000001-61602749" "13:1-12500000" "13:12500001-25000000" "13:25000001-37500000" "13:37500001-50000000" "13:50000001-62500000" "13:62500001-75000000" "13:75000001-87500000" "13:87500001-100000000" "13:100000001-112500000" "13:112500001-125000000" "13:125000001-137500000" "13:137500001-150000000" "13:150000001-162500000" "13:162500001-175000000" "13:175000001-187500000" "13:187500001-200000000" "13:200000001-208334590" "14:1-12500000" "14:12500001-25000000" "14:25000001-37500000" "14:37500001-50000000" "14:50000001-62500000" "14:62500001-75000000" "14:75000001-87500000" "14:87500001-100000000" "14:100000001-112500000" "14:112500001-125000000" "14:125000001-137500000" "14:137500001-141755446" "15:1-12500000" "15:12500001-25000000" "15:25000001-37500000" "15:37500001-50000000" "15:50000001-62500000" "15:62500001-75000000" "15:75000001-87500000" "15:87500001-100000000" "15:100000001-112500000" "15:112500001-125000000" "15:125000001-137500000" "15:137500001-140412725" "16:1-12500000" "16:12500001-25000000" "16:25000001-37500000" "16:37500001-50000000" "16:50000001-62500000" "16:62500001-75000000" "16:75000001-79944280" "17:1-12500000" "17:12500001-25000000" "17:25000001-37500000" "17:37500001-50000000" "17:50000001-62500000" "17:62500001-63494081" "18:1-12500000" "18:12500001-25000000" "18:25000001-37500000" "18:37500001-50000000" "18:50000001-55982971")
for loc in ${myarr[@]}
do
chr=`echo $loc | cut -d ":" -f1`
num=`echo $loc | cut -d "-" -f2`
int=`expr $num / 12500001 + 1`

sbatch --export=LOC="$loc",Chr="$chr",Order="$int" --output=log/gatk/haplotypeCaller/Chr"$chr"."$int".haplotypecaller.sbatch.log --error=log/gatk/haplotypeCaller/Chr"$chr"."$int".haplotypecaller.sbatch.err /mnt/home/quanjian/quan/RNA/script/haplotypeCaller.batch
done

#3.4 merge the haplotypecaller results of each region from the same chromasome
cd /mnt/scratch/quanjian/RNA_seq/GATK/HaplotypeCaller/
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
	cat Chr"$chr".1.allSample.vcf > Chr"$chr".allSample.vcf
	num=`ls Chr"$chr".[0-9]*.vcf|wc -l`
	for ((i=2; i<="$num"; i=i+1))
	do
		tail -n+640 Chr"$chr"."$i".allSample.vcf >> Chr"$chr".allSample.vcf
	done
done

#3.5 merge the autosome
cd /mnt/scratch/quanjian/RNA_seq/GATK/HaplotypeCaller/
cat Chr1.allSample.vcf > ChrAll.allSample.vcf
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
	tail -n+640 Chr"$chr".allSample.vcf >> AllAuto.allSample.vcf
done

#3.6 generate the ibs matrix of all individuals according to the variants infomation in autosome
/mnt/research/qgg/software/plink-v1.90b5.3/plink --vcf /mnt/scratch/quanjian/RNA_seq/GATK/HaplotypeCaller/Autosome.allSample.vcf --geno 0.2 --maf 0.05 --recode --double-id --out /mnt/scratch/quanjian/RNA_seq/plink/Autosome.allSample.plink
/mnt/research/qgg/software/plink-v1.90b5.3/plink --file /mnt/scratch/quanjian/RNA_seq/plink/Autosome.allSample.plink --make-bed --out Autosome.allSample
cd mnt/scratch/quanjian/RNA_seq/plink/
/mnt/research/qgg/software/plink-v1.90b5.3/plink --file /mnt/scratch/quanjian/RNA_seq/plink/Autosome.allSample.plink --distance square
mv plink.dist Autosome.allSample.ibs.matrix.txt
mv plink.dist Autosome.allSample.ibs.matrixID.txt


#4. generate the bigwig file used for IGV visualization
while read line
do
sample=$(echo $line | cut -d " " -f 1)
env=/mnt/ufs18/home-081/quanjian/resources/resourceHome.env

sbatch --export=sample="$sample",env="$env" /mnt/home/quanjian/quan/RNA/script/bamToBW.batch

done < /mnt/home/quanjian/quan/RNA/ori.alignment.OnePairSampleList.txt




#-----------------------------------------------------------------------------------------
#                                        B_Modified Genome                               #
#-----------------------------------------------------------------------------------------
#---------------------------------------0. generate lifted bed file -----------------------

#0.1 For Maternal genome
#"chr,start,end,ref,Mat,id"
path="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome"
for file in `echo [LD]B*_AutoX.bed`
do
	tail -n+2 "$file"| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$7,$5}' > "$path"/`echo "$file" | cut -d "_" -f 2 `_prelift.bed
done


#0.2 For Paternal genome
#"chr,start,end,ref,Pat,id"
path="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome"
for file in `echo [LD]B*_AutoX.bed`
do
	tail -n+2 "$file"| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' > "$path"/`echo "$file" | cut -d "_" -f 1 `_prelift.bed
done


#0.3 lift gtf and produce new reference based on above bed file
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



#------------------- 1. build index ----------------
#1.1 extract the splice and exon in lifted gtf file
#==================================================#
for file in DB1 LS1 DB2 LS2 DB3 LS3 DB4 LS4 LB1 DS1 LB2 DS2 LB3 DS3 LB4 DS4
do
	sbatch --export=mem=100G,env=/mnt/home/quanjian/quan/WGS/scripts/resource.env,dir=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome,gtf="$file"_lift.gtf,sample="$file",fasta="$file"_lift.fa,out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedIndex --output=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedIndex/"$file"/log/"$file".spliceExon.sbatch.log --error=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedIndex/"$file"/log/"$file".spliceExon.err /mnt/ufs18/home-081/quanjian/quan/RNA/script/gMG_Hisat2_extract_spliceExon.sbatch
done


#1.2 Build index of each parental individual
#=========================================#
for file in DB1 LS1 DB2 LS2 DB3 LS3 DB4 LS4 LB1 DS1 LB2 DS2 LB3 DS3 LB4 DS4
do
	sbatch --export=mem=300G,env=/mnt/home/quanjian/quan/WGS/scripts/resource.env,dir=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome,sample="$file",fasta="$file"_lift.fa,out=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedIndex /mnt/ufs18/home-081/quanjian/quan/RNA/script/gMG_Hisat2_BuildIndex.sbatch
done


#----------------------------- 2.alignment --------------------------------
#
while read line
do
	sample=$(echo $line | cut -d " " -f 1)
	fq1=$(echo $line | cut -d " " -f 2)
	fq2=$(echo $line | cut -d " " -f 3)
	ref=$(echo $line | cut -d " " -f 6)
	Refdir="/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedIndex"
	out="/mnt/scratch/quanjian/RNA_seq/gAlignment"

	sbatch --export=sample="$sample",env=/mnt/home/quanjian/quan/WGS/scripts/resource.env,fq1="$fq1",fq2="$fq2",\
	ref="$ref",out="$out",Refdir="$Refdir" \
	--output="$out"/"$sample"/log/"$sample".ori.sbatch.log \
	--error="$out"/"$sample"/log/"$sample".ori.sbatch.err /mnt/ufs18/home-081/quanjian/quan/RNA/script/gHisat2_Align_OnePair.batch
done < /mnt/home/quanjian/quan/RNA/g.alignment.OnePairSampleList.txt



#--------------------- 3.assign reads to maternal and paternal genome --------------------------
#
while read line
do
	sample=$(echo $line|cut -d ' ' -f1)
	#paternal genome
	ref1=$(echo $line|cut -d ' ' -f2)
	#maternal genome
	ref2=$(echo $line|cut -d ' ' -f3)
	#the path save parental modified genome 
	ParentPath=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome

	p1bed="$ParentPath"/"$ref1"/"$ref1"_lift.bed
	p2bed="$ParentPath"/"$ref2"/"$ref2"_lift.bed
	p1gtf="$ParentPath"/"$ref1"/"$ref1"_lift.gtf
	p2gtf="$ParentPath"/"$ref2"/"$ref2"_lift.gtf
	p1ref="$ParentPath"/"$ref1"/"$ref1"_lift.fa
	p2ref="$ParentPath"/"$ref2"/"$ref2"_lift.fa
	p1bam=/mnt/scratch/quanjian/RNA_seq/gAlignment/$sample/"$sample"_"$ref1".align.bam
	p2bam=/mnt/scratch/quanjian/RNA_seq/gAlignment/$sample/"$sample"_"$ref2".align.bam
	dir=/mnt/scratch/quanjian/RNA_seq/gReadsCount/$sample
	env=/mnt/ufs18/home-081/quanjian/resources/resourceHome.env
	
	sbatch --export=p1bed="$p1bed",p2bed="$p2bed",dir="$dir",env="$env",p1bam="$p1bam",p2bam="$p2bam",p1ref="$p1ref",p2ref="$p2ref",p1gtf="$p1gtf",p2gtf="$p2gtf" --out=/mnt/scratch/quanjian/RNA_seq/gReadsCount/log/"$sample".rnaCounts.out --error=/mnt/scratch/quanjian/RNA_seq/gReadsCount/log/"$sample".rnaCounts.err /mnt/research/qgg/software/scripts/rnaCounts.sbatch

done < /mnt/ufs18/home-081/quanjian/quan/RNA/g.sample.pat.mat.txt



#---------------- 4.Generate gene expression matrix (both maternal and paternal) -----------------
#
#4.1 unique gene list
#
cd /mnt/scratch/quanjian/RNA_seq/gReadsCount
for file in `realpath */*gene.count`
do 
	cat $file  >> allGene.list
done && cut -f1 allGene.list |sort|uniq > allUniqGene.list.txt && rm allGene.list


#4.2 join each sample with unique gene list
#
cd /mnt/scratch/quanjian/RNA_seq/gReadsCount
for  path in `realpath *|grep -v 'log'|grep -v 'allUniqGene.list.txt'`
do
	join -t $'\t' -a 1 -a 2 -e '0' -o '0,2.2' allUniqGene.list.txt <(sort -k1,1 $path/p1.gene.count) > $path/p1.allUniqueGene.count
	join -t $'\t' -a 1 -a 2 -e '0' -o '0,2.2' allUniqGene.list.txt <(sort -k1,1 $path/p2.gene.count) > $path/p2.allUniqueGene.count
done


#4.3 combine all joined samples
#
cd /mnt/scratch/quanjian/RNA_seq/gReadsCount

#all gene expression matrix
paste */*allUniqueGene.count|awk '{for(i=3;i<=NF;i+=2) $i="" }1'| tr ' ' '\t' > AllGeneCountFinal.matrix.txt

#autosome gene expression matrix
grep -v -f /mnt/home/quanjian/resources/NoAutosome.gene.list.bed allUniqGene.list.txt > Autosome.UniqGene.list.txt
grep -f Autosome.UniqGene.list.txt AllGeneCountFinal.matrix.txt > AutosomeGeneCountFinal.matrix.txt



##------------------------------ 5. generate bedgraph file ------------------------------
#5.1 build offset.bed file
source /mnt/research/qgg/resource/ase/resource.env
cd /mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome
for file in DB1 DS1 LB1 LS1 DB2 DS2 LB2 LS2 DB3 DS3 LB3 LS3 DB4 DS4 LB4 LS4
do
	#Chr NewStart NewEnd REF ALT OFFSET
	awk 'length($4) != length($5)' ./"$file"/"$file"_lift.bed | perl $OFFSET > ./"$file"/"$file".offset.bed
done


#5.2 intersect the mapping bam file under parental genome with offset.bed
while read line
do
	sample=$(echo $line|cut -d ' ' -f1)
	ref=$(echo $line|cut -d ' ' -f2)
	#the path save parental modified genome 
	ParentPath=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome
	bam=/mnt/scratch/quanjian/RNA_seq/gAlignment/$sample/"$sample"_"$ref".align.bam
	dir=/mnt/scratch/quanjian/RNA_seq/gReadsCount/$sample
	env=/mnt/ufs18/home-081/quanjian/resources/resourceHome.env
	sbatch --export=dir="$dir",env="$env",bam="$bam",ParentPath="$ParentPath",ref="$ref",sample="$sample" /mnt/ufs18/home-081/quanjian/quan/RNA/script/bamToCoverage.batch
done < /mnt/ufs18/home-081/quanjian/quan/RNA/Samples_ref.txt


#5.3 parental origin reads were used for analyzing
#p1 refer reads from paternal genome，p2 refer reads from meternal genome
while read line
do
	sample=$(echo $line|cut -d ' ' -f1)
	sbatch --export=sample="$sample" /mnt/ufs18/home-081/quanjian/quan/RNA/script/GetFinalReads.batch
done < /mnt/ufs18/home-081/quanjian/quan/RNA/g.sample.pat.mat.txt


#5.4 extract the p1 and p2 reads from bed file, and generate bedgraph file
while read line
do
	sample=$(echo $line|cut -d ' ' -f1)
	#paternal genome
	ref1=$(echo $line|cut -d ' ' -f2)
	#maternal genome
	ref2=$(echo $line|cut -d ' ' -f3)
	sbatch --export=sample="$sample",ref1="$ref1",ref2="$ref2" /mnt/ufs18/home-081/quanjian/quan/RNA/script/GetBedgraph.batch
done < /mnt/ufs18/home-081/quanjian/quan/RNA/g.sample.pat.mat.txt

