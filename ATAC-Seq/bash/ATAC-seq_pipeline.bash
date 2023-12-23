##---------------------------------------##
#			ATAC-Seq data analysis       
##---------------------------------------##

##------------------------------------------------------ Based on original genome -------------------------------------------------
##1)The statistics of raw sequencing data
while read line
do
sample=$(echo $line | cut -f1)
sbatch --export=sample="$sample" /mnt/home/quanjian/quan/ATAC/script/rawdata_stat.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##2)The quality control of raw data
while read line
do
sample=$(echo $line | cut -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/preprocess/fastqc/"$sample".sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/preprocess/fastqc/"$sample".sbatch.err /mnt/home/quanjian/quan/ATAC/script/fastqc.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##3)Index generate of original reference
sbatch BWA_BuildIndex_Ori.sbatch


##4)Mapped the clean data to original Reference
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/alignment/"$sample"_OriRef.fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/alignment/"$sample"_OriRef.fq2gvcf.sbatch.err /mnt/home/quanjian/quan/ATAC/script/align_original.bwa.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##5)Redundant sequences are removed from the bam file, and mitochondrial DNA and low-quality Reads are also removed
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/picard/removeduplication/"$sample"_ReDup.fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/picard/removeduplication/"$sample"_ReDup.fq2gvcf.sbatch.err /mnt/home/quanjian/quan/ATAC/script/removeDuplication_OriRef.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##6)Alignment result statistic
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/stat/"$sample"_Stat.fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/stat/"$sample"_Stat.fq2gvcf.sbatch.err /mnt/home/quanjian/quan/ATAC/script/Statistic.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##7)Insert Size distribution
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/InsertSizes/"$sample"_Stat.fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/InsertSizes/"$sample"_Stat.fq2gvcf.sbatch.err /mnt/home/quanjian/quan/ATAC/script/CollectInsetSize.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##8)Bam files of three biological duplicate samples of the same developmental stage, same tissue, and same sex were merged
while read line
do
merge_id=$(echo $line | cut -d ' ' -f1)
rep1=$(echo $line | cut -d ' ' -f2)
rep2=$(echo $line | cut -d ' ' -f3)
rep3=$(echo $line | cut -d ' ' -f4)
sbatch --export=merge_id="$merge_id",rep1="$rep1",rep2="$rep2",rep3="$rep3" --output=/mnt/home/quanjian/quan/ATAC/log/bamMerge/"$sample"_mergeBam.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/bamMerge/"$sample"_mergeBam.sbatch.err /mnt/home/quanjian/quan/ATAC/script/mergeBamFiles_sex_tissue_stage.batch
done < /mnt/ufs18/home-081/quanjian/quan/ATAC/bash_data/merge_3rep.txt


##9)Convert the merged bam file to bed file and call the peaks according to bed file
while read line
do
merge_id=$(echo $line | cut -d ' ' -f1)
sbatch --export=merge_id="$merge_id" --output=/mnt/home/quanjian/quan/ATAC/log/peakcalling/"$sample"_Peakcalling.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/peakcalling/"$sample"_Peakcalling.sbatch.err /mnt/ufs18/home-081/quanjian/quan/ATAC/script/peakcalling_OriRef.merged.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/merge_id.txt


##10)Merge all peak files to generate union peak list
cd /mnt/scratch/quanjian/ATAC/peakCalling
/mnt/home/quanjian/software/bedtools2-2.30.0/bin/bedtools merge -i <(cat *_[FM]_*narrowPeak|cut -f1,2,3|sort -k1,1 -k2,2n) |awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t"$3-$2"\t""Peak_"NR}'|grep -v ^[AFXY].* > /mnt/scratch/quanjian/ATAC/peakCombine/AllPeakList.merged.txt


##11)Convert the union peak file with bed format to saf format
awk 'BEGIN{FS=OFS="\t"}{print $1,"PeakID","exon",$2+1,$3,".","+",".","gene_id ""\""$4"\"""\;"" ""transcript_id ""\""$6"\"""\;"}' /mnt/scratch/quanjian/ATAC/peakCombine/AllPeakList.merged.txt > /mnt/scratch/quanjian/ATAC/peakCombine/AllPeakList.merged.gtf
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $6, $1, $2+1, $3, "."}' /mnt/scratch/quanjian/ATAC/peakCombine/AllPeakList.merged.txt > /mnt/scratch/quanjian/ATAC/peakCombine/AllPeakList.merged.saf


##12)According to the final peaklist, the chromatin accessibility of each sample was counted
sbatch featureCount.batch


##13)Calculated the FRiP
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/FRiP/"$sample"_OriRef.sbatch.AllSample.log --error=/mnt/home/quanjian/quan/ATAC/log/FRiP/"$sample"_OriRef.sbatch.AllSample.err /mnt/ufs18/home-081/quanjian/quan/ATAC/script/FRiP_calculated.AllSamplemerged.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt


##14)Analyze the ATAC-Seq data
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/ATAC/log/deeptools_BamTobw/"$sample"_OriRef.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/deeptools_BamTobw/"$sample"_OriRef.sbatch.err /mnt/ufs18/home-081/quanjian/quan/ATAC/script/deeptools_analysis_ori.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/rawdata_list.txt






##-------------------------------------------------------- Based on modified genome ---------------------------------------------
##15)Build bwa index of each parental genome
for file in DB1 LS1 DB2 LS2 DB3 LS3 DB4 LS4 LB1 DS1 LB2 DS2 LB3 DS3 LB4 DS4
do
sbatch --export=mem=300G,env=/mnt/home/quanjian/quan/WGS/scripts/resource.env,dir=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome,sample="$file",fasta="$file"_lift.fa,out=/mnt/scratch/quanjian/ATAC/bwaIndex /mnt/home/quanjian/quan/ATAC/script/gMG_BWA_BuildIndex.sbatch
done


##16)Alignment based on modifited personlized genome
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
ref=$(echo $line | cut -d ' ' -f2)
sbatch --export=sample="$sample",ref="$ref" --output=/mnt/home/quanjian/quan/ATAC/log/alignment/"$sample"_"$ref".fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/alignment/"$sample"_"$ref".fq2gvcf.sbatch.err /mnt/home/quanjian/quan/ATAC/script/align.mod.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/sample_ref.txt2


##17)Check the integrity of bam files
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
ref=$(echo $line | cut -d ' ' -f2)
sbatch --export=sample="$sample",ref="$ref" --output=/mnt/home/quanjian/quan/ATAC/log/alignment/"$sample"_"$ref".fq2gvcf.check.sbatch.log --error=/mnt/home/quanjian/quan/ATAC/log/alignment/"$sample"_"$ref".fq2gvcf.check.sbatch.err /mnt/home/quanjian/quan/ATAC/script/AlignmentCheck.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/sample_ref.txt


##18)Statistical alignment result
cd /mnt/home/quanjian/quan/ATAC/bash_data/
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
ref=$(echo $line | cut -d ' ' -f2)
sbatch --export=sample="$sample",ref="$ref" /mnt/home/quanjian/quan/ATAC/script/Alignstastics.batch
done < sample_ref.txt


##19)Assigned reads count to each parental genome
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
p1bam=$(echo $line | cut -d ' ' -f2)
p2bam=$(echo $line | cut -d ' ' -f3)

#p1 maternal, p2 paternal
ref1=$(echo $line | cut -d ' ' -f4)
ref2=$(echo $line | cut -d ' ' -f5)
ParentPath=/mnt/scratch/quanjian/WGS/orginial_rawdata/fq2gvcf/personalizedGenome
mkdir -p /mnt/scratch/quanjian/ATAC/atacCount/$sample

sbatch --export=p1bed="$ParentPath"/"$ref1"/"$ref1"_lift.bed,p2bed="$ParentPath"/"$ref2"/"$ref2"_lift.bed,dir=/mnt/scratch/quanjian/ATAC/atacCount/$sample,env=/mnt/research/qgg/resource/ase/resource.env,p1bam="$p1bam",p2bam="$p2bam" --out=/mnt/scratch/quanjian/ATAC/atacCount/log/"$sample".out --error=/mnt/scratch/quanjian/ATAC/atacCount/log/"$sample".err  /mnt/ufs18/home-081/quanjian/quan/ATAC/script/atacCounts.sbatch
done < /mnt/home/quanjian/quan/ATAC/bash_data/atacCount.txt


##20)Annotated the peaks according to the coordinate of peaks that called out based on original genome
while read line
do
sample=$(echo $line | cut -d ' ' -f1)
sbatch --export=sample="$sample" /mnt/ufs18/home-081/quanjian/quan/ATAC/script/RefReadInGenome.batch
done < /mnt/home/quanjian/quan/ATAC/bash_data/atacCount.txt


##21)join the reads number of each peak to union peak list
while read line
do
sample=$(echo $line |cut -d " " -f 1)
dir="/mnt/gs21/scratch/quanjian/ATAC/atacCount"
peak="/mnt/gs21/scratch/quanjian/ATAC/peakCombine"

join -a1 -a2 -1 4 -2 4 -e '0' -o 0,1.1,1.2,1.3,2.5 <(sort -k4,4 "$peak"/AllPeakList.merged.txt) <(sort -k4,4 "$dir"/$sample/p1.ref.Mu168.count) > "$dir"/$sample/p1.PeakReadNumber.union.txt
join -a1 -a2 -1 4 -2 4 -e '0' -o 0,1.1,1.2,1.3,2.5 <(sort -k4,4 "$peak"/AllPeakList.merged.txt) <(sort -k4,4 "$dir"/$sample/p2.ref.Mu168.count) > "$dir"/$sample/p2.PeakReadNumber.union.txt
done < /mnt/home/quanjian/quan/ATAC/bash_data/atacCount.txt


##22)merged all joined file to get parental chromoatin accessibility matrix
paste <(cut -d ' ' -f1 ./DB115_Br1/p1.PeakReadNumber.union.txt) <(paste ./*/*PeakReadNumber.union.txt|awk 'FS=" "{for (i = 5; i <= NF; i+=5) printf("%s ", $i); printf("\n")}')|tr ' ' '\t' > /mnt/gs21/scratch/quanjian/ATAC/atacCount/AllGenomePeakReadNumber.txt
cat sample_ref AllGenomePeakReadNumber.txt > MG_AllGenomePeakReadNumber_alltitle.txt


##23)Make bedgraph using reads original coordinate for visualization
while read line
do
	sample=$(echo $line | cut -d ' ' -f1)

	#p1 maternal, p2 paternal
	ref1=$(echo $line | cut -d ' ' -f4)
	ref2=$(echo $line | cut -d ' ' -f5)
	
	sbatch --export=sample="$sample",ref1="$ref1",ref2="$ref2" --out=/mnt/scratch/quanjian/ATAC/atacCount/log/"$sample".bedgraph.out --error=/mnt/scratch/quanjian/ATAC/atacCount/log/"$sample".bedgraph.err  /mnt/home/quanjian/quan/ATAC/script/GetBedgraph.batch
done < <(grep '70' /mnt/home/quanjian/quan/ATAC/bash_data/atacCount.txt)
