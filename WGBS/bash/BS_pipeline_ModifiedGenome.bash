#----------------------------------------------------------------
#              0.1 personal bismark database prepration           #
#----------------------------------------------------------------
##构建bismark软件比对所需要的索引文件(bismar指定bowtie2的参数变了，已经由--path_to_bowtie改为了--path_to_aligner)
for sample in DB2 DS2 LB2 LS2
do
	sbatch --export=sample="$sample" --output=/mnt/home/quanjian/quan/BS/log/BuildIndex/"$sample".fq2gvcf.sbatch.log --error=/mnt/home/quanjian/quan/BS/log/BuildIndex/"$sample".fq2gvcf.sbatch.err /mnt/home/quanjian/quan/BS/script/bismark_index_alignment/BismarkGenomePreparation.batch
done


#-----------------------------------------------------------------
#                1.Alignment,Deduplication,Extractor             #
#-----------------------------------------------------------------
#1)using Bismark to mapped reads based on modified genome
##要在scrath上提交任务，会在提交任务的文件夹产生大量临时文件。
cd /mnt/scratch/quanjian/BS/Alignment
while read line
do
	sample=$(echo $line|cut -f2|cut -d '/' -f7|cut -d '_' -f1,2)
	r1=$(echo $line | cut -d " " -f2)
	r2=$(echo $line | cut -d " " -f3)
	ref=$(echo $line | cut -d " " -f1)
	
	cd /mnt/scratch/quanjian/BS/Alignment/"$ref"
	sbatch --export=sample="$sample",r1="$r1",r2="$r2",ref="$ref" --output=/mnt/home/quanjian/quan/BS/log/Alignment/fq2gvcf/"$sample"_Ref"$ref".fq2gvcf2.sbatch.log --error=/mnt/home/quanjian/quan/BS/log/Alignment/fq2gvcf/"$sample"_Ref"$ref".fq2gvcf2.sbatch.err /mnt/ufs18/home-081/quanjian/quan/BS/script/bismark_index_alignment/BismarkAlignment.batch

done < /mnt/home/quanjian/quan/BS/runfile/BSCutadaptTrimmedData.txt


#2)BismarkDeduplication.batch
cd /mnt/home/quanjian/quan/BS/runfile
while read line
do
sample=$(echo $line|cut -d " " -f1)
ref=$(echo $line | cut -d " " -f3)

sbatch --export=sample="$sample",ref="$ref" --output=/mnt/home/quanjian/quan/BS/log/Deduplication/"$sample"_Ref"$ref".sbatch.log --error=/mnt/home/quanjian/quan/BS/log/Deduplication/"$sample"_Ref"$ref".sbatch.err /mnt/home/quanjian/quan/BS/script/bismark_index_alignment/BismarkDeduplication.batch
done < /mnt/home/quanjian/quan/BS/runfile/BSAlignedSamples.txt


#---------------------------------------------------------------------------------
#                            2. Assign reads to parental genome                  #
#---------------------------------------------------------------------------------
#0)methylation count程序输入文件的构建 
#
paste <(awk '$2 ~ /S2$/' BSDeduplicated.txt|awk '{print substr($0,1,8)}') <(awk '$2 ~ /S2$/' BSDeduplicated.txt|cut -f2) <(awk '$2 ~ /B2$/' BSDeduplicated.txt|cut -f2)|tr ' ' '\t' > Samples.methylCounts.txt

#1)父本、母本的甲基化位点统计(p1为母本，p2为父本)
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

sbatch --mem=256G --export=p1bed="$p1bed",p2bed="$p2bed",dir="$dir"/"$sample",env="$env",p1bam="$p1bam",p2bam="$p2bam" --out="$dir"/"$sample"/methyl.out --error="$dir"/"$sample"/methyl.err /mnt/ufs18/rs-015/qgg/software/scripts/methylCounts.sbatch
done < <(grep -f /mnt/scratch/quanjian/BS/err.txt /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt)


#2)合并p1.site.meth.count和p2.site.meth.count文件中各自相同的c位点，
#且合并父本、母本的甲基化位点统计(p1为母本，p2为父本)
while read line
do
	sample=$(echo $line|cut -d ' ' -f1)
	sbatch --export=sample="$sample" /mnt/ufs18/home-081/quanjian/quan/BS/script/p1_p2.site.meth.count.join.batch
done < /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt


#3)将不同染色体的C碱基位点分别合并，然后sort|uniq，最终再将不同的染色体合并起来
seq 3 18|while read line
do
	chr=$(echo $line)
	sbatch --export=chr="$chr" /mnt/ufs18/home-081/quanjian/quan/BS/script/p1_p2.site.meth.count.join.chrUnion.batch
done


#4)将所有样本的p1_p2.site.meth.count.join.groupby与union文件join到一起
while read line
do
sample=$(echo $line|cut -d ' ' -f1)
sbatch --export=sample="$sample" /mnt/ufs18/home-081/quanjian/quan/BS/script/p1_p2.site.meth.count.join.chrUnion.batch
done < /mnt/home/quanjian/quan/BS/runfile/Samples.methylCounts.txt


#5)合并所有样本的p1_p2.site.meth.count.join.union.join文件已获得全基因组样本的甲基化矩阵
for tissue in B L M P
do
	tissue=$(echo $tissue)
	sbatch --export=tissue="$tissue" /mnt/ufs18/home-081/quanjian/quan/BS/script/Union.allChr.Csites.Num.batch
done



#-----------------------------------------------------------------------------------
#    3. Compare read number that were mapped to paternal and maternal genome       
#-----------------------------------------------------------------------------------
#1)统计全基因组Promoter, Utr5, 1stExon, geneBody, intron, utr3的父本、母本甲基化水平
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


##2)获取POE基因与AGE基因在gene body、1st exon 以及 gene promoter 的父、母本甲基化水平
#1. POE基因列表
#windows系统的换行符与Unix的换行符不一致，不然后续查找无输出，需要先进行换行符格式转化
#
cd /mnt/home/quanjian/resources
dos2unix POE.finalList.txt
cat ~/quan/BS/runfile/Samples.methylCounts.txt|cut -f1|sort|uniq|while read id
do
	grep -f /mnt/ufs18/home-081/quanjian/resources/POE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/gene1stExon.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/POE.gene1stExon.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/POE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/geneBody.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/POE.geneBody.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/POE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/genePromoter.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/POE.genePromoter.bed
done 


#2. AGE基因列表
dos2unix AGE.finalList.txt
cat ~/quan/BS/runfile/Samples.methylCounts.txt|cut -f1|sort|uniq|while read id
do
	grep -f /mnt/ufs18/home-081/quanjian/resources/AGE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/gene1stExon.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/AGE.gene1stExon.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/AGE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/geneBody.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/AGE.geneBody.bed
	
	grep -f /mnt/ufs18/home-081/quanjian/resources/AGE.finalList.txt /mnt/scratch/quanjian/BS/methylCounts/$id/genePromoter.parentalCom.meth_unMeth.bed > /mnt/scratch/quanjian/BS/methylCounts/$id/AGE.genePromoter.bed
done 


###3.统计POE与AGE基因列表
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


###构建POE与AGE的union list 甲基化矩阵
realpath */AGE.genePromoter.bed|cut -d '/' -f8|while read sample
do
	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqPOE.geneBody.list.txt) <(sort -k1,1 ./$sample/POE.geneBody.bed) |tr ' ' '\t' > ./$sample/POE.geneBody.Union.bed

	join -a 1 -e 'O' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqPOE.genePromoter.list.txt) <(sort -k1,1 ./$sample/POE.genePromoter.bed) |tr ' ' '\t' > ./$sample/POE.genePromoter.Union.bed

	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqPOE.gene1stExon.list.txt) <(sort -k1,1 ./$sample/POE.gene1stExon.bed) |tr ' ' '\t' > ./$sample/POE.gene1stExon.Union.bed
	
	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqAGE.geneBody.list.txt) <(sort -k1,1 ./$sample/AGE.geneBody.bed) |tr ' ' '\t' > ./$sample/AGE.geneBody.Union.bed

	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqAGE.genePromoter.list.txt) <(sort -k1,1 ./$sample/AGE.genePromoter.bed) |tr ' ' '\t' > ./$sample/AGE.genePromoter.Union.bed

	join -a 1 -e '0' -o '1.1,2.2,2.3,2.4,2.5' <(sort -k1,1 allUniqAGE.gene1stExon.list.txt) <(sort -k1,1 ./$sample/AGE.gene1stExon.bed) |tr ' ' '\t' > ./$sample/AGE.gene1stExon.Union.bed
done


###