#1)fastqc was used to evaluate the rowdata
cd /mnt/home/quanjian/quan/BS/runfile/
while read line
do
sample=$(echo $line | cut -d " " -f 1)
r1=$(echo $line | cut -d " " -f 2)
r2=$(echo $line | cut -d " " -f 3)
sbatch --export=sample="$sample",r1="$r1",r2="$r2" --output=/mnt/home/quanjian/quan/BS/log/fastqcCheck/"$sample".sbatch.log --error=/mnt/home/quanjian/quan/BS/log/fastqcCheck/"$sample".sbatch.err /mnt/home/quanjian/quan/BS/script/quality_Check_trim/fastqcCheck.batch
done < BSrawdata.txt


#2)cutadapt was used to cut adapter sequence and low quality bases
cd /mnt/home/quanjian/quan/BS/runfile/
while read line
do
sample=$(echo $line | cut -d " " -f 1)
r1=$(echo $line | cut -d " " -f 2)
r2=$(echo $line | cut -d " " -f 3)
sbatch --export=sample="$sample",r1="$r1",r2="$r2" --output=/mnt/home/quanjian/quan/BS/log/Cutadapt/"$sample".sbatch.log --error=/mnt/home/quanjian/quan/BS/log/Cutadapt/"$sample".sbatch.err /mnt/home/quanjian/quan/BS/script/quality_Check_trim/cutadapt.batch
done < BSrawdata.txt


#3)fastqc was used to evaluate the trimmed reads quality
cd /mnt/home/quanjian/quan/BS/runfile/
while read line
do
sample=$(echo $line | cut -d " " -f 1)
r1=$(echo $line | cut -d " " -f 2)
r2=$(echo $line | cut -d " " -f 3)
sbatch --export=sample="$sample",r1="$r1",r2="$r2" --output=/mnt/home/quanjian/quan/BS/log/fastqcCheck_trimmed/"$sample".sbatch.log --error=/mnt/home/quanjian/quan/BS/log/fastqcCheck_trimmed/"$sample".sbatch.err /mnt/home/quanjian/quan/BS/script/quality_Check_trim/fastqcCheck_trimmed.batch
done


#4)Original reference index building
sbatch BismarkOriRefPreparation.batch


#5)reads mapping based on original reference
while read line
do
sample=$(echo $line|cut -d " " -f1|cut -d '/' -f9|cut -d '_' -f1,2)
r1=$(echo $line |cut -d " " -f1)
r2=$(echo $line |cut -d " " -f2)
sbatch --export=sample=$sample,r1=$r1,r2=$r2 /mnt/home/quanjian/quan/BS/script/bismark_index_alignment/BismarkOriRefAlignment.batch
done < /mnt/home/quanjian/quan/BS/runfile/BSCutadaptTrimmedDataOriRef.txt


#6)Duplication
while read line
do
 sample=$(echo $line |cut -d " " -f1)
 sbatch --export=sample=$sample /mnt/ufs18/home-081/quanjian/quan/BS/script/bismark_index_alignment/BismarkOriRefReDup.batch
done < <(cut -f1 /mnt/home/quanjian/quan/BS/runfile/BSAlignedSamples.txt|sort|uniq)


#7)CpG site extractor
while read line
do
 sample=$(echo $line |cut -d " " -f1)
 sbatch --export=sample=$sample /mnt/home/quanjian/quan/BS/script/bismark_index_alignment/BismarkMethylationExtractor_RefOri.batch
done < <(cut -f1 /mnt/home/quanjian/quan/BS/runfile/BSAlignedSamples.txt|sort|uniq)
