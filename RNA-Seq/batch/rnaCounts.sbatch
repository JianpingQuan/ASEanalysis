#!/bin/bash

# Job name:
#SBATCH --job-name=rna-counts
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=4
#
# Memory per node:
#SBATCH --mem=16G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

# ============================================================

source "$env"

echo -e "arguments given:\ndir=$dir\np1bed=$p1bed\np2bed=$p2bed\nenv=$env\np1bam=$p1bam\np2bam=$p2bam\np1ref=$p1ref\np2ref=$p2ref\np1gtf=$p1gtf\np2gtf=$p2gtf"

# 1. using the two lift.bed, get informative position in bed format
# the format is
# chr, start, end (relative to the lifted genome), base in lifted genome,
# original ref base, position in ref genome, snp alleles if mismatch to the
# other parent, for example, 1_111_SNP_C_T means
# when mapping to this lifted genome, reads from the other genome
# would have a C->T snp
# ============================================================

paste $p1bed $p2bed | awk '$4 ~ /[ATCGN]/ && $10 ~ /[ATCGN]/ && $4 != $10 && $6 == $12' | perl -wne 'chomp $_;
  @line = split /\t/, $_;
  if (length($line[3]) == length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1], "_SNP_", $line[3], "_", $line[9], "\n";
  } elsif (length($line[3]) > length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[9]), "_DEL_";
    print substr($line[3], length($line[9])), "_", "-" x (length($line[3]) - length($line[9])), "\n";
  } else {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[3]), "_INS_";
    print "-" x (length($line[9]) - length($line[3])), "_", substr($line[9], length($line[3])), "\n";
	}' > $dir/p1.informative.bed &

paste $p2bed $p1bed | awk '$4 ~ /[ATCGN]/ && $10 ~ /[ATCGN]/ && $4 != $10 && $6 == $12' | perl -wne 'chomp $_;
  @line = split /\t/, $_;
  if (length($line[3]) == length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1], "_SNP_", $line[3], "_", $line[9], "\n";
  } elsif (length($line[3]) > length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[9]), "_DEL_";
    print substr($line[3], length($line[9])), "_", "-" x (length($line[3]) - length($line[9])), "\n";
  } else {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[3]), "_INS_";
    print "-" x (length($line[9]) - length($line[3])), "_", substr($line[9], length($line[3])), "\n";
	}' > $dir/p2.informative.bed &
wait

# 2. subset bam files to output only those intersect with informative regions
# then output the read with variants in the present alignment
# only consider unique alignment "samtools view -q 50"
# ============================================================

$BEDTOOLS intersect -a $p1bam -b $dir/p1.informative.bed -split | $SAMTOOLS view -q 50 | perl $SAM2BED | $BEDTOOLS intersect -a - -b $dir/p1.informative.bed -split -wa -wb -bed | awk -F "\t" '{print $4"||"$9"\t"$18"||"$19}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2 -o distinct | sed 's/||/\t/' | sort -k1,1 > $dir/p1.var.read &
$BEDTOOLS intersect -a $p2bam -b $dir/p2.informative.bed -split | $SAMTOOLS view -q 50 | perl $SAM2BED | $BEDTOOLS intersect -a - -b $dir/p2.informative.bed -split -wa -wb -bed | awk -F "\t" '{print $4"||"$9"\t"$18"||"$19}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2 -o distinct | sed 's/||/\t/' | sort -k1,1 > $dir/p2.var.read &
wait

# 3. join the reads and process them to assign reads
# ============================================================

join -a 1 -a 2 -t $'\t' -o '0,1.2,1.3,2.2,2.3' -e '-' $dir/p1.var.read $dir/p2.var.read | perl $ASSIGNREADS  > $dir/read.assign.txt
# summarize assignments, the majority should be assigned becasue we
# have in the previous step subset the bam files
cut -f 2 $dir/read.assign.txt | sort | uniq -c > $dir/read.assign.sum

# get the names of the assigned reads, do not distinguish between 2 reads of the same pair
awk '$2 == 1 {print $1}' $dir/read.assign.txt | awk -F "/" '{print $1}' | uniq | sort | uniq > $dir/p1.read &
awk '$2 == 2 {print $1}' $dir/read.assign.txt | awk -F "/" '{print $1}' | uniq | sort | uniq > $dir/p2.read &
wait

# 4. generate constitutive, non-overlapping exon bed files
# ============================================================

$SAMTOOLS faidx $p1ref &
$SAMTOOLS faidx $p2ref &
wait

awk '$3 == "exon"' $p1gtf | perl -wne 'chomp $_; @line = split /\t/, $_; 
if ($line[8] =~ m/gene_id \"(.*?)\".*transcript_id \"(.*?)\";/)
{ print $line[0], "\t", $line[3] - 1, "\t", $line[4], "\t", $1, "\t", $2, "\t", $line[6], "\n";}' > $dir/p1.exon.bed &
awk '$3 == "exon"' $p2gtf | perl -wne 'chomp $_; @line = split /\t/, $_; 
if ($line[8] =~ m/gene_id \"(.*?)\".*transcript_id \"(.*?)\";/)
{ print $line[0], "\t", $line[3] - 1, "\t", $line[4], "\t", $1, "\t", $2, "\t", $line[6], "\n";}' > $dir/p2.exon.bed &
wait

# find constitutive exons
# ============================================================

for parent in p1 p2
do
  while read line
  do
    gene=`echo $line | awk '{print $1}'`
    txs=`echo $line | awk '{print $2}' | sed 's/,/\n/g'`
    txcount=`echo $line | awk '{print $3}'`
    strand=`echo $line | awk '{print $4}'`
    if [[ $txcount -eq 1 ]]
    then
      awk '$4 == "'$gene'"' $dir/$parent.exon.bed
	  else
      $BEDTOOLS coverage -a <(awk '$4 == "'$gene'"' $dir/$parent.exon.bed | sort -k1,1 -k2,2n | $BEDTOOLS merge -i - -c 4 -o distinct) -b <(awk '$4 == "'$gene'"' $dir/$parent.exon.bed) -d | awk '$6 == "'$txcount'"' | awk '{print $1"\t"$2+$5-1"\t"$2+$5}' | $BEDTOOLS merge -i - | awk '{print $0"\t'$gene'\t.\t'$strand'"}'
    fi
  done < <(cut -f 4,5,6 $dir/$parent.exon.bed | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2,2,3 -o distinct,count_distinct,distinct) > $dir/$parent.const.exon.bed &
done
wait

# remove overlapping part
# ============================================================

cat <(awk '$6 == "+"' $dir/p1.const.exon.bed | sort -k1,1 -k2,2 | $BEDTOOLS genomecov -bg -i stdin -g $p1ref.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3"\t.\t.\t+"}') <(awk '$6 == "-"' $dir/p1.const.exon.bed | sort -k1,1 -k2,2 | $BEDTOOLS genomecov -bg -i stdin -g $p1ref.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3"\t.\t.\t-"}')  > $dir/p1.multiCov.bed &

cat <(awk '$6 == "+"' $dir/p2.const.exon.bed | sort -k1,1 -k2,2 | $BEDTOOLS genomecov -bg -i stdin -g $p1ref.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3"\t.\t.\t+"}') <(awk '$6 == "-"' $dir/p2.const.exon.bed | sort -k1,1 -k2,2 | $BEDTOOLS genomecov -bg -i stdin -g $p1ref.fai | awk '$4 > 1 {print $1"\t"$2"\t"$3"\t.\t.\t-"}')  > $dir/p2.multiCov.bed &
wait

$BEDTOOLS subtract -a $dir/p1.const.exon.bed -b $dir/p1.multiCov.bed -s | sort -k1,1 -k2,2n > $dir/p1.const.exon.no.overlap.bed &
$BEDTOOLS subtract -a $dir/p2.const.exon.bed -b $dir/p2.multiCov.bed -s | sort -k1,1 -k2,2n > $dir/p2.const.exon.no.overlap.bed &
wait

# 5. assign fragments for paired end reads and assign to genes
# a read must be 80% fall within constitutive exons
# so need a read length data
# ============================================================

$BEDTOOLS intersect -a $p1bam -b $dir/p1.informative.bed -split | $SAMTOOLS view -b -q 50 | $BEDTOOLS bamtobed -i - -split | awk '{print $4"\t"$3-$2}' | $BEDTOOLS groupby -g 1 -c 2 -o sum | awk -F "/" '{print $1"\t"$2}' | sort -k1,1 | join -t $'\t' - $dir/p1.read | awk '{print $1"/"$2"\t"$3}' | sort -k1,1 > $dir/p1.read.len &
$BEDTOOLS intersect -a $p2bam -b $dir/p2.informative.bed -split | $SAMTOOLS view -b -q 50 | $BEDTOOLS bamtobed -i - -split | awk '{print $4"\t"$3-$2}' | $BEDTOOLS groupby -g 1 -c 2 -o sum | awk -F "/" '{print $1"\t"$2}' | sort -k1,1 | join -t $'\t' - $dir/p2.read | awk '{print $1"/"$2"\t"$3}' | sort -k1,1 > $dir/p2.read.len &
wait

$BEDTOOLS intersect -a $p1bam -b $dir/p1.informative.bed -split | $SAMTOOLS view -b -q 50 -f 0x40 | $BEDTOOLS intersect -a stdin -b $dir/p1.const.exon.no.overlap.bed -bed -S -wao | awk '$19 > 0 {print $4"\t"$16"\t"$19}' | sort -k1,1 | join -t $'\t' - $dir/p1.read.len | awk '{print $1"::::"$2"\t"$3"\t"$4}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2,3 -o sum,distinct | awk '$2/$3 >= 0.8 {print $1}' | awk -F "::::" '{print $1"\t"$2}' | awk -F "/" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$3}' | sort -k1,1 > $dir/p1.first.read.gene &
$BEDTOOLS intersect -a $p2bam -b $dir/p2.informative.bed -split | $SAMTOOLS view -b -q 50 -f 0x40 | $BEDTOOLS intersect -a stdin -b $dir/p2.const.exon.no.overlap.bed -bed -S -wao | awk '$19 > 0 {print $4"\t"$16"\t"$19}' | sort -k1,1 | join -t $'\t' - $dir/p2.read.len | awk '{print $1"::::"$2"\t"$3"\t"$4}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2,3 -o sum,distinct | awk '$2/$3 >= 0.8 {print $1}' | awk -F "::::" '{print $1"\t"$2}' | awk -F "/" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$3}' | sort -k1,1 > $dir/p2.first.read.gene &

$BEDTOOLS intersect -a $p1bam -b $dir/p1.informative.bed -split | $SAMTOOLS view -b -q 50 -f 0x80 | $BEDTOOLS intersect -a stdin -b $dir/p1.const.exon.no.overlap.bed -bed -s -wao | awk '$19 > 0 {print $4"\t"$16"\t"$19}' | sort -k1,1 | join -t $'\t' - $dir/p1.read.len | awk '{print $1"::::"$2"\t"$3"\t"$4}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2,3 -o sum,distinct | awk '$2/$3 >= 0.8 {print $1}' | awk -F "::::" '{print $1"\t"$2}' | awk -F "/" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$3}' | sort -k1,1 > $dir/p1.second.read.gene &
$BEDTOOLS intersect -a $p2bam -b $dir/p2.informative.bed -split | $SAMTOOLS view -b -q 50 -f 0x80 | $BEDTOOLS intersect -a stdin -b $dir/p2.const.exon.no.overlap.bed -bed -s -wao | awk '$19 > 0 {print $4"\t"$16"\t"$19}' | sort -k1,1 | join -t $'\t' - $dir/p2.read.len | awk '{print $1"::::"$2"\t"$3"\t"$4}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2,3 -o sum,distinct | awk '$2/$3 >= 0.8 {print $1}' | awk -F "::::" '{print $1"\t"$2}' | awk -F "/" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$3}' | sort -k1,1 > $dir/p2.second.read.gene &
wait

# 6. join to produce gene counts
# both reads must fall within the gene to an allele to qualify
# for the fragment to be counted
# ============================================================

join -t $'\t' -a 1 -a 2 -e '-' -o '0,1.2,2.2' $dir/p1.first.read.gene $dir/p1.second.read.gene | awk '{if ($2 == "-") { print $3; } else if ($3 == "-") { print $2; } else if ($2 == $3) { print $2; } }' | sort | uniq -c | awk '{print $2"\t"$1}' > $dir/p1.gene.count &
join -t $'\t' -a 1 -a 2 -e '-' -o '0,1.2,2.2' $dir/p2.first.read.gene $dir/p2.second.read.gene | awk '{if ($2 == "-") { print $3; } else if ($3 == "-") { print $2; } else if ($2 == $3) { print $2; } }' | sort | uniq -c | awk '{print $2"\t"$1}' > $dir/p2.gene.count &
wait

rm $dir/p*.informative.bed $dir/p*.read $dir/p*.exon.bed $dir/p*.multiCov.bed $dir/p*.const.exon.no.overlap.bed $dir/p*.read.len $dir/p*.first.read.gene $dir/p*.second.read.gene

echo $(date) done.main.process
