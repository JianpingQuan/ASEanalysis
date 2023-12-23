#!/bin/bash

# Job name:
#SBATCH --job-name=atac-counts
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=2
##SBATCH --qos=scavenger
# Memory per node:
#SBATCH --mem=32G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

# ============================================================

source "$env"

echo -e "arguments given:\ndir=$dir\np1bed=$p1bed\np2bed=$p2bed\nenv=$env\np1bam=$p1bam\np2bam=$p2bam"

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

# 4. get reads and lift them
# ============================================================

$BEDTOOLS intersect -a $p1bam -b $dir/p1.informative.bed -split | $SAMTOOLS view -q 50 | perl $SAM2BED | $BEDTOOLS intersect -a - -b $dir/p1.informative.bed -split -wa -wb -bed | sed 's/\//\t/' | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' | sort -k1,1 -k5,5 | uniq | join $dir/p1.read - | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5}' > $dir/p1.read.bed &
$BEDTOOLS intersect -a $p2bam -b $dir/p2.informative.bed -split | $SAMTOOLS view -q 50 | perl $SAM2BED | $BEDTOOLS intersect -a - -b $dir/p2.informative.bed -split -wa -wb -bed | sed 's/\//\t/' | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$5}' | sort -k1,1 -k5,5 | uniq | join $dir/p2.read - | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5}' > $dir/p2.read.bed &
wait

# 5. now lift reads back to the original reference, using the lift.bed
# ============================================================

awk 'length($4) != length($5)' $p1bed | perl $OFFSET > $dir/p1.offset.bed &
awk 'length($4) != length($5)' $p2bed | perl $OFFSET > $dir/p2.offset.bed &
wait

# 6. lift by intersect
# ============================================================

awk '{print $1"\t"$2"\t"$2+1"\t"$3-$2"\t"$4"\t"$5}' $dir/p1.read.bed | $BEDTOOLS intersect -a $dir/p1.offset.bed -b - -wa -wb -bed | awk '{print $5"\t"$6 + $4"\t"$6 + $4 + $8"\t"$9"/"$10}' > $dir/p1.ref.bed &
awk '{print $1"\t"$2"\t"$2+1"\t"$3-$2"\t"$4"\t"$5}' $dir/p2.read.bed | $BEDTOOLS intersect -a $dir/p2.offset.bed -b - -wa -wb -bed | awk '{print $5"\t"$6 + $4"\t"$6 + $4 + $8"\t"$9"/"$10}' > $dir/p2.ref.bed &
wait

# clean up
# ============================================================

rm $dir/p*.informative.bed $dir/p*.var.read* $dir/read.assign.txt $dir/p*.read.bed $dir/p*.offset.bed

echo $(date) done.main.process
