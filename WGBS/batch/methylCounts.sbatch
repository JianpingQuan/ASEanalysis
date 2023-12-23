#!/bin/bash

# Job name:
#SBATCH --job-name=methyl-counts
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
#SBATCH --mem=96G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=12:00:00

# ============================================================

source "$env"

echo -e "arguments given:\ndir=$dir\np1bed=$p1bed\np2bed=$p2bed\nenv=$env\np1bam=$p1bam\np2bam=$p2bam\n"

# 1. using the two lift.bed, get informative position in bed format
# the format is
# chr, start, end (relative to the lifted genome), base in lifted genome,
# original ref base, position in ref genome, snp alleles if mismatch to the
# other parent, for example, 1_111_SNP_C_T means
# when mapping to this lifted genome, reads from the other genome
# would have a C->T snp
# ============================================================

echo $(date) "begin informative beds."

paste $p1bed $p2bed | awk '$4 ~ /[ATCGN]/ && $10 ~ /[ATCGN]/ && $4 != $10 && $6 == $12' | perl -wne 'chomp $_;
  @line = split /\t/, $_;
  if (length($line[3]) == length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1], "_SNP_", $line[3], "_", $line[9], "\n";
  } elsif (length($line[3]) > length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[9]), "_DEL_";
    print substr($line[3], length($line[9])), "_", "-" x (length($line[3]) - length($line[9])), "\n";
  } else {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[3]), "_INS_";
    print "-" x (length($line[9]) - length($line[3])), "_", substr($line[9], length($line[3])), "\n"; }' > $dir/p1.informative.bed &

paste $p2bed $p1bed | awk '$4 ~ /[ATCGN]/ && $10 ~ /[ATCGN]/ && $4 != $10 && $6 == $12' | perl -wne 'chomp $_;
  @line = split /\t/, $_;
  if (length($line[3]) == length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1], "_SNP_", $line[3], "_", $line[9], "\n";
  } elsif (length($line[3]) > length($line[9])) {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[9]), "_DEL_";
    print substr($line[3], length($line[9])), "_", "-" x (length($line[3]) - length($line[9])), "\n";
  } else {
    print join("\t", @line[0..5]), "\t", $line[0], "_", $line[1] + length($line[3]), "_INS_";
    print "-" x (length($line[9]) - length($line[3])), "_", substr($line[9], length($line[3])), "\n"; }' > $dir/p2.informative.bed &
wait

echo $(date) "end informative beds."

# 2. subset bam files to output only those intersect with informative regions
# then output the read with variants in the present alignment
# it appears that mapq has nothing to do with uniqueness
# so just use all reads
# for DNA-seq, no "N" in CIGAR
# ============================================================

echo $(date) "begin var.read."

$BEDTOOLS intersect -a $p1bam -b $dir/p1.informative.bed -split | $SAMTOOLS view | perl $BSSAM2BED | $BEDTOOLS intersect -a - -b $dir/p1.informative.bed -split -wa -wb -bed | awk -F "\t" '{print $4"||"$9"\t"$18"||"$19}' | $BEDTOOLS groupby -g 1 -c 2 -o distinct | sed 's/||/\t/' > $dir/p1.var.read &
$BEDTOOLS intersect -a $p2bam -b $dir/p2.informative.bed -split | $SAMTOOLS view | perl $BSSAM2BED | $BEDTOOLS intersect -a - -b $dir/p2.informative.bed -split -wa -wb -bed | awk -F "\t" '{print $4"||"$9"\t"$18"||"$19}' | $BEDTOOLS groupby -g 1 -c 2 -o distinct | sed 's/||/\t/' > $dir/p2.var.read &
wait

echo $(date) "end var.read."

echo $(date) "begin sort var.read."
sort -k1,1 $dir/p1.var.read > $dir/p1.var.read.sorted &
sort -k1,1 $dir/p2.var.read > $dir/p2.var.read.sorted &
echo $(date) "end sort var.read."
wait

# 3. join the reads and process them to assign reads
# ============================================================

echo $(date) "begin assign read."
join -a 1 -a 2 -t $'\t' -o '0,1.2,1.3,2.2,2.3' -e '-' $dir/p1.var.read.sorted $dir/p2.var.read.sorted | perl $ASSIGNBISREADS  > $dir/read.assign.txt
# summarize assignments, the majority should be assigned becasue we
# have in the previous step subset the bam files
cut -f 2 $dir/read.assign.txt | sort | uniq -c > $dir/read.assign.sum
echo $(date) "end assign read."

# 4. get methylation status at each site
# ============================================================

echo $(date) "begin get methylation status."

awk '$2 == 1 && $3 != "-" && $3 != ""' $dir/read.assign.txt | perl -wne 'chomp $_; @line = split /\t/, $_; @mts = split /,/, $line[2]; for (my $i = 0; $i <= $#mts; $i++) { @this = split /_/, $mts[$i]; if ($this[2] eq "MT") { print $this[0], "_", $this[1], "\t", $this[5], "\n"; }}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2 -o freqdesc > $dir/p1.site &

awk '$2 == 2 && $3 != "-" && $3 != ""' $dir/read.assign.txt | perl -wne 'chomp $_; @line = split /\t/, $_; @mts = split /,/, $line[2]; for (my $i = 0; $i <= $#mts; $i++) { @this = split /_/, $mts[$i]; if ($this[2] eq "MT") { print $this[0], "_", $this[1], "\t", $this[5], "\n"; }}' | sort -k1,1 | $BEDTOOLS groupby -g 1 -c 2 -o freqdesc > $dir/p2.site &
wait
echo $(date) "end get methylation status."

# 5. convert the methylation status back to the original coordinate
# ============================================================

awk 'length($4) != length($5)' $p1bed | perl $OFFSET > $dir/p1.offset.bed &
awk 'length($4) != length($5)' $p2bed | perl $OFFSET > $dir/p2.offset.bed &
wait

# 6. There are bases that are part of the inserted sequences in either genomes
# these are non-existent in the reference genome and need to be removed
# these need to be done in both informative and uninformative regions
# ============================================================

awk 'length($4) > length($5) { print $1"\t"$2+length($5)"\t"$3 }' $p1bed > $dir/p1.inserted.bed
awk 'length($4) > length($5) { print $1"\t"$2+length($5)"\t"$3 }' $p2bed > $dir/p2.inserted.bed

# 7. count
# ============================================================

perl -wne 'chomp $_; @line = split /\t/, $_; @info = split /,/, $line[1]; $me = 0; $ume = 0; $good = 0; $type = ""; if ($#info == 0) { @count = split /:/, $info[0]; if ($count[0] =~ m/[ZXHU]/) { $good = 1; $me = $count[1]; $ume = 0; $type = $count[0]; } elsif ($count[0] =~ m/[zxhu]/) { $good = 1; $ume = $count[1]; $me = 0; $type = uc($count[0]); } } elsif ($#info == 1) { @count0 = split /:/, $info[0]; @count1 = split /:/, $info[1]; if ($count0[0] =~ m/[ZXHU]/ && $count1[0] eq lc($count0[0])) { $good = 1; $me = $count0[1]; $ume = $count1[1]; $type = $count0[0]; } elsif ($count0[0] =~ m/[zxhu]/ && $count1[0] eq uc($count0[0]) ) { $good = 1; $me = $count1[1]; $ume = $count0[1]; $type = $count1[0]; } } if ($good == 1) { @pos = split /_/, $line[0]; print $pos[0], "\t", $pos[1], "\t", $pos[1] + 1, "\t", $type, "\t", $me, "\t", $ume, "\n"; }' $dir/p1.site | $BEDTOOLS intersect -a - -b $dir/p1.inserted.bed -v -wa | $BEDTOOLS intersect -a $dir/p1.offset.bed -b - -wa -wb -bed | awk '{print $5"_"$6+$4"_"$7+$4"\t"$8":"$9":"$10}' | sort -k1,1 > $dir/p1.site.meth.count &

perl -wne 'chomp $_; @line = split /\t/, $_; @info = split /,/, $line[1]; $me = 0; $ume = 0; $good = 0; $type = ""; if ($#info == 0) { @count = split /:/, $info[0]; if ($count[0] =~ m/[ZXHU]/) { $good = 1; $me = $count[1]; $ume = 0; $type = $count[0]; } elsif ($count[0] =~ m/[zxhu]/) { $good = 1; $ume = $count[1]; $me = 0; $type = uc($count[0]); } } elsif ($#info == 1) { @count0 = split /:/, $info[0]; @count1 = split /:/, $info[1]; if ($count0[0] =~ m/[ZXHU]/ && $count1[0] eq lc($count0[0])) { $good = 1; $me = $count0[1]; $ume = $count1[1]; $type = $count0[0]; } elsif ($count0[0] =~ m/[zxhu]/ && $count1[0] eq uc($count0[0]) ) { $good = 1; $me = $count1[1]; $ume = $count0[1]; $type = $count1[0]; } } if ($good == 1) { @pos = split /_/, $line[0]; print $pos[0], "\t", $pos[1], "\t", $pos[1] + 1, "\t", $type, "\t", $me, "\t", $ume, "\n"; }' $dir/p2.site | $BEDTOOLS intersect -a - -b $dir/p2.inserted.bed -v -wa | $BEDTOOLS intersect -a $dir/p2.offset.bed -b - -wa -wb -bed | awk '{print $5"_"$6+$4"_"$7+$4"\t"$8":"$9":"$10}' | sort -k1,1 > $dir/p2.site.meth.count &

wait

# 8. join
# ============================================================

join -t $'\t' -o '0,1.2,2.2' -e '-' $dir/p1.site.meth.count $dir/p2.site.meth.count > $dir/meth.count

# clean up
# ============================================================

#rm $dir/p*.informative.bed $dir/p*.var.read* $dir/read.assign.txt $dir/p*.site $dir/p*.offset.bed $dir/p*.inserted.bed $dir/p*.site.meth.count

echo $(date) done.main.process
