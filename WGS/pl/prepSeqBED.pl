#! /usr/bin/perl -w
# process VCF file to perform
# 1) call consensus
# 2) remove clusters of snps that overlap with deletions
# 3) imputation of missing genotypes based on consensus

use Getopt::Long;
use strict;
use warnings;
use List::Util qw(first sum);

my $vcf;
my $geno = 0.8; # non-missing rate
my $freq = 0.8; # frequency of the major allele
my $group1; # sample name for group1, comma separated, must come from the same breed, regardless of sire/dam
my $group2; # sample name for group2, comma separated, must come from the same breed, regardless of sire/dam

# output also a cross file bed to indicate informative sites and the genotypes expected
# ============================================================

GetOptions("vcf=s" => \$vcf,
           "geno=f" => \$geno,
					 "freq=f" => \$freq,
           "group1=s" => \$group1,
           "group2=s" => \$group2);

my @samples;
my @group1ids = split /,/, $group1;
my @group2ids = split /,/, $group2;

open VCF, "<$vcf";

while ( <VCF> ) {
  
  chomp $_;
  if ($_ =~ m/^\#\#/) {
    next;
  } elsif ($_ =~ m/^\#CHROM/) {
    my @line = split /\t/, $_;
    @samples = @line[9..$#line];
    last;
  }
  
}

# process the samples based on groups
my @group1idx;
my @group2idx;
my $idx;

for (my $i = 0; $i <= $#group1ids; $i++) {
  
  $idx = first { $samples[$_] eq $group1ids[$i] } 0..$#samples;
  if (defined($idx)) {
    $group1idx[$i] = $idx;
  } else {
    $group1idx[$i] = -1;
  }
  
}

for (my $i = 0; $i <= $#group2ids; $i++) {
  
  $idx = first { $samples[$_] eq $group2ids[$i] } 0..$#samples;
  if (defined($idx)) {
    $group2idx[$i] = $idx;
  } else {
    $group2idx[$i] = -1;
  }
  
}

my $pre_chr = 0;
my $pre_start = 0;
my $pre_end = 0;
my @pre_line = ();

print "chr\tstart\tend\tref\tid\t", join("\t", @group1ids), "\t", join("\t", @group2ids), "\n";

my $firstline = <VCF>;
chomp $firstline;
my @firstlines = split /\t/, $firstline;

$pre_chr = $firstlines[0];
$pre_start = $firstlines[1] - 1;
$pre_end = $firstlines[1] - 1 + length($firstlines[3]);
@pre_line = @firstlines;

while (<VCF>) {
  
  # filter VCF for overlapping variants with deletions
  # only keep the deletion when overlapping
  
  chomp $_;
  my @line = split /\t/, $_;
  
  if (($line[0] ne $pre_chr || $line[1] - 1 >= $pre_end) && !($line[4] =~ m/\*/)) {
    
      print $pre_line[0], "\t", $pre_start, "\t", $pre_start + length($pre_line[3]), "\t", $pre_line[3], "\t", $pre_line[0], "_", $pre_start, "_", $pre_end, "\t";
      print join("\t", processLine($geno, $freq, join(",", @group1idx), @pre_line)), "\t";
      print join("\t", processLine($geno, $freq, join(",", @group2idx), @pre_line)), "\n";
      @pre_line = @line;
      $pre_start = $line[1] - 1;
      $pre_end = $line[1] - 1 + length($line[3]);
      $pre_chr = $line[0];
    
  } else {
    
    if ($line[1] - 1 + length($line[3]) > $pre_end) {
      
      $pre_end = $line[1] - 1 + length($line[3]);
      
    }
    
  }

}

close VCF;

print $pre_line[0], "\t", $pre_start, "\t", $pre_start + length($pre_line[3]), "\t", $pre_line[3], "\t", $pre_line[0], "_", $pre_start, "_", $pre_end, "\t";
print join("\t", processLine($geno, $freq, join(",", @group1idx), @pre_line)), "\t";
print join("\t", processLine($geno, $freq, join(",", @group2idx), @pre_line)), "\n";

# process lines
# rules
# 1) impute homo major based on $freq and $geno, if don't meet criteria, non-informative
# 2) for het snp, use the first allele for sequence, flag as non-informative
# use upper case for informative, lower for non-informative

sub processLine {
  
  my $geno = $_[0];
	my $freq = $_[1];
  my @groupidx = split /,/, $_[2];
  my $chr = $_[3];
  my $pos = $_[4];
  my $ref = $_[6];
  my $alt = $_[7];
  my @alleles = split /,/, $ref.",".$alt;
  my @samples = @_[12..$#_];
  my %allele_count = ();
  my @groupgeno = ();
  my @allele_sort = ();
	my $total_count = 0;
	my $miss_count = 0;
	my $nomiss = 0;
	my $total_allele_count = 0;
  
  # process genotype
  
  for (my $i = 0; $i <= $#groupidx; $i++) {
    
    if ($groupidx[$i] >= 0) {
      
      my @this_sample_info = split /:/, $samples[$groupidx[$i]];
      
      if ($this_sample_info[0] eq "./." || $this_sample_info[0] eq ".|." || $this_sample_info[0] eq ".") {
        
        push(@groupgeno, ".");
				$miss_count += 1;
				$total_count += 1;        
        
      } else {
				
				$total_count += 1;
        
        my @this_alleles = split /[\/\|]/, $this_sample_info[0];
        
        for (my $j = 0; $j <= $#this_alleles; $j++) {
          
          if (defined($allele_count{$alleles[$this_alleles[$j]]})) {
            
            $allele_count{$alleles[$this_alleles[$j]]} += 1;
            
          } else {
            
            $allele_count{$alleles[$this_alleles[$j]]} = 1;
            
          }
					
					$total_allele_count += 1;
          
        }
        
        if ($#this_alleles == 0) {
          
          push(@groupgeno, $alleles[$this_alleles[0]]);
          
        } else {
        
          if ($this_alleles[0] eq $this_alleles[1]) {
          
            push(@groupgeno, $alleles[$this_alleles[0]]);
        
          } else {
        
            push(@groupgeno, lc($alleles[$this_alleles[0]]));
        
          }
          
        }
        
      }
      
    } else {
      
      push(@groupgeno, ".");
      
    }
    
  }
	
  # imputation or output genotype
  @allele_sort = sort { $allele_count{$b} <=> $allele_count{$a} } keys %allele_count;
	$nomiss = 1 - $miss_count/$total_count;
	
	for (my $i = 0; $i <= $#groupgeno; $i++) {
		
    if ($groupgeno[$i] eq ".") {
			
			if ($#allele_sort >= 0) {
							
				if ($allele_count{$allele_sort[0]}/$total_allele_count >= $freq && $nomiss >= $geno) {
					
					$groupgeno[$i] = $allele_sort[0];
					
				} else {
					
					$groupgeno[$i] = lc($allele_sort[0]);
					
				}
				
			} else {
				
				$groupgeno[$i] = lc($ref);
				
			}
			
		}
		
	}
			
  
  return(@groupgeno);
  
}

