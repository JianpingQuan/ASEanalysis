#! /usr/bin/perl -w
# ============================================================

use Getopt::Long;
use strict;
use warnings;

# read bed file
# ============================================================

my %bed_start = ();
my %bed_end = ();
my %ref = ();
my %alt = ();

while (<>) {
  
    chomp $_;
    my @this_line = split /\t/, $_;
    push(@{ $bed_start{$this_line[0]} }, $this_line[1]);
    push(@{ $bed_end{$this_line[0]} }, $this_line[2]);
    push(@{ $ref{$this_line[0]} }, $this_line[3]);
    push(@{ $alt{$this_line[0]} }, $this_line[4]);
  
}

# process one chromosome at a time
# ============================================================

foreach my $chr (sort(keys(%bed_start))) {
  
    my @starts = @{ $bed_start{$chr} };
    my @ends = @{ $bed_end{$chr} };
    my @refs = @{ $ref{$chr} };
    my @alts = @{ $alt{$chr} };
    
    my @idx = sort { $starts[$a] <=> $starts[$b] } 0..$#starts;
    my @starts_sorted = @starts[@idx];
    my @ends_sorted = @ends[@idx];
    my @refs_sorted = @refs[@idx];
    my @alts_sorted = @alts[@idx];
    
    my $offset = 0;
    my $preend = 0;
    
    for (my $i = 0; $i <= $#starts_sorted; $i++) {
	
	print $chr, "\t", $preend, "\t", $starts_sorted[$i] + length($refs_sorted[$i]), "\t", $offset, "\n";
	$offset += length($alts_sorted[$i]) - length($refs_sorted[$i]);
	$preend = $ends_sorted[$i];
	
    }
    
    print $chr, "\t", $preend, "\t", $preend + 10000000, "\t", $offset, "\n";
    
}


