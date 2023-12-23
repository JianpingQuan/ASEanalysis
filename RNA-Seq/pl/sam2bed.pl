#! /usr/bin/perl -w

use strict;
use warnings;
use List::Util qw(sum);

my @line;
my $read = 0;
my $strand;
my $nm = 0;
my $thismd;
my $thiscigar;
my $chr;
my $thisstart;
my @parse;
my $blocksum;
my @var;

while (<>) {
	
	chomp $_;
	@line = split /\t/, $_;
	
	if ($line[1] & 0x40) {
		
		$read = 1;
		
	} else {
		
		$read = 2;
		
	}
	
	
	if ($line[1] & 0x10) {
		$strand = "-";
	} else {
		$strand = "+";
	}
	
	$thiscigar = $line[5];
	$thisstart = $line[3] - 1;
	$chr = $line[2];
	
	if ($_ =~ m/NM:i:(\d+)\s/) {
		
		$nm = $1;
		
		if ($nm > 0) {
			
			if ($_ =~ m/MD:Z:(.*?)\s/) {
				
				$thismd = $1;
				
			}
			
		} else {
			
			$thismd = "";
			
		}
				
	}
		
	@parse = processCIGAR($thiscigar, $thismd, $line[9], $thisstart);
	$blocksum = $parse[0][$#{$parse[0]}] + $parse[1][$#{$parse[1]}];
	
	if ($#{$parse[2]} < 0) {
		
		@var = ("");
		
	} else {
	
		for (my $i = 0; $i <= $#{$parse[2]}; $i++) {
		
			push(@var, $chr."_".$parse[2][$i]."_".$parse[3][$i]."_".$parse[4][$i]."_".$parse[5][$i]); 
		
		}
		
	}
	
	print $line[2], "\t", $thisstart, "\t", $thisstart + $blocksum, "\t",
		$line[0], "/", $read, "\t.\t", $strand, "\t", $thisstart, "\t", $thisstart + $blocksum,
		"\t", $nm, ":", join(",", @var), "\t", $#{$parse[0]} + 1, "\t", join(",", @{$parse[1]}), "\t", join(",", @{$parse[0]}), "\n";
	
	@var = ();
	
}

sub processCIGAR {
	
	my $cigar = $_[0];
	my $md = $_[1];
	my $seq = $_[2];
	my $start = $_[3];
	my @blocksize = ();
	my @blockstart = ();
	my @readblocksize = ();
	my @readblockstart = ();
	my @varpos = ();
	my @varbase = ();
	my @vartype = ();
	my @refbase = ();
	my $thisoff = 0;
	my $preblocksize = 0;
	my $preblockstart = 0;
	my $prereadblocksize = 0;
	my $prereadblockstart = 0;
	my $readoffset = 0;
	
	# get rid of soft clip at the ends, the 
	if ($cigar =~ m/(\d+)S$/) {
		
		$cigar = substr $cigar, 0, length($cigar) - length($1) - 1;
		$seq = substr $seq, 0, length($seq) - $1;
		
	}
  
  if ($cigar =~ m/(\d+)H$/) {
    
    $cigar = substr $cigar, 0, length($cigar) - length($1) - 1;
    
  }
  
  if ($cigar =~ m/^(\d+)H/) {
    
		$cigar = substr $cigar, length($1) + 1;
    
  }
  
	
	if ($cigar =~ m/^(\d+)S/) {
		
		$seq = substr $seq, $1;
		$cigar = substr $cigar, length($1) + 1;
		
	}

	while (!($cigar =~ m/^$/)) {
		
		if ($cigar =~ m/^(\d+)M/) {
			
	 		# match, extend block;
			$preblocksize += $1;
			$cigar = substr $cigar, length($1) + 1;
			$prereadblocksize += $1;
			next;

		} elsif ($cigar =~ m/^(\d+)I/) {

			# insertion, also extend block in read but not genome
			
			# get insertion sequence and position
			push(@varpos, $start + $preblockstart + $preblocksize);
			push(@vartype, "INS");
			push(@varbase, substr($seq, $prereadblockstart + $prereadblocksize - $readoffset, $1));
			push(@refbase, "-" x $1);
			
			# get rid of inserted sequence
			substr $seq, $prereadblockstart + $prereadblocksize - $readoffset, $1, "";
			$readoffset += $1;
			$prereadblocksize += $1;
			$cigar = substr $cigar, length($1) + 1;
			
			next;

		} elsif ($cigar =~ m/^(\d+)D/) {

			# deletion, extend block in genome but not read
			
			# add sequence
			substr $seq, $prereadblockstart + $prereadblocksize - $readoffset, 0, "-" x $1;
			$readoffset -= $1;
			
			$preblocksize += $1;
			$cigar = substr $cigar, length($1) + 1;
			next;

		} elsif ($cigar =~ m/^(\d+)N/) {

			# when seeing N, put everything to the blocks
			$cigar = substr $cigar, length($1) + 1;
			push(@readblocksize, $prereadblocksize);
			push(@readblockstart, $prereadblockstart);
			$prereadblockstart = $prereadblockstart + $prereadblocksize;
			$prereadblocksize = 0;
			push(@blocksize, $preblocksize);
			push(@blockstart, $preblockstart);
			$preblockstart = $preblockstart + $preblocksize + $1;
			$preblocksize = 0;
			next;

		}

	}

	push(@readblocksize, $prereadblocksize);
	push(@readblockstart, $prereadblockstart);
	push(@blocksize, $preblocksize);
	push(@blockstart, $preblockstart);
	
	# now process the MD tag
	
	if ($md ne "") {

		my $preend = 0; # preend is the end of the previous MD segment
		my $base;
		my $pos;
		my $thisint = 0;
	
		# calculate cumulative sum
		# this is the cumulative length for blocks (reference)
		my @cumsum = (0);
	
		for (my $i = 1; $i <= $#blocksize; $i++) {
		
			$cumsum[$i] = $cumsum[$i - 1] + $blocksize[$i - 1];
		
		}

		while (!($md =~ m/^$/) && !($md =~ m/^\d+$/)) {

			if ($md =~ m/^(\d+)([ATCGN\^]+)/) {

				$pos = $1;
				$base = $2;

				for (my $i = $thisint; $i <= $#cumsum; $i++) {

					if ($preend + $pos >= $cumsum[$i] && $preend + $pos < $cumsum[$i] + $blocksize[$i]) {
					
						$thisint = $i;
						last;

					}

				}
						
				# deletion
				if ($base =~ m/^\^/) {
				
					push(@vartype, "DEL");
					push(@varbase, substr($seq, $preend + $pos, length($base) - 1));
					push(@refbase, substr($base, 1));
					push(@varpos, $start + $blockstart[$thisint] + $preend + $pos - $cumsum[$thisint]);
					$preend += $pos + length($base) - 1;
				
				} else {
				
					if ($pos == 0) {
									
						if ($#refbase >= 0) {
						
							if ($vartype[$#vartype] eq "SNP") { 
						
								$refbase[$#refbase] .= $base;
								$varbase[$#varbase] .= substr $seq, $preend + $pos, length($base);
						
							} else {
							
								push(@vartype, "SNP");
								push(@varbase, substr($seq, $preend + $pos, length($base)));
								push(@refbase, $base);
								push(@varpos, $start + $blockstart[$thisint] + $preend + $pos - $cumsum[$thisint]);
							
							}
						
						} else {
						
							push(@vartype, "SNP");
							push(@varbase, substr($seq, $preend + $pos, length($base)));
							push(@refbase, $base);
							push(@varpos, $start + $blockstart[$thisint] + $preend + $pos - $cumsum[$thisint]);
						
						}
					
					} else {
					
						push(@vartype, "SNP");
						push(@varbase, substr($seq, $preend + $pos, length($base)));
						push(@refbase, $base);
						push(@varpos, $start + $blockstart[$thisint] + $preend + $pos - $cumsum[$thisint]);
					
					}
				
					$preend += $pos + length($base);
				
				}
			
				substr $md, 0, length($pos) + length($base), "";
						
			}

		}
	
	}
	
	# print $seq, "\n";
	# print "varpos:", join(",", @varpos), "\n";
	# print "vartype:", join(",", @vartype), "\n";
	# print "refbase:", join(",", @refbase), "\n";
	# print "varbase:", join(",", @varbase), "\n";
	return(\@blockstart, \@blocksize, \@varpos, \@vartype, \@refbase, \@varbase);
	
}


