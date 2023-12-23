#! /usr/bin/perl -w
# assign reads and also tabulate methylated base versus non-methylated base
# ============================================================
# z     unmethylated C in CpG context
# Z     methylated C in CpG context
# x     unmethylated C in CHG context
# X     methylated C in CHG context
# h     unmethylated C in CHH context
# H     methylated C in CHH context
# u     unmethylated C in Unknown context (CN or CHN)
# U     methylated C in Unknown context (CN or CHN)
# ============================================================

use strict;
use warnings;
use List::Util qw(sum);

my @line;
my $nm1 = 0;
my $nm2 = 0;
my $snp1;
my $snp2;
my $score1 = 0;
my $score2 = 0;
my @snps1 = ();
my @snps2 = ();
my @tags1 = ();
my @tags2 = ();
my @tagname1 = ();
my @tagid1 = ();
my @tagname2 = ();
my @tagid2 = ();
my $preread = "";
my $preg = 0;
my @presnps = ();
my @thissnps = ();
my $thisg = 0;
my $thisread = "";
my $predir = 1;
my $thisdir = 1;

while (<>) {
    
    chomp $_;
    @line = split /\t/, $_;
    ($thisread, $thisdir) = split /\//, $line[0];
    
    $score1 = "-";
    $score2 = "-";
    $nm1 = "-";
    $nm2 = "-";
    @thissnps = ("-") x 3;
    
    if ($line[1] ne "-") {
	
	($nm1, $snp1) = split /:/, $line[1];
	@snps1 = split /,/, $snp1;
	@tags1 = split /,/, $line[2];
	@tagname1 = ();
	@tagid1 = ();
	$thissnps[1] = $snp1;
	
	$score1 = 0;
	for (my $i = 0; $i <= $#tags1; $i++) {
	    ($tagname1[$i], $tagid1[$i]) = split(/\|\|/, $tags1[$i]);
	}
	
	for (my $i = 0; $i <= $#snps1; $i++) {
	    for (my $j = 0; $j <= $#tagid1; $j++) {
		if ($tagid1[$j] eq $snps1[$i]) {
		    $score1 -= 1;
		}
	    }
	}
	$score1 /= ($#tagid1 + 1);
    }
    
    if ($line[3] ne "-") {
	
	($nm2, $snp2) = split /:/, $line[3];
	@snps2 = split /,/, $snp2;
	@tags2 = split /,/, $line[4];
	@tagname2 = ();
	@tagid2 = ();
	$thissnps[2] = $snp2;
	
	$score2 = 0;
	for (my $i = 0; $i <= $#tags2; $i++) {
	    ($tagname2[$i], $tagid2[$i]) = split(/\|\|/, $tags2[$i]);
	}
	
	for (my $i = 0; $i <= $#snps2; $i++) {
	    for (my $j = 0; $j <= $#tagid2; $j++) {
		if ($tagid2[$j] eq $snps2[$i]) {
		    $score2 -= 1;
		}
	    }
	}
	$score2 /= ($#tagid2 + 1);
    }
    
    # rule to check direction, could be 0 (not sure), 1, or 2.
    
    if ($score1 eq "-") {
	
	if ($score2 < 0) {
	    
	    $thisg = 1;
	    
	} elsif ($score2 == 0 && $nm2 <= 3) {
	    
	    $thisg = 2;
	    
	} else {
	    
	    $thisg = 0;
	    
	}
	
    } elsif ($score2 eq "-") {
	
	if ($score1 < 0) {
	    
	    $thisg = 2;
	    
	} elsif ($score1 == 0 && $nm1 <= 3) {
	    
	    $thisg = 1;
	    
	} else {
	    
	    $thisg = 0;
	    
	}
	
    } else {
	
	if ($score1 < $score2) {
	    
	    $thisg = 2;
	    
	} elsif ($score1 > $score2) {
	    
	    $thisg = 1;
	    
	} elsif ($score1 == $score2) {
	    
	    if ($nm1 > $nm2) {
		
		$thisg = 2;
		
	    } elsif ($nm1 < $nm2) {
		
		$thisg = 1;
		
	    } else {
		
		$thisg = 0;
		
	    }
	    
	}
	
    }
    
    
    if ($preread eq $thisread) {
	
	if ($thisg == $preg) {
	    
	    print $preread, "/", $predir, "\t", $preg, "\t", $presnps[$preg], "\n";
	    
	} elsif ($thisg != $preg) {
	    
	    if ($thisg == 0) {
		
		print $preread, "/", $predir, "\t", $preg, "\t", $presnps[$preg], "\n";
		$thisg = $preg;
		
	    } elsif ($preg == 0) {
		
		print $preread, "/", $predir, "\t", $thisg, "\t", $presnps[$thisg], "\n";
		
	    } else {
		
		print $preread, "/", $predir, "\t", "0", "\t", "-", "\n";
		$thisg = 0;
		
	    }
	    
	}
	
    } else {
	
	if ($preread ne "") {
	    
	    print $preread, "/", $predir, "\t", $preg, "\t", $presnps[$preg], "\n";
	    
	}
	
    }
    
    $preg = $thisg;
    $preread = $thisread;
    $predir = $thisdir;
    @presnps = @thissnps;
    
}

print $preread, "/", $predir, "\t", $preg, "\t", $presnps[$preg], "\n";

