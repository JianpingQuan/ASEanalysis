#! /usr/bin/perl -w

use Getopt::Long;
use strict;
use warnings;

my $fasta;
my $gtf;
my $bed;
my $gtf_out;
my $fasta_out;
my $var_out;

GetOptions("gtf=s" => \$gtf,
           "bed=s" => \$bed,
           "fasta=s" => \$fasta,
           "gtfout=s" => \$gtf_out,
           "varout=s" => \$var_out,
           "fastaout=s" => \$fasta_out);
           
# get fasta file if it's given
# ============================================================

my %seq = ();

if (defined($fasta)) {
  
  open FASTA, "<$fasta";
  my $chr;
  
  while (<FASTA>) {
    
    chomp $_;
    
    if ($_ =~ m/^>(.*?)(\s|$)/) {
     
      $chr = $1;
      $seq{$chr} = "";
      next;
      
    }
    
    $seq{$chr} .= uc($_);
    
  }
  
}

close FASTA;

# read gtf file
# ============================================================

my %gtf_line = ();
open GTF, "<$gtf";

while (<GTF>) {
  
  chomp $_;
  my @this_line = split /\t/, $_;
  push(@{ $gtf_line{$this_line[0]} }, $_);
  
}

close GTF;

# read bed file
# ============================================================

my %bed_line = ();
open BED, "<$bed";

while (<BED>) {
  
  chomp $_;
  my @this_line = split /\t/, $_;
  push(@{ $bed_line{$this_line[0]} }, $_);
  
}

close BED;

# process one chromosome at a time
# ============================================================

open NEWREF, ">$fasta_out";
print NEWREF "";
close NEWREF;

open NEWGTF, ">$gtf_out";
print NEWGTF "";
close NEWGTF;

open VAR, ">$var_out";
print VAR "";
close VAR;

my $chr_gtf_line_ref;
my $chr_bed_line_ref;
my $chr_seq;

foreach my $chr (sort(keys(%seq))) {
  
	$chr_gtf_line_ref = \@{ $gtf_line{$chr} };
  $chr_bed_line_ref = \@{ $bed_line{$chr} };
  $chr_seq = $seq{$chr};
	
  my $time = localtime();
  print STDERR "$time processing chromosome $chr.\n";
  liftChr($chr, $chr_gtf_line_ref, $chr_bed_line_ref, $chr_seq, $gtf_out, $fasta_out, $var_out);
  
}


# subroutine to process one chromosome at a time
# ============================================================

sub liftChr {
  
  # array reference for the arrays
  # ============================================================
  
  my ($chr_ref, $gtf_ref, $bed_ref, $seq_ref, $gtf_out_ref, $fasta_out_ref, $var_out_ref) = @_;
  my @chr_gtf_line = @{ $gtf_ref };
  my @chr_bed_line = @{ $bed_ref };
  my @chr_seq = $seq_ref;
  
  # get gtf information
  
  my @chr_gtf_index = ();
  my @chr_gtf_pos = ();
  my @chr_gtf_drtn = (); # 0 = start; 1 = end
  
  for (my $i = 0; $i <= $#chr_gtf_line; $i++) {
    
    my @this_gtf_line_split = split /\t/, $chr_gtf_line[$i];
    push(@chr_gtf_index, ($i, $i));
    push(@chr_gtf_pos, ($this_gtf_line_split[3] - 1, $this_gtf_line_split[4]));
    push(@chr_gtf_drtn, (0, 1));
    
  }
  
  # get bed information, also process sequences and
  # prepare for variant print outs
  
  my @chr_bed_pos = ();
  my @chr_bed_drtn = ();
  my @chr_bed_offset = ();
  my @chr_bed_index = ();
  my @chr_bed_replace = ();
	my $rep_index = 0;
  my @chr_var_start = ();
  my @chr_var_offset = ();
  my @chr_var_end = ();
  my @chr_var_ref = ();
  my @chr_var_alt = ();
  my @chr_var_info = ();
    
  for (my $i = 0; $i <= $#chr_bed_line; $i++) {
    
    my @this_bed_line_split = split /\t/, $chr_bed_line[$i];
    
    push(@chr_var_start, $this_bed_line_split[1]);
    push(@chr_var_end, $this_bed_line_split[2]);
    push(@chr_var_ref, $this_bed_line_split[3]);
    push(@chr_var_alt, $this_bed_line_split[4]);
    push(@chr_var_offset, length($this_bed_line_split[4]) - length($this_bed_line_split[3]));
    push(@chr_var_info, join("\t", @this_bed_line_split[(5..$#this_bed_line_split)]));
        
    if (length($this_bed_line_split[3]) == length($this_bed_line_split[4])) {
    
      push(@chr_bed_pos, ($this_bed_line_split[1], $this_bed_line_split[2]));
	    push(@chr_bed_replace, $this_bed_line_split[4]);
			push(@chr_bed_index, ($rep_index, $rep_index));
    	push(@chr_bed_drtn, (0, 1));
	    push(@chr_bed_offset, 0);
			$rep_index++;
			
    } elsif (length($this_bed_line_split[3]) > length($this_bed_line_split[4])) {
			
			# add part of the sequence in case they are different from ref
			push(@chr_bed_replace, $this_bed_line_split[4]);
    	push(@chr_bed_pos, ($this_bed_line_split[1], $this_bed_line_split[1] + length($this_bed_line_split[4])));
			push(@chr_bed_index, ($rep_index, $rep_index));
    	push(@chr_bed_drtn, (0, 1));
	    push(@chr_bed_offset, 0);
			$rep_index++;
			
			# add the deleted part
			push(@chr_bed_replace, "");
      push(@chr_bed_pos, ($this_bed_line_split[1] + length($this_bed_line_split[4]), $this_bed_line_split[2]));
			push(@chr_bed_index, ($rep_index, $rep_index));
    	push(@chr_bed_drtn, (0, 1));
	    push(@chr_bed_offset, length($this_bed_line_split[4]) - length($this_bed_line_split[3]));
			$rep_index++;      
    
    } else {
    	
			# add part of the sequence in case they are different from ref
			push(@chr_bed_replace, substr($this_bed_line_split[4], 0, length($this_bed_line_split[3])));
    	push(@chr_bed_pos, ($this_bed_line_split[1], $this_bed_line_split[1] + length($this_bed_line_split[3])));
			push(@chr_bed_index, ($rep_index, $rep_index));
    	push(@chr_bed_drtn, (0, 1));
	    push(@chr_bed_offset, 0);
			$rep_index++;
			
			# add the inserted part
			push(@chr_bed_replace, substr($this_bed_line_split[4],length($this_bed_line_split[3])));
      push(@chr_bed_pos, ($this_bed_line_split[2], $this_bed_line_split[2]));
			push(@chr_bed_index, ($rep_index, $rep_index));
	    push(@chr_bed_drtn, (0, 1));
	    push(@chr_bed_offset, length($this_bed_line_split[4]) - length($this_bed_line_split[3]));
			$rep_index++;
			
    }
    
  }
  
  # merge and sort
  my @chr_merge_pos = (@chr_gtf_pos, @chr_bed_pos);
  my @chr_merge_drtn = (@chr_gtf_drtn, @chr_bed_drtn);
  my @chr_merge_file = (("gtf") x ($#chr_gtf_pos + 1), ("bed") x ($#chr_bed_pos + 1));
  my @chr_merge_file_index = (@chr_gtf_index, @chr_bed_index);
	
  my @chr_merge_sort_index = sort { $chr_merge_pos[$a] <=> $chr_merge_pos[$b] || $chr_merge_file[$a] cmp $chr_merge_file[$b] || $chr_merge_file_index[$a] <=> $chr_merge_file_index[$b] } 0..$#chr_merge_pos;
  my @chr_merge_sort_pos = @chr_merge_pos[@chr_merge_sort_index];
  my @chr_merge_sort_drtn = @chr_merge_drtn[@chr_merge_sort_index];
  my @chr_merge_sort_file = @chr_merge_file[@chr_merge_sort_index];
  my @chr_merge_sort_file_index = @chr_merge_file_index[@chr_merge_sort_index];
  
  my @chr_var_sort_start_index = sort { $chr_var_start[$a] <=> $chr_var_start[$b] } 0..$#chr_var_start;
  
  # print variant set
  open VAR, ">>$var_out_ref";
  my $var_offset = 0;
  for (my $i = 0; $i <= $#chr_var_sort_start_index; $i++) {

    print VAR $chr_ref, "\t", $chr_var_start[$chr_var_sort_start_index[$i]] + $var_offset, "\t", $chr_var_start[$chr_var_sort_start_index[$i]] + $var_offset + length($chr_var_alt[$chr_var_sort_start_index[$i]]), "\t", $chr_var_alt[$chr_var_sort_start_index[$i]], "\t", $chr_var_ref[$chr_var_sort_start_index[$i]], "\t", $chr_var_info[$chr_var_sort_start_index[$i]], "\n";
    $var_offset += $chr_var_offset[$chr_var_sort_start_index[$i]];
    
  }
  close VAR;
  
  # print join(",", @chr_bed_offset), "\n";
  # print join(",", @chr_bed_pos), "\n";
  # print join(",", @chr_bed_index), "\n";
  # print join(",", @chr_merge_sort_pos), "\n";
  # print join(",", @chr_merge_sort_file), "\n";
  # print join(",", @chr_merge_sort_drtn), "\n";

  my $offset = 0;
  my $var_start = 0;
  my $new_seq = "";
  my $pre_end = 0;
  my $pre_offset = 0;

  for (my $i = 0; $i <= $#chr_merge_sort_pos; $i++) {
  
    if ($chr_merge_sort_file[$i] eq "gtf") {
      
      # print $var_start, ",", $offset, ",", $chr_merge_sort_pos[$i], ",", $chr_bed_offset[$chr_merge_sort_file_index[$i]], "\n";
    
      if ($var_start == 0 || $pre_offset == 0) {
        
        $chr_merge_sort_pos[$i] += $offset;
      
      } else { # variant initiated
              
        $chr_merge_sort_pos[$i] = $offset + $var_start;
        
      }
      
    } else {
    
      if ($var_start == 0) {
      
        $var_start = $chr_merge_sort_pos[$i];
        # handle sequence        
        $new_seq .= substr($chr_seq, $pre_end, $var_start - $pre_end);
        $pre_offset = $chr_bed_offset[$chr_merge_sort_file_index[$i]];
      
      } else {
        
        if ($chr_merge_sort_drtn[$i] != 1 && $var_start > $chr_merge_sort_pos[$i]) {
          
          print STDERR "overlapping variant at $var_start.\n";
          exit;
					
        }
      
        $offset += $chr_bed_offset[$chr_merge_sort_file_index[$i]];
        $var_start = 0;
        $pre_end = $chr_merge_sort_pos[$i];
        $new_seq .= $chr_bed_replace[$chr_merge_sort_file_index[$i]];
        
      }
    
    }
  
  }
  
  $new_seq .= substr($chr_seq, $pre_end);
  
  open NEWREF, ">>$fasta_out_ref";
  
  print NEWREF ">$chr_ref\n";
  
  for (my $i = 0; $i <= int(length($new_seq)/50); $i++) {
  
    print NEWREF substr($new_seq, $i*50, 50), "\n";
    
  }
  
  close NEWREF;
  
  open NEWGTF, ">>$gtf_out_ref";
  
  my %print_buffer = ();

  for (my $i = 0; $i <= $#chr_merge_sort_pos; $i++) {
  
    if ($chr_merge_sort_file[$i] eq "bed") {
    
      next;
    
    }
  
    if ($chr_merge_sort_drtn[$i] == 0) {
    
      $print_buffer{$chr_merge_sort_file_index[$i]} = $chr_merge_sort_pos[$i];
    
    } else {
    
      if ($chr_merge_sort_pos[$i] == $print_buffer{$chr_merge_sort_file_index[$i]}) {
      
        delete $print_buffer{$chr_merge_sort_file_index[$i]};
      
      } else {
      
        my @this_gtf_line = split /\t/, $chr_gtf_line[$chr_merge_sort_file_index[$i]];
        $this_gtf_line[3] = $print_buffer{$chr_merge_sort_file_index[$i]} + 1;
        $this_gtf_line[4] = $chr_merge_sort_pos[$i];
      
        delete $print_buffer{$chr_merge_sort_file_index[$i]};
      
        print NEWGTF join("\t", @this_gtf_line), "\n";
      
      }
    
    }
  
  }
  
  close NEWGTF;
  
}
