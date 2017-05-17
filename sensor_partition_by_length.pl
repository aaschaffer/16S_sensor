#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to partition a fasta file into two fasta files
# depending on whether the length of the sequence is within a
# prescribed range or not
# Usage: sensor_partition_by_length.pl <input fasta file>  <output file in range> <output file out of range> <output_summary> <lower bound> <upper bound>

use strict;
use warnings;

my $input_file; #input fasta file 
my $middle_output_file; #output fasta file for sequences in the middle range
my $extreme_output_file; #output fasta file for sequences outside the middle range
my $output_summary_file; #output file with classifications
my $nextline; #one line of the file
my $total_length; #total length of a sequence
my $this_length; #length of one sequence
my $lower_bound; #lower bound on lengths in the middle range
my $upper_bound; #upper bound on lengths in the middle range
my $query_name; #identifier for one query
my $line_ctr; #counter of lines in sequence file 

my @length_array; #array of sequence lengths
my $sequence_index; #index on sequences;
my $in_range; #is the length of the sequence in the prescribed range

my $TOO_LONG = -2;
my $TOO_SHORT = -1;
my $NO_MATCH = 0;
my $FULL_MATCH = 1;
my $PARTIAL = 2;

my $usage = "sensor_partition_by_length.pl <input fasta file>  <output file in range> <output file out of range> <output_summary><lower bound> <upper bound>\n";

$input_file = $ARGV[0];
$middle_output_file = $ARGV[1];
$extreme_output_file = $ARGV[2];
$output_summary_file = $ARGV[3];
$lower_bound = $ARGV[4];
$upper_bound = $ARGV[5];


if (@ARGV != 6) {
    print "The correct number of arguments is 6\n";
    die $usage;
}

open(INPUT, "<$input_file") or die "Cannot open 1 $input_file\n"; 
open(MIDDLE_OUTPUT, ">$middle_output_file") or die "Cannot open 2 $middle_output_file\n"; 
open(EXTREME_OUTPUT, ">$extreme_output_file") or die "Cannot open 3 $extreme_output_file\n"; 
open(OUTPUT_SUMMARY, ">$output_summary_file") or die "Cannot open 4 $output_summary_file\n"; 

$sequence_index  = 0;
$total_length = 0;
$line_ctr = 0;
while(defined($nextline = <INPUT>)) {
    chomp($nextline);
    $line_ctr++;
    if (($nextline =~m/^>/)) {
#	if ($total_length > 0) { # EPN: previously this line checked that $total_length > 0, 
#                                # but that meant that length 0 sequences had no output line in the OUTPUT_SUMMARY file, now they do
      if($line_ctr > 1) { 
	    $length_array[$sequence_index] = $total_length;
	    $sequence_index++;
      }
      $total_length = 0;
    }
    else {
	$this_length = length($nextline);
	# print "Adding $this_length to $total_length\n";
	$total_length += $this_length;
    }
}
$length_array[$sequence_index] = $total_length;
close(INPUT);

$sequence_index  = 0;
$in_range = 1;
open(INPUT, "<$input_file") or die "Cannot open 4 $input_file\n"; 
while(defined($nextline = <INPUT>)) {
    chomp($nextline);
    if (($nextline =~m/^>/)) {
        ($query_name) = ($nextline =~m/^>(\S+)/);
	if (($length_array[$sequence_index] >= $lower_bound) &&
            ($length_array[$sequence_index] <= $upper_bound)) {
	    $in_range = 1;
	    print OUTPUT_SUMMARY "$query_name\t0\t";
	    print OUTPUT_SUMMARY "$length_array[$sequence_index]";
	    print OUTPUT_SUMMARY "\n";
	}
	else {
	    if ($length_array[$sequence_index] < $lower_bound) {
		$in_range = $TOO_SHORT;
	    }
	    else {
		$in_range = $TOO_LONG;
	    }
	    print OUTPUT_SUMMARY "$query_name\t";
	    print OUTPUT_SUMMARY "$in_range";
	    print OUTPUT_SUMMARY "\t$length_array[$sequence_index]";
	    print OUTPUT_SUMMARY "\n";
	}
	$sequence_index++;
    }
    if (1 == $in_range) {
	print MIDDLE_OUTPUT $nextline;
	print MIDDLE_OUTPUT "\n";
    }
    else {
	print EXTREME_OUTPUT $nextline;
	print EXTREME_OUTPUT "\n";
    }
}
close (MIDDLE_OUTPUT);
close (EXTREME_OUTPUT);
close (OUTPUT_SUMMARY);
close(INPUT);
