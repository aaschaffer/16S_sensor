#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to classify whether a sequence is a 16S RNA or not.
# possible classifications are: yes, no, partial, too long, or too short
# the separation of the too long an too short classes is done in an earlier
# program and passed through this program
# Usage: classification_v3.pl <length_summary file> <BLAST results file in tab-delimited format 6> <output> <identity threshold> <E-value threshold>

use strict;
use warnings;

my $length_input_file; #input file with query name, length, and preliminary classification
my $BLAST_input_file; #input file with BLAST results
my $output_file; #output file;
my @query_names; #array of tokens that are treated as the names of the queries; query names are
                 #assumed to be unique in the set
my @query_classifications; #array of query classifications, partly based on preprocessing for length
my @query_strand_match; #array of + or - indicating whether best match, if any is on the same strand 
my @query_num_alignments_best_match; #array indicating the number of local alignments to the best match  (> 1 is anomalous)
my @query_lengths; #array of query lengths;
my %query_indices; #hash mapping names to indices;
my $float_coverage; #query coverage in the best alignment to a query as a number between 0 and 1
my @int_coverage; #array of query coverages in the best alignment to a query as an integer between 0 and 100
my $num_queries;
my $query_index;
my $nextline;
my @fields;
my $this_index;

my $match_E_threshold;


my $QUERY_COLUMN = 0;
my $E_COLUMN = 10;
my $IDENTITY_COLUMN = 2;
my $QLENGTH_COLUMN = 3;
my $QUERY_START_COLUMN = 6;
my $QUERY_END_COLUMN = 7;
my $SUBJECT_START_COLUMN = 8;
my $SUBJECT_END_COLUMN = 9;

my $TOO_LONG = -2;
my $TOO_SHORT = -1;
my $NO_MATCH = 0;
my $FULL_MATCH = 1;
my $PARTIAL = 2;
my $IMPERFECT = 3;

my $lower_length = 900;
my $upper_length = 1800;
my $identity_threshold;

my @output_strings;
$output_strings[0] = "too long\tNA\tNA\tNA";
$output_strings[1] = "too short\tNA\tNA\tNA";
$output_strings[2] = "no\tNA\tNA\tNA";
$output_strings[3] = "yes\t";
$output_strings[4] = "partial\t";
$output_strings[5] = "imperfect_match\t";

my $usage = "classification_v2.pl <length_summary file> <BLAST results file in tab-delimited format 6> <output> <identity threshold> <Evalue threshold";

if (@ARGV != 5) {
	print "The correct number of arguments is 5\n";
	die $usage;
}

$length_input_file = $ARGV[0];
$BLAST_input_file = $ARGV[1];
$output_file = $ARGV[2];
$identity_threshold = $ARGV[3];
$match_E_threshold = $ARGV[4];

open(LENGTH_INPUT, "<$length_input_file") or die "Cannot open 1 $length_input_file\n"; 

$query_index  = 0;
while(defined($nextline = <LENGTH_INPUT>)) {
    chomp($nextline);
    @fields = split /\t/, $nextline;
    $query_names[$query_index] = $fields[0];
    $query_classifications[$query_index] = $fields[1];
    $query_lengths[$query_index] = $fields[2]; 
    $query_indices{$fields[0]} = $query_index;
    $query_index++;
}
$num_queries = $query_index;

close(LENGTH_INPUT); 

open(BLAST_INPUT, "<$BLAST_input_file") or die "Cannot open 2 $BLAST_input_file\n"; 
open(OUTPUT, ">$output_file") or die "Cannot open 3 $output_file\n";





while(defined($nextline = <BLAST_INPUT>)) {
    chomp($nextline);
    @fields = split /\t/, $nextline;
    if (defined($query_indices{$fields[$QUERY_COLUMN]})) {
	$this_index = $query_indices{$fields[$QUERY_COLUMN]};
	$query_num_alignments_best_match[$this_index]++; 
	if ((($fields[$E_COLUMN] == 0.0) ||
	     ($fields[$E_COLUMN] < $match_E_threshold)) &&
	    ($fields[$IDENTITY_COLUMN] > $identity_threshold)) {
	    $query_classifications[$this_index] = $FULL_MATCH;
	    $float_coverage = ($fields[$QUERY_END_COLUMN] - $fields[$QUERY_START_COLUMN])/$query_lengths[$this_index];
	    $int_coverage[$this_index] = int((100 * $float_coverage)+0.5);
#	    if (($fields[$QLENGTH_COLUMN] < $lower_length) ||
#		($fields[$QLENGTH_COLUMN] > $upper_length)) {
#		$query_classifications[$this_index] = $PARTIAL;
#	    }
	}
	else {
	    if (!($query_classifications[$this_index] eq $FULL_MATCH) &&
		!($query_classifications[$this_index] eq $PARTIAL) && !($query_classifications[$this_index] eq $IMPERFECT)) {
		$query_classifications[$this_index] = $IMPERFECT;
		$float_coverage = ($fields[$QUERY_END_COLUMN] - $fields[$QUERY_START_COLUMN])/$query_lengths[$this_index];
		$int_coverage[$this_index] = int((100 * $float_coverage)+0.5);
	    }
	}
	if ($fields[$SUBJECT_START_COLUMN] < $fields[$SUBJECT_END_COLUMN]) {
	    if ((!defined($query_strand_match[$this_index])) || ("plus" eq $query_strand_match[$this_index])) {
		$query_strand_match[$this_index] = "plus";
	    }
	    else {
		$query_strand_match[$this_index] = "mixed";
	    }
	}
	else {
	    if ((!defined($query_strand_match[$this_index])) || ("minus" eq $query_strand_match[$this_index])) {
		$query_strand_match[$this_index] = "minus";
	    }
	    else {
		$query_strand_match[$this_index] = "mixed";
	    }
	}
    }
    else {
	printf "Cannot find index for query $fields[0]\n";
    }

}

for($query_index = 0; $query_index < $num_queries; $query_index++) {
    print OUTPUT "$query_names[$query_index]";
    print OUTPUT "\t";
    print OUTPUT $output_strings[$query_classifications[$query_index]+2];
    if ($query_classifications[$query_index] > 0) {
	if (defined($int_coverage[$query_index])) { 
	    print OUTPUT "$query_strand_match[$query_index]\t$query_num_alignments_best_match[$query_index]\t$int_coverage[$query_index]";
	}
	else {
	    print OUTPUT "$query_strand_match[$query_index]\t$query_num_alignments_best_match[$query_index]\tNA";
	}
    }
    print OUTPUT "\n";
}

close(OUTPUT);
