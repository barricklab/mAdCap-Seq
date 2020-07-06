#!/usr/bin/env perl

use strict;
use warnings;

my ($input, $output) = @ARGV;

open INPUT, "<$input";
open OUTPUT, ">$output";

print OUTPUT "#=GENOME_DIFF 1.0\n";

while (<INPUT>) {
	chomp $_;
	print $_ . "\n";
	my @split_line = split /,|\t/, $_;
	my @split_name = split /\./, $split_line[1];
	
	#print "$#split_line" . "\n";
	#print "@split_name" . "\n";
	
	next if ($#split_name < 7);

	if ($split_name[1] eq "JC") {
		my $ref1 = $split_name[0];
		$ref1 =~ s/_/-/;
		my $ref2 = $split_name[4];
		$ref2 =~ s/_/-/;
		print OUTPUT join("\t", "JC", ".", ".", $ref1, $split_name[2], (($split_name[3] eq "n1") ? 
		"-1" : "1"), $ref2, $split_name[5], (($split_name[6] eq "n1") ? 
		"-1" : "1") ,"0") . "\n"; 
	}
	elsif ($split_name[1] eq "RA") {
		my $ref = $split_name[0];
		$ref =~ s/_/-/;
		my $pos = $split_name[2];
		my $from_base = $split_name[3];
		my $to_base = $split_name[4];
		my $insert_pos = $split_name[5];
		print OUTPUT join("\t", "RA", ".", ".", $ref, $pos, $insert_pos, $from_base, $to_base, "frequency=1") . "\n"; 
	}
}
