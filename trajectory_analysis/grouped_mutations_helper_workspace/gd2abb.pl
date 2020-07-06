#!/usr/bin/env perl

use strict;
use warnings;

my ($input, $output) = @ARGV;

open INPUT, "<$input";
open OUTPUT, ">$output";

my @input_lines = <INPUT>;
@input_lines = grep /^MOB/, @input_lines;

foreach my $_ (@input_lines) {
	chomp $_;
	print $_ . "\n";
	my @split_line = split /\t/, $_;
  
	$split_line[3] =~ s/-/_/;

  die "Unexpected number of items in line." if (scalar @split_line != 8);
  
  print OUTPUT join(".", (
    "REL606",
    "MOB",
    $split_line[5],
    ($split_line[6] < 0) ? "n1" : "1",
      $split_line[3],
    $split_line[4],
    $split_line[7],
    "REF",
    "1"
  )) . "\n";
}
