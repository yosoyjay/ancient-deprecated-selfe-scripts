#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

# in_file is assumed to start at day 1 (sday) and end at or after the last day (eday)

if ($#ARGV!=1) {
    print "$0 <bp_file> <csv_file>\n";
    exit 0;
}

my $in_file = $ARGV[0];
my $out_file = $ARGV[1];
my $i;
my @vals;

open(BP, $in_file) or die "could not open $in_file for reading\n";
my @bp = <BP>;
close(BP);

open(CSV, ">$out_file") or die "coul dnot open $out_file for writing\n";

for ($i=2; $i<=$#bp; ++$i) {
    @vals = split(/\s+/, $bp[$i]);
    print CSV "$vals[1], $vals[2], 0, 0\n";
}

close(CSV);

