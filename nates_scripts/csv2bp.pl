#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

# in_file is assumed to start at day 1 (sday) and end at or after the last day (eday)

if ($#ARGV!=1) {
    print "$0 <csv_file> <bp_file>\n";
    exit 0;
}

my $in_file = $ARGV[0];
my $out_file = $ARGV[1];
my $i;
my @vals;

open(CSV, $in_file) or die "could not open $in_file for reading\n";
my @csv = <CSV>;
close(CSV);

open(BP, ">$out_file") or die "coul dnot open $out_file for writing\n";

print BP "bp file from $in_file\n";
print BP "".scalar(@csv)."\n";

for ($i=0; $i<=$#csv; ++$i) {
    @vals = split(/\s+/, $csv[$i]);
    $vals[0] =~ /(\d+\.\d*),/;
    my $x = $1;
    $vals[1] =~ /(\d+\.\d*),/;
    my $y = $1;   
    print BP "".($i+1)." $x $y\n";
}

close(BP);

