#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

if ($#ARGV!=1) {
    print "$0 <in file> <out file>\n";
    exit 0;
}

my @command;
my $source_file = $ARGV[0];
my $target_file = $ARGV[1];

my @input;
my $line;

my $i;

open(INPUT, $source_file) or die "could not open $source_file\n";
@input = <INPUT>;
close(INPUT);
open(OUTPUT, ">$target_file") or die "could not open $target_file for writing\n";

my @x;
my @y;
my @names;

foreach $line (@input) {
    if ($line =~ /^(\d+\.\d+),\s+(\d+\.\d+),\s+(\S+)/) {
        $x[$#x+1] = $1;
        $y[$#y+1] = $2;
        $names[$#names+1] = $3;
	chomp($names[$#names]);
    }
}

print OUTPUT "".scalar(@x)."\n";
for ($i=0; $i<=$#x; ++$i) {
    print OUTPUT "1\t$names[$i]\n";
    print OUTPUT "$x[$i]\t$y[$i]\t0 0\n";
}

close(OUTPUT);
