#!/usr/bin/perl

use strict;
use CORIE;

if ($#ARGV!=1) {
    print "$0 <in_file> <out_file>\n";
}

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(FLUX, $infile) or die "could not open $infile for reading\n";
my @flux = <FLUX>;
close(FLUX);

my $i;

open(OUT, ">$outfile") or die "could not open $outfile for writing\n";

for ($i=0; $i<=$#flux; ++$i) {
    my @vals = split(/\s+/, $flux[$i]);
    if ($#vals ==3) {
	print OUT ($vals[1]-86400*3)." $vals[2] $vals[3]\n";
    } else {
	print OUT ($vals[0]-86400*3)." $vals[1] $vals[2]\n";
    }
}
