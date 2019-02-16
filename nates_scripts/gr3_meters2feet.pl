#!/usr/bin/perl

use strict;
use CORIE;

if ($#ARGV!=2) {
    print "<source hgrid> <out hgrid> <0|1: 0 to go m to ft, 1 to go ft to m>\n";
    exit 0;
}


my $hgrid = $ARGV[0];
my $out_file = $ARGV[1];
my $factor = 1/0.3048;
if ($ARGV[2] == 1) {
    $factor = 0.3048;
}

open(HGRID, $hgrid) or die "could not open $hgrid for reading\n";
my @in = <HGRID>;
close(HGRID);

my ($ne, $nn) = split(/\s+/, $in[1]);
my $i;
my @vals;

open(OUT, ">$out_file") or die "could not open $out_file for writing\n";
if ($ARGV[2] == 1) {
    print(OUT "ft to m\n");
} else {
    print(OUT "m to ft\n");     
}

print(OUT $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print OUT "$vals[0] ".($vals[1]*$factor)." ".($vals[2]*$factor)."  $vals[3]\n";
}

for ($i=0; $i<$ne; ++$i) {
    print OUT $in[$nn+2+$i];
}

close(OUT);
