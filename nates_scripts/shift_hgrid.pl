#!/usr/bin/perl

use strict;
use CORIE;

if ($#ARGV!=3) {
    print "<source hgrid> <out hgrid> <x shift> <y shift>\n";
    exit 0;
}


my $hgrid = $ARGV[0];
my $out_hgrid = $ARGV[1];
my $xshift = $ARGV[2];
my $yshift = $ARGV[3];

open(HGRID, $hgrid) or die "could not open $hgrid for reading\n";
my @in = <HGRID>;
close(HGRID);

my ($ne, $nn) = split(/\s+/, $in[1]);
my $i;
my @vals;

open(SHIFT, ">$out_hgrid") or die "could not open $out_hgrid for writing\n";
print(SHIFT $in[0]);
print(SHIFT $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print SHIFT "$vals[0] ".($vals[1]+$xshift)." ".($vals[2]+$yshift)." $vals[3]\n";
}

my $curr_l = $nn+2;
for ($i=$curr_l; $i<=$#in; ++$i) {
    print SHIFT $in[$i];
}

close(SHIFT);
