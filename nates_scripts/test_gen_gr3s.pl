#!/usr/bin/perl

use strict;
use CORIE;

if ($#ARGV!=1) {
    print "<source hgrid> <out_dir>\n";
    exit 0;
}


my $hgrid = $ARGV[0];
my $out_dir = $ARGV[1];

open(HGRID, $hgrid) or die "could not open $hgrid for reading\n";
my @in = <HGRID>;
close(HGRID);

my ($ne, $nn) = split(/\s+/, $in[1]);
my $i;
my @vals;


open(TMP, ">$out_dir/tmp.txt") or die "could not open $out_dir/tmp.txt for writing\n";
for ($i=2; $i<=$nn+1; ++$i) {
    print(TMP $in[$i]);
}
close(TMP);
open(SPCS_IN, "| perl /home/users/hyde/bin/spcs.pl");
print SPCS_IN "tmp.txt\n";
print SPCS_IN "ll.tmp\n";
print SPCS_IN "2\n";
print SPCS_IN "1\n";
print SPCS_IN "8\n";
print SPCS_IN "1\n";
print SPCS_IN "3601\n";
print SPCS_IN "wo\n";
print SPCS_IN "2\n";
close SPCS_IN;

open(TMP, "ll.tmp") or die "could not open ll.tmp, an intermediate file in hgrid.ll generation\n";
my @ll = <TMP>;
close(TMP);
open(HGRIDLL, ">$out_dir/hgrid.ll") or die "could not opne hgrid.ll for writing";
print(HGRIDLL "hgrid.ll\n$ne $nn\n");
for ($i=2; $i<=$nn+1; ++$i) {
    print(HGRIDLL $ll[$i]);
}
for ($i=0; $i<$ne; ++$i) {
    print HGRIDLL $in[$nn+2+$i];
}
close(HGRIDLL);
