#!/usr/bin/perl

use strict;
use CORIE;

if ($#ARGV!=5) {
    print "$0 <source dir> <target dir> <source start> <source end> <n stacks> <target start>\n";
    exit 1;
}

my $source_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $source_start = $ARGV[2];
my $source_end = $ARGV[3];
my $nstacks = $ARGV[4];
my $target_start = $ARGV[5];

my ($i, $j, $k);

my @link_files = ("conc.63", "elev.61", "hvel.64", "pres.61", "salt.63", "srad.61", "tdff.63", "temp.63", "vdff.63", "vert.63", "wind.62", "zcor.63");

for ($i=$source_start; $i<=$source_end; ++$i) {
    my $target_num = $target_start+$i-$source_start;
    print "doing $i -> $target_num\n";

    for ($j=0; $j<=$#link_files; ++$j) {
	for ($k=0; $k<$nstacks; ++$k) {
	    my $source_fname = sprintf("%d_%04d_%s", $i, $k, $link_files[$j]);
	    my $target_fname = sprintf("%d_%04d_%s", $target_num, $k, $link_files[$j]);
	    runcommand("cp", "$source_dir/$source_fname $target_dir/$target_fname", 1);
	    runcommand("cp", "$source_dir/$source_fname $target_dir/$target_fname", 0);
	}
    }
}

