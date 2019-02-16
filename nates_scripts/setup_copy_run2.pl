#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

if ($#ARGV!=7) {
    print "$0 <reference_dir> <source_dir> <target_dir> <year> <julian start day> <julian end day> <temp_dir> <new_ts>\n";
    exit 0;
}

my @command;
my $ref_dir = $ARGV[0];
my $source_dir = $ARGV[1];
my $target_dir = $ARGV[2];
my $year = $ARGV[3];
my $sday = $ARGV[4];
my $eday = $ARGV[5];
my $temp_dir = $ARGV[6];
my $new_ts = $ARGV[7];

my ($m1, $d1, $m2, $d2);

my @dirs;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    
    if (!-e "$target_dir$dirs[$i]") {
	mkdir("$target_dir$dirs[$i]");
    }
    if (!-e "$target_dir$dirs[$i]/run") {
	mkdir("$target_dir$dirs[$i]/run");
    }

    runcommand("rm -r", "$target_dir$dirs[$i]/run/hdf", 0);
    mkdir("$target_dir$dirs[$i]/run/hdf");


    chdir("$target_dir$dirs[$i]/run/hdf");

    runcommand("/home/workspace/ccalmr/mazulauf/scripts/make_wind_links_gfs_global.csh","$year ".($sday+$i)." ".($sday+$i+2), 0);
    runcommand("cp","$ref_dir/zelfe1_3m_sflux7d_zhang_1 $target_dir/$dirs[$i]/run", 0);
}
