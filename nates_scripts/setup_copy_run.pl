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

    runcommand("cp","$source_dir/$dirs[$i]/run/flux.th $target_dir/$dirs[$i]/run/", 0);
    runcommand("cp","$source_dir/$dirs[$i]/run/temp.th $target_dir/$dirs[$i]/run/", 0);
#
#    runcommand("ln","-sf $source_dir/$dirs[$i]/run/flux.th $target_dir/$dirs[$i]/run/", 0);
#    runcommand("ln","-sf $source_dir/$dirs[$i]/run/temp.th $target_dir/$dirs[$i]/run/", 0);
#
#    runcommand("ln","-sf $source_dir/$dirs[$i]/run/elev3D.th $target_dir/$dirs[$i]/run/", 0);
#    runcommand("ln","-sf $source_dir/$dirs[$i]/run/temp3D.th $target_dir/$dirs[$i]/run/", 0);
#    runcommand("ln","-sf $source_dir/$dirs[$i]/run/salt3D.th $target_dir/$dirs[$i]/run/", 0);
#    runcommand("ln","-sf $source_dir/$dirs[$i]/run/uv3D.th $target_dir/$dirs[$i]/run/", 0);

#    runcommand("cp","$ref_dir/copy/hgrid.gr3 $target_dir/$dirs[$i]/run/", 0);
#    runcommand("cp","$ref_dir/copy/* $target_dir/$dirs[$i]/run/", 0);

#    @command=("cp","$source_dir/$dirs[$i]/run/*.th", "$target_dir/$dirs[$i]/run");
#    system(@command);
#    runcommand("cp","$source_dir/$dirs[$i]/run/*.in $target_dir/$dirs[$i]/run/", 0);
#    runcommand("ln","-s $source_dir/$dirs[$i]/run/hdf $target_dir/$dirs[$i]/run/hdf", 0);

#    runcommand("perl", "/home/workspace/ccalmr38/hindcasts/scripts/change_flux_dt.pl $source_dir/$dirs[$i]/run/flux.th $target_dir/$dirs[$i]/run/flux.th $new_ts 86400", 0);
#    runcommand("perl", "/home/workspace/ccalmr38/hindcasts/scripts/change_flux_dt.pl $source_dir/$dirs[$i]/run/temp.th $target_dir/$dirs[$i]/run/temp.th $new_ts 86400", 0);

    #seems to be some files getting zipped...
#    runcommand("rm","$temp_dir/*.in", 0);
#    runcommand("cp","$source_dir/$dirs[$i]/run/*_nu.in.gz $temp_dir", 0);
#    runcommand("gunzip", "$temp_dir/*.gz", 0);
#    runcommand("cp","$temp_dir/*.th $target_dir/$dirs[$i]/run", 0);
#    runcommand("cp","$temp_dir/*_nu.in $target_dir/$dirs[$i]/run", 0);
#    runcommand("rm","$temp_dir/*", 0);
#
#    ($m1, $d1) = getmonthday($year, $sday+$i);
#    ($m2, $d2) = getmonthday($year, $sday+$i+1);
#
#    chdir("$target_dir/$dirs[$i]/run");
#
#    runcommand("/home/workspace/ccalmr/mazulauf/amb10xx/netcdf/cvs_stuff/forecasts/bin/atmos_nc/scripts/make_sflux_links.csh","1 $year $m1 $d1 $year $m2 $d2", 0);
#
#    runcommand("cp","$ref_dir/* $target_dir/$dirs[$i]/run", 0);
}
