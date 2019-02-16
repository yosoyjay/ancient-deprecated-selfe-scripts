#!/usr/bin/perl -w
#
# copyhdf.pl for PROTO forecasts
#
# Author: hyde, Apr 2008, adapted from copyhdf.pl, pturner Apr, 2005
#
#

use strict;
use CORIE;
use DBI;

my $debug = 0;
my $err;

if ($#ARGV != 2) {
    print "usage: $0 <year> <julian start day> <julian end day>\n";
    exit 0;
}

my $yr = $ARGV[0];
my $sday = $ARGV[1];
my $eday = $ARGV[2];

my $y = sprintf("%02d", $yr-2000);

# HDF file creation
my $base_dir = "/home/workspace/ccalmr38/hindcasts/";

chdir "$base_dir"."hdf_files/";

my $hdf_dir = "/home/workspace/ccalmr/elcirc/inputs/wind_files/hdf/";
my $i;

for ($i=$sday; $i<=$eday; ++$i) {
    my @md = getmonthday($yr, $i);
    my $m = sprintf("%02d", $md[0]);
    my $d = sprintf("%02d", $md[1]);
    my $jd = sprintf("%03d", $i);

    my $month_dir = "$hdf_dir$y$m/";
    # define the forcing data files
    my $rad_file = "eta_rad.$y$m$d.local.hdf";
    my $air_file = "eta_air.$y$m$d.local.hdf";
    # create the hdf directories

    my $outdir = "$base_dir"."hdf_files/$yr-$jd/";

    if (!(-d $outdir)) {
	$err = runcommand("mkdir", "$outdir", $debug);
    }
    $err = runcommand("/usr/bin/scp", "-pq $month_dir$rad_file $outdir"."flux_file_1.001.hdf", $debug);
    if ($err) {
	print "Unable to copy HDF files for $month_dir$rad_file, please check\n";
    }
    $err = runcommand("/usr/bin/scp", "-pq $month_dir$air_file $outdir"."wind_file_1.001.hdf", $debug);
    if ($err) {
	print "Unable to copy HDF files for $month_dir$air_file, please check\n";
    }
}
## Get GFS
#$rad_file = "gfs_rad.$ymd.hdf";
#$air_file = "gfs_air.$ymd.hdf";
## copy these files to the run directory
#$err = runcommand("/usr/bin/scp", "-pq $month_dir$rad_file gfs_flux_file_1.001.hdf", $debug);
#if ($err) {
#    print "Unable to copy HDF files for $month_dir$rad_file, please check\n";
#}
#$err = runcommand("/usr/bin/scp", "-pq $month_dir$air_file gfs_wind_file_1.001.hdf", $debug);
#if ($err) {
#    print "Unable to copy HDF files for $month_dir$air_file, please check\n";
#}
