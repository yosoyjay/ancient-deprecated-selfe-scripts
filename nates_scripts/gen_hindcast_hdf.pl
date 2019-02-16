#!/usr/bin/perl -w
#
# copyhdf.pl for PROTO forecasts
#
# Author: hyde, Apr 2008
#
#

use strict;
use CORIE;
use DBI;

my $debug = 0;
my $err;

if ($#ARGV != 4) {
    print "usage: $0 <base_dir> <forecast_name> <year> <julian start day> <julian end day>\n";
    exit 0;
}

my $base_dir = $ARGV[0];
my $fname = $ARGV[1];
my $yr = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];

#my $y = sprintf("%02d", $yr-2000);

#chdir "$base_dir"."hdf_files/";

my $i;

for ($i=$sday; $i<=$eday; ++$i) {
    my $jd = sprintf("%03d", $i);
    my $outdir = $base_dir."/$fname/$yr-$jd/hdf/";

    if (!(-d $outdir)) {
	$err = runcommand("mkdir", "$outdir", $debug);
    }
    chdir "$outdir";
    $err = runcommand("/home/workspace/ccalmr/mazulauf/scripts/make_wind_links_eta.csh", "$yr $i ".($i+1), $debug);
    if ($err) {
	print "Unable to comply\n";
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