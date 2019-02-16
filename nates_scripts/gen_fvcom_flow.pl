#! /usr/bin/perl -w
#
# uses files from /home/workspace/ccalmr/data/river_data/verified/, such as columbia_discharge.2004 [corie_day discharge(m^3/s)]
#
# input file is:
#   river node numbers
#   level percentages
#  e.g. (for nodes 1-4, 11 vertical levels - or 10 vertical sides):
   #1 2 3 4     
   #0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1

use strict;
use CORIE;

my $debug=0;

if ($#ARGV !=5) {
    print "$0 <base_dir> <year> <julian start day> <julian end day> <model dt> <input args file> <output_file>\n";
    exit 1; 
}
my $base_dir = $ARGV[0];
my $year = $ARGV[1];
my $sday = $ARGV[2];
my $eday = $ARGV[3];
my $dt = $ARGV[4];
my $fin = $ARGV[5];
my $fout = $ARGV[6];

my $flow_dir = "/home/workspace/ccalmr/data/river_data/verified/";
my $flow_file = "columbia_discharge";
my $temp_file = "bonneville_temperature";

my @dirs;

my ($i, $j);

##### get the input data #####
open(DIS, "$flow_dir$flow_file."."$year") or die "could not open $flow_dir$flow_file."."$year\n";
my @discharge = <DIS>;
close(DIS);
open(TMP, "$flow_dir$temp_file."."$year") or die "could not open $flow_dir$temp_file."."$year\n";
my @temp = <TMP>;
close(TMP);

##### get fvcom args and start casename_riv.dat file ####
open(OUT,">$fout") or die "could not open $fout for writing\n";



for ($i=$sday; $i<=$eday; ++$i) {

    my ($m, $d) = getmonthday($year, $i);
    
    for ($j=$dt; $j<=86400; $j=$j+$dt) {
	
    }


    
}
