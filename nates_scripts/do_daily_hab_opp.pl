#!/usr/bin/perl

### currently set to use neg_hab_opp_by_node.c, which produces file "regional_nodal_hab_opp_percent.dat"
###     neg_hab_opp_by_node requires a max depth argument in addition to regular hab_opp args

use strict;
use CORIE;

my $debug=0;

if ($#ARGV!=5) {
    print "$0 <script_dir> <target_dir> <regions file> <year> <julian start day> <julian end day>\n";
    exit 0;
}

my @command;
my $source_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $region_file = $ARGV[2];
my $year = $ARGV[3];
my $sday = $ARGV[4];
my $eday = $ARGV[5];

my $mdepth = 2;   ### hard coded value maybe should be in command line args...
#my $mdepth = 50;   ### hard coded value maybe should be in command line args...

# quick tasks to every directory.  Change tasks as needed

my @dirs;
my @input;

my $i;
my $j;
my @vars;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    if (!-e "$target_dir/$dirs[$i]/process") {
	mkdir("$target_dir/$dirs[$i]/process");
    }
    runcommand("$source_dir/neg_hab_opp_by_node_z2","$target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/ $region_file 1 $mdepth" , 1);
    runcommand("$source_dir/neg_hab_opp_by_node_z2","$target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/ $region_file 1 $mdepth" , 0);
#    runcommand("$source_dir/hab_opp_by_node","$target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/ $region_file 1 $mdepth" , 0);
}

