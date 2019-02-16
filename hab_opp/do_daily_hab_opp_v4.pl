#!/usr/bin/perl

### currently set to use neg_hab_opp_by_node.c, which produces file "regional_nodal_hab_opp_percent.dat"
###     neg_hab_opp_by_node requires a max depth argument in addition to regular hab_opp args

use strict;
use CORIE;

my $debug=0;

if ($#ARGV!=6) {
    print "$0 <script_dir> <filter> <target_dir> <regions file> 
<year> <julian start day> <julian end day>\n
	habopp_code options:
		0 - salmon
		1 - chinook_s
		2 - chinook_a
		3 - coho_s
		4 - coho_a
		5 - sockeye_s
		6 - sockeye_a
		7 - steelhead_s
		8 - steelhead_a
		9 - spawning
		10 - burrowing
		11 - metamorphose\n";
    exit 0;
}

my @command;
my $source_dir = $ARGV[0];
my $filter = $ARGV[1];
my $target_dir = $ARGV[2];
my $region_file = $ARGV[3];
my $year = $ARGV[4];
my $sday = $ARGV[5];
my $eday = $ARGV[6];

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

my $blah = $filter eq 'sdddddn';
print "$filter $blah \n";
for ($i=0; $i<=$#dirs; ++$i) {
    if (!-e "$target_dir/$dirs[$i]/process") {
	mkdir("$target_dir/$dirs[$i]/process");
    }
    
    if ($filter eq "salmon") {
        if (!-e "$target_dir/$dirs[$i]/process/salmon") {
        mkdir("$target_dir/$dirs[$i]/process/salmon");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/salmon $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/salmon $region_file 1" , 0);
    }
    elsif ($filter eq "chinook_s") {
        if (!-e "$target_dir/$dirs[$i]/process/chinook_s") {
        mkdir("$target_dir/$dirs[$i]/process/chinook_s");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/chinook_s $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/chinook_s $region_file 1" , 0);
    }
    elsif ($filter eq "chinook_a") {
        if (!-e "$target_dir/$dirs[$i]/process/chinook_a") {
        mkdir("$target_dir/$dirs[$i]/process/chinook_a");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/chinook_a $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/chinook_a $region_file 1" , 0);
    }
    elsif ($filter eq "coho_s") {
        if (!-e "$target_dir/$dirs[$i]/process/coho_s") {
        mkdir("$target_dir/$dirs[$i]/process/coho_s");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/coho_s $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/coho_s $region_file 1" , 0);
    }
    elsif ($filter eq "coho_a") {
        if (!-e "$target_dir/$dirs[$i]/process/coho_a") {
        mkdir("$target_dir/$dirs[$i]/process/coho_a");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/coho_a $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/coho_a $region_file 1" , 0);
    }
    elsif ($filter eq "sockeye_s") {
        if (!-e "$target_dir/$dirs[$i]/process/sockeye_s") {
        mkdir("$target_dir/$dirs[$i]/process/sockeye_s");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/sockeye_s $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/sockeye_s $region_file 1" , 0);
    }
    elsif ($filter eq "sockeye_a") {
        if (!-e "$target_dir/$dirs[$i]/process/sockey_a") {
        mkdir("$target_dir/$dirs[$i]/process/sockeye_a");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/sockeye_a $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/sockeye_a $region_file 1" , 0);
    }
    elsif ($filter eq "steelhead_s") {
        if (!-e "$target_dir/$dirs[$i]/process/steelhead_s") {
        mkdir("$target_dir/$dirs[$i]/process/steelhead_s");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/steelhead_s $region_file 2" , 
1);
        runcommand("$source_dir/hab_opp","$$filter target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/steelhead_s $region_file 1" , 
0);
    }
    elsif ($filter eq "steelhead_a") {
        if (!-e "$target_dir/$dirs[$i]/process/steelhead_a") {
        mkdir("$target_dir/$dirs[$i]/process/steelhead_a");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/steelhead_a $region_file 2" , 
1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/steelhead_a $region_file 1" , 
0);
    }
    elsif ($filter eq "spawning") {
        if (!-e "$target_dir/$dirs[$i]/process/spawning") {
        mkdir("$target_dir/$dirs[$i]/process/spawning");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/spawning $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/spawning $region_file 1" , 0);
    }
    elsif ($filter eq "burrowing") {
        if (!-e "$target_dir/$dirs[$i]/process/burrowing") {
        mkdir("$target_dir/$dirs[$i]/process/burrowing");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/burrowing $region_file 2" , 1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/burrowing $region_file 1" , 0);
    }
    else {
        if (!-e "$target_dir/$dirs[$i]/process/metamorphose") {
        mkdir("$target_dir/$dirs[$i]/process/metamorphose");
        }
	runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/metamorphose $region_file 2" , 
1);
        runcommand("$source_dir/hab_opp","$filter $target_dir/$dirs[$i]/run/ $target_dir/$dirs[$i]/process/metamorphose $region_file 1" , 
0);    
    }    
}

