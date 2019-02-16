#!/usr/bin/perl

# assumes weekly model directories have been set up already
# also assumes it is starting from the beginning.  If not, manually copy 
# [casename]_restart.dat files from previous day

use strict;

my $debug=0;

if ($#ARGV!=6) {
    print "$0 <base_dir> <run name> <year> <start day> <end day> <nprocs> <external_steps_per_day>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0]."/";
my $runname = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $nprocs = $ARGV[5];
my $secpday = $ARGV[6];
my $qsub = "qsub -pe orte $nprocs run.qsub";

my @dirs;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1] = sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    chdir("$home_dir$dirs[$i]");

    print "doing $dirs[$i]\n";

    @command=("$qsub");
    print "\n$qsub from $home_dir$dirs[$i]\n";
    system(@command);

    print "stuff 1\n";

    my $restart = sprintf("re_%06d.dat", $i+1);
    my $restart_wd = sprintf("re_%08d_wd", ($i+1)*$secpday);
    while (! -e "$home_dir$dirs[$i]/$restart") {
        print "going to sleep...\n";
	sleep(1800);
    }
    sleep(60); # make sure the executable is done writing
    
    chdir("$home_dir"); 
    if ($i<$#dirs) { 
	print "inside1\n";

	if (-e "$home_dir$dirs[$i+1]/$runname"."_restart.dat") {
	    unlink("$home_dir$dirs[$i+1]/$runname"."_restart.dat");
	}
	if (-e "$home_dir$dirs[$i+1]/$runname"."_restart_wd.dat") {
	    unlink("$home_dir$dirs[$i+1]/$runname"."_restart_wd.dat");
	}

	@command=("ln", "-sf", "$home_dir$dirs[$i]/$restart", "$home_dir$dirs[$i+1]/$runname"."_restart.dat");
	print "@command\n";
	system(@command);
	@command=("ln", "-sf", "$home_dir$dirs[$i]/$restart_wd", "$home_dir$dirs[$i+1]/$runname"."_restart_wd.dat");
	print "@command\n";
	system(@command);
    }

    print "did $dirs[$i]\n";
}
