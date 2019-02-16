#!/usr/bin/perl

# assumes weekly model directories have been set up already
# this script assumes cron, or run manually 
# hard-coded vals rather than command-line args to facilitate easier cron | manual runs

use strict;

#    print "$0 <base_dir> <run name> <year> <start day> <end day> <nprocs> <external_steps_per_day>\n";

my @command;
my $home_dir = "/home/hyde/ccalmr46fvcom/credb11/credb11_6/run6/";
my $runname = "credb11";
my $year = "2006";
my $sday = "1";
my $eday = "48";
my $nprocs = "16";
my $secpday = "57600";
my $qsub = "qsub -pe orte $nprocs run.qsub";

my @dirs;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1] = sprintf("%d-%03d", $year, $i);
}


my $done = 0;
for ($i=0; $i<=$#dirs && $done==0; ++$i) {
    my $restart = sprintf("re_%06d.dat", $i+1);
    my $restart_wd = sprintf("re_%08d_wd", ($i+1)*$secpday);

    if (! -e "$home_dir$dirs[$i]/$restart") {
        if ($i>0) {
	    @command=("ln", "-sf", "$home_dir$dirs[$i]/$restart", "$home_dir$dirs[$i+1]/$runname"."_restart.dat");
	    print "@command\n";
	    system(@command);
	    @command=("ln", "-sf", "$home_dir$dirs[$i]/$restart_wd", "$home_dir$dirs[$i+1]/$runname"."_restart_wd.dat");
	    print "@command\n";
	    system(@command);
	}

	chdir("$home_dir$dirs[$i]");

	@command=("$qsub");
	print "\n$qsub from $home_dir$dirs[$i]\n";
	system(@command);
    }
}

if (!$done) {
    print "Did not run anything!\n";
}
