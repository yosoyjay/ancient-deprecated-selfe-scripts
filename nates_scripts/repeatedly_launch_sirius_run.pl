#!/usr/bin/perl

use strict;

#### code will try to launch a process using 16 cpus.  Resource issues tend to crush these runs, but repeated attempts 
#### eventually work.  This should do the legwork

#### written with time pressure.  Very stupid but hopefully workable algorithm and assumptions
####  -- not finished!!!!!!!!!!!!



my $ntries = 10;
my $sleep_secs = 3;


my $command = "qsub -pe orte 16 run.qsub";
my $check_command = "qstat > tmp.qstat";
my $i;

for ($i=1; $i<=$ntries; ++$i) {
    system($check_command);

    my @fcheck;
    open(IN, "tmp.qstat");
    @fcheck = <IN>;

    print "length fcheck = ".scalar(@fcheck)."\n";
}
