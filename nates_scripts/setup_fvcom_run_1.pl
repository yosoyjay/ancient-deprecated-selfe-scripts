#!/usr/bin/perl

use strict;
use CORIE;
use POSIX qw(ceil floor);

my $debug=0;

if ($#ARGV!=8) {
    print "usage: $0 <base_dir> <year> <start day> <end day> <copy dir> <casename> <time step> <isplit> <ramp-up period>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0];
my $year = $ARGV[1];
my $sday = $ARGV[2];
my $eday = $ARGV[3];
my $ref_dir = $ARGV[4];
my $casename = $ARGV[5];
my $tstep = $ARGV[6];
my $isplit = $ARGV[7];
my $ramp = $ARGV[8];

my @dirs;
my @days;
my ($m1, $m2, $d1, $d2);

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1] = sprintf("%d-%03d", $year, $i);
    $days[$#days+1] = $i;
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    
    if (!-e "$home_dir$dirs[$i]") {
	mkdir("$home_dir/$dirs[$i]");
    }

    chdir("$home_dir/$dirs[$i]");
    
    runcommand("cp", "$ref_dir/* $home_dir/$dirs[$i]/", 1);
    runcommand("cp", "$ref_dir/* $home_dir/$dirs[$i]/", 0);

    adjust_rundat("$casename"."_run.dat", $tstep, $isplit, $ramp, $days[$i]);
    adjust_runqsub("run.qsub", $home_dir, $dirs[$i], $casename);
}



##############################################

sub adjust_rundat() {
    my $runfile = $_[0];
    my $dt      = $_[1];
    my $isplit  = $_[2];
    my $ramp    = $_[3];
    my $day     = $_[4];

    my $jj;

    open(RUNFILE, "$runfile") or die "could not open $runfile for reading\n";
    my @run = <RUNFILE>;
    close(RUNFILE);

    open(RUNFILE, ">$runfile") or die "could not open $runfile for writing\n";

    for ($jj=0; $jj<=$#run; ++$jj) {
	if ($run[$jj] =~ /DTE[\s|=]/) {
	    print RUNFILE "DTE    =    $dt\n";
	}
	elsif ($run[$jj] =~ /ISPLIT[\s|=]/) {
	    print RUNFILE "ISPLIT =   $isplit\n";
	}
	elsif ($run[$jj] =~ /NSTEPS[\s|=]/) {
            print RUNFILE "NSTEPS = ".($day*86400/($isplit*$dt))."\n";
	}
	else {
	    print RUNFILE "$run[$jj]";
	}
    }
    close(RUNFILE);
}

sub adjust_runqsub() {
    my $runfile  = $_[0];
    my $basedir  = $_[1];
    my $dir      = $_[2];
    my $casename = $_[3];

    my $exec_line = "/usr/mpi/gcc/openmpi-1.2.5/bin/mpirun --mca btl_tcp_if_include eth0 -np \$NSLOTS $basedir/$dir/fvcom_1 $casename\n";
    my $jj;

    open(RUNFILE, "$runfile") or die "could not open $runfile for reading\n";
    my @run = <RUNFILE>;
    close(RUNFILE);

    open(RUNFILE, ">$runfile") or die "could not open $runfile for writing\n";

    for ($jj=0; $jj<=$#run; ++$jj) {
	chomp($run[$jj]);
	if ($run[$jj]=~/mpirun/) {
	    print RUNFILE "$exec_line\n";
	}
	else {
	    print RUNFILE "$run[$jj]\n";
	}
    }
    close(RUNFILE);
}

