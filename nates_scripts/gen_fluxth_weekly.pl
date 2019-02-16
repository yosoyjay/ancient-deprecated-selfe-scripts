#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

# in_file is assumed to start at day 1 (sday) and end at or after the last day (eday)

if ($#ARGV!=4) {
    print "$0 <in_file> <out_dir> <ts> <outfile base name> <% change>\n";
    exit 0;
}

my @command;
my $in_file = $ARGV[0];
my $target_dir = $ARGV[1];
my $ts = $ARGV[2];
my $out_name = $ARGV[3];
my $factor = 1+$ARGV[4];

my @input;

my $i;
my $j;
my @vars;
my $curr_week;
my $last_week;
my $curr_time;

open(FLUX, $in_file) or die "could not open $in_file for reading\n";
my @flux = <FLUX>;
close(FLUX);

$last_week = -1;
for ($i=0; $i<=$#flux; ++$i) {
    print "$i of $#flux\n";

    $curr_week = int($i/7)+1;
    if ($last_week != $curr_week) {
	if ($i!=0) {
	    close(FLUXTH);
	}
	open(FLUXTH, ">$target_dir/$out_name"."_$curr_week.th") or die "could not open $target_dir/$out_name"."_$curr_week.th for writing\n";
        $curr_time = $ts;
	$last_week = $curr_week;
    }
    
    for ($j=$ts; $j<=86400; $j=$j+$ts) {
	print FLUXTH "$curr_time\t-".($flux[$i]*$factor)."\n";
        $curr_time = $curr_time + $ts;
    }
}

