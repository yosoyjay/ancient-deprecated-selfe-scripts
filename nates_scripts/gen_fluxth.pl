#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

# in_file is assumed to start at day 1 (sday) and end at or after the last day (eday)

if ($#ARGV!=7) {
    print "$0 <in_file> <out_dir> <year> <julian start day> <julian end day> <ts> <final time> <outfile name>\n";
    exit 0;
}

my @command;
my $in_file = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $ts = $ARGV[5];
my $last_time = $ARGV[6];
my $out_name = $ARGV[7];

my @dirs;
my @input;

my $i;
my $j;
my @vars;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

open(FLUX, $in_file) or die "could not open $in_file for reading\n";
my @flux = <FLUX>;
close(FLUX);

if ($#flux < $#dirs) {
    print "more days than flux entries!\n";
    exit(1);
}

for ($i=0; $i<=$#dirs; ++$i) {
    open(FLUXTH, ">$target_dir/$dirs[$i]/run/$out_name") or die "could not open $target_dir/$dirs[$i]/run/$out_name for writing\n";
    
    for ($j=$ts; $j<=$last_time; $j=$j+$ts) {
	chomp($flux[$i]);
	print FLUXTH "$j\t$flux[$i]\n";
    }

    close(FLUXTH);
}

