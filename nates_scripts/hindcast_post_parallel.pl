#!/usr/bin/perl

############ Note, this script is not complete ############

use CORIE;
use strict;

if ($#ARGV != 5) {
    print "usage: perl $0 <outputs_dir> <hindcast_base_dir> <n processors> <day 1 year> <day 1 day of year>\n";
    exit 1;
}

my $outputs = $ARGV[0];
my $hindcast_dir = $ARGV[1];
my $nprocs = $ARGV[2];
my $day1_yr = $ARGV[3];
my $day1_day = $ARGV[4];

my @files = ('airt.61', 'elev.61', 'hvel.64', 'salt.63', 'temp.63', 'wind.62', 'zcor.63');
my $i;


opendir(DIR, $outputs);
my $file;

my @to_do = ();

while ( $file = readdir(DIR)) {
    next if $file eq "." or $file eq "..";
    if ($file ~= /1_0000_$files[$i]/) {
	
    }
}
closedir(DIR);
