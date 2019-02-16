#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

# in_file is assumed to start at day 1 (sday) and end at or after the last day (eday)

if ($#ARGV!=1) {
    print "$0 <source dir> <target dir>\n";
    exit 0;
}

my $in_dir = $ARGV[0];
my $out_dir = $ARGV[1];
my $i;
my @vals;

my %clims = (
	     '14' => 'evmp',
             '16' => 'evme'
 	     );

my %prods = ( 
	      'MIN' => 'min',
              'MAX' => 'max',
	      'AVG' => 'avg',
	      'AMAX' => 'amax',
	      'AMIN' => 'amin',
	      'AAVG' => 'aavg',
	      'RANGE' => 'rng',
	      'ARANGE' => 'arng'
	     );
my %source = (
	      '14' => 'db14',
              '16' => 'db16'
 	      );
my %interval = (
	      '01' => 'jan',
	      '02' => 'feb',
	      '03' => 'mar',
	      '04' => 'apr',
	      '05' => 'may',
	      '06' => 'jun',
	      '07' => 'jul',
	      '08' => 'aug',
	      '09' => 'sep',
	      '10' => 'oct',
	      '11' => 'nov',
	      '12' => 'dec',
	      'summer' => 'sum', 
	      'spring' => 'spr',
	      'fall' => 'fal',
	      'winter' => 'win'
               );


opendir(DIR, $in_dir);
my $file;
while ( $file = readdir(DIR)) {
    next if $file eq "." or $file eq "..";
    @vals = split(/-/, $file);

    $vals[$#vals] =~ /(.*)\.gif/;

    my $new_file = $clims{$1}."_".$vals[0]."_".$prods{$vals[1]}."_".$source{$1}."_";

    my $yr1 = $vals[2];
    my $yr2 = $vals[3];
    my $yr = "all";
    if ($yr1 == $yr2) {
	$yr = $yr1;
    }
    my $int = "all";
    if ($#vals == 5) {
	$int = $interval{$vals[4]};
    }

    $new_file = $new_file.$int."_$yr.gif";

    print "$file -> $new_file\n";

    system("cp", "$in_dir/$file", "$out_dir/$new_file");
#    system("mv temp.tmp $file");
}
closedir(DIR);
