#!/usr/bin/perl

use strict;

my $debug=0;

if ($#ARGV!=4) {
    print "$0 <file> <target_dir> <year> <julian start day> <julian end day>\n";
    exit 0;
}

my @command;
my $fname = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];

my ($m1, $d1, $m2, $d2);

my @dirs;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
#    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    
    if (!-e "$target_dir$dirs[$i]") {
	mkdir("$target_dir$dirs[$i]");
    }
    if (!-e "$target_dir$dirs[$i]/run") {
	mkdir("$target_dir$dirs[$i]/run");
    }

    system("cp $fname $target_dir/$dirs[$i]/run/");
}
