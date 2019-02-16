#!/usr/bin/perl

use strict;
use POSIX qw(ceil floor);

if ($#ARGV!=5) {
    print "$0 <base dir> <target dir> <year> <julian start day> <julian end day> <db>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $db = $ARGV[5];
my ($finame, $foname);

my ($i, $j, $k);

my @link_files = ("2Dvac", "2Dvcc", "bvac", "bvcc", "bvmag", "davac", "davcc", "davmag", "svac", "svcc", "svmag");

for ($i=$sday; $i<=$eday; ++$i) {
    my $dir=sprintf("%d-%03d", $year, $i);
    
    for ($j=0; $j<=$#link_files; ++$j) {
	$foname = sprintf("%d-%03d-%s_%s.61", $year, $i, $db, $link_files[$j]);
	runcommand("ln", "-sf $base_dir/$dir/process/1_$link_files[$j].61 $target_dir/$foname", 0);
    }
}

########################################################

sub runcommand {
    my ( $e, $a, $debug ) = @_;

    if ($debug) { # in case $debug is undef
    } else {
       $debug = 0;
    }

    # invoke the real system call here
    my $err = 0;
    if ($debug) {
        print "$e $a\n";
    } else {
        $err = system("$e $a");
    }

    return $err;
}
