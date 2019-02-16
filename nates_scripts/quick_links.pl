#!/usr/bin/perl

use strict;
use POSIX qw(ceil floor);

if ($#ARGV!=5) {
    print "$0 <base dir> <target dir> <year> <julian start day> <julian end day> <source db>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $db = $ARGV[5];

my ($i, $j, $k);

my @link_files = ("elev.61", "hvel.64", "temp.63", "salt.63", "zcor.63", "hgrid.gr3");

for ($i=$sday; $i<=$eday; ++$i) {
    my $dir=sprintf("%d-%03d", $year, $i);

    my $w = floor(($i-1)/7)+1;
    my $local_day = ($i-1) % 7 + 1;

    if (!-e "$target_dir/$dir") {
	mkdir("$target_dir/$dir");
    }
    if (!-e "$target_dir/$dir/process") {
	mkdir("$target_dir/$dir/process");
    }

#    runcommand("ln", "-s $base_dir/$dir/run/ $target_dir/$dir/", 1);
    runcommand("ln", "-s $base_dir/$dir/run/ $target_dir/$dir/", 0);
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
