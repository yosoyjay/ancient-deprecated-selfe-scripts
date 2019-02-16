#!/usr/bin/perl

use strict;
use POSIX qw(ceil floor);

if ($#ARGV!=4) {
    print "$0 <base dir> <target dir> <year> <julian start day> <julian end day>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];

my ($i, $j, $k);

for ($i=$sday; $i<=$eday; ++$i) {
    my $dir=sprintf("%d-%03d", $year, $i);

#    runcommand("ln", "-s $base_dir/$dir/run/ $target_dir/$dir/", 1);
    runcommand("rm", "$target_dir/$dir/process", 0);
    runcommand("ln", "-s $base_dir/$dir/process $target_dir/$dir/process", 0);

#    runcommand("mkdir", "$target_dir/$dir/process", 0);
#    runcommand("mv", "$base_dir/$dir/process/* $target_dir/$dir/process/", 0);
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
