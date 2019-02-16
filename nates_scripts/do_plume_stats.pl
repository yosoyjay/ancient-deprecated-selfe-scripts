#!/usr/bin/perl

use strict;

if ($#ARGV!=4) {
    print "$0 <base dir> <year> <julian start day> <julian end day> <plume psu cutoff>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $year = $ARGV[1];
my $sday = $ARGV[2];
my $eday = $ARGV[3];
my $plume_psu = $ARGV[4];

my $ex = "/home/users/hyde/models/plume/compute_plumevol_SELFE";
my $plume_def = "/home/users/hyde/models/plume/db22.flg";

# ./compute_plumevol_SELFE db14.flg /home/workspace/ccalmr/hindcasts/2000-10-14/run/ /home/users/hyde/models/plume/data/ 1 $plume_psu

my ($i, $j, $k);

for ($i=$sday; $i<=$eday; ++$i) {
    my $dir=sprintf("%d-%03d", $year, $i);

    if (!-e "$base_dir/$dir/process") {
	mkdir("$base_dir/$dir/process");
    }

    runcommand("$ex $plume_def $base_dir/$dir/run/ $base_dir/$dir/process/ 1 $plume_psu");
}

#####################

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
