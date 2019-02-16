#!/usr/bin/perl

use strict;

my $db="db16";
my $base_dir = "/home/workspace/ccalmr38/hindcasts/columbia/";
my $shapes_dir = "/home/workspace/ccalmr38/hindcasts/columbia/shapes/";
my $out_dir = "/home/workspace/ccalmr32/hyde/columbia/db16/processed/";
# my $executable = "/home/users/hyde/ourestuaries/pythonlib/cmop/do_daily_process_transects.py";
# my $executable = "/home/users/hyde/ourestuaries/pythonlib/cmop/do_daily_process_cs.py";
my $executable = "/home/users/hyde/ourestuaries/pythonlib/cmop/do_daily_process_pts.py";

#nohup python /home/users/hyde/ourestuaries/pythonlib/cmop/do_daily_process_transects.py db16 /home/workspace/ccalmr38/hindcasts/columbia/ /home/workspace/ccalmr38/hindcasts/columbia/shapes/ 1999 101 100 /home/workspace/ccalmr32/hyde/columbia/db16/processed/

my $debug=0;

my %years;

#@{ $years{"1999"} } = (0, 50, 100, 150, 200, 250, 300, 365);
#@{ $years{"2000"} } = (0, 50, 100, 150, 200, 250, 300, 365);
#@{ $years{"2001"} } = (0, 50, 100, 150, 200, 250, 300, 365);
#@{ $years{"2002"} } = (0, 50, 100, 150, 200, 250, 300, 365);
#@{ $years{"2003"} } = (0, 50, 100, 150, 200, 250, 300, 365);
#@{ $years{"2004"} } = (0, 50, 100, 150, 200, 250, 300, 365);
#@{ $years{"2005"} } = (0, 50, 100, 150, 200, 250, 300, 364);
#@{ $years{"2006"} } = (0, 50, 100, 150, 200, 250, 300, 364);

@{ $years{"2005"} } = (100, 150, 200);
#@{ $years{"2006"} } = (300, 364);

my @keyst = keys (%years);
my @keys = sort(@keyst);

my ($i,$j,$k);
my $key;
my $day;

foreach $key (@keys) {
    print "$key\n";

    my @days = @{ $years{$key} };
    for ($j=0; $j<$#days; ++$j) {
	if ($j%2 == 0) {
	    runcommand("nohup nice", "python $executable $db $base_dir $shapes_dir $key ".($days[$j]+1)." ".($days[$j+1]-$days[$j])." $out_dir &", 0);
	}
	else {
	    runcommand("nohup nice", "python $executable $db $base_dir $shapes_dir $key ".($days[$j]+1)." ".($days[$j+1]-$days[$j])." $out_dir", 0);
	}
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