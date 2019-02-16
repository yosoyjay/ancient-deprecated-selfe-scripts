#!/usr/bin/perl -w
#
# do_intrusionlength.pl
#
# From original do_intrusionlength.pl modified by pturner
#
# Modified to work with non standard products, output plotted in Matlab 
#
# Example:
# do_intrusionlength.pl path/to/station.sta 4 path/to/run/dir 
#
#
use strict;
use CORIE;

my $debug = 0;
my $err   = 0;

if ( scalar @ARGV != 4 ) {
    print "Error: Number of command line arguments is incorrect, must be 3\n";
    print "Usage: $0 station_file ndays run_dir psu\n";
    exit(1);
}

my $sta      = $ARGV[0];
my $ndays    = $ARGV[1];
my $path     = $ARGV[2];
my $psu      = $ARGV[3];

if ( !( -e $sta ) ) {
    print STDERR "Error: station_file does not exist: $sta\n";
    exit(2)
}

for (my $i = $psu; $i <= $psu; $i++) {
  my @out_file = ("");
  my $cutoff      = $i;
  my $target_file = "intrusion_length$cutoff.dat";
  for (my $j = 5; $j <= $ndays; $j++) {
    my $salt = "$path/outputs/$j"."_salt.63";
    $out_file[$j] = "$j"."_intrusion_length$cutoff";
    $err = runcommand("/home/workspace/users/lopezj/scripts/selfe_post/intrusion_length/compute_intrusionlength", "$sta $salt $out_file[$j] $cutoff ", $debug);
  }
  $err = runcommand("/bin/rm", "-f $target_file ", $debug);
  $err = runcommand("/bin/cat", "@out_file > $target_file ", $debug);
}
