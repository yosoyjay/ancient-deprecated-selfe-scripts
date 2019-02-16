#!/usr/bin/perl

use strict;

if (@ARGV != 4) {
    die "usage: $0 <infile> <outfile> <interval> <last time>\n";
}

open(INFILE, "$ARGV[0]");
open(OUTFILE, ">$ARGV[1]");
my @input = <INFILE>;

my $interval = $ARGV[2];
my $last_t = $ARGV[3];

my ($n, $i, $ii, $j);
my $in_ts=0;
my $last_in_ts;
my @last_in_data;
my @vals;
my @in_data;
my $in_interval;

@vals = split(/\s+/, $input[0]);
$in_interval = $vals[0];

#print "in_interval: $in_interval\n";

if ($interval <= $in_interval) {
    @input = (($input[0]), @input);
}

for ($n = 1; $n<=$#vals; ++$n) {
    $in_data[$n-1] = $vals[$n];
}

print "input interval: $in_interval : input data: ";
for ($n = 0; $n<=$#in_data; ++$n) {
    print "$in_data[$n]\t";
}
print "\n";

for ($i=$interval; $i<=$last_t; $i=$i+$interval) {

    print "$i $in_ts\n";

    if ($in_ts == 0 || $in_ts < $i) {
	while ($in_ts < $i && $in_ts < $last_t) {
#	    print "$in_ts < $i ??? $in_ts < $last_t\n";
#	    sleep(10);

  	    $last_in_ts = $in_ts;
  	    @last_in_data = @in_data;  	  
 
  	    ++$ii;
  	    $in_ts = $in_ts + $in_interval;
  	    @vals = split(/\s+/, $input[$ii]);
  	    for ($n = 1; $n<=$#vals; ++$n) {
  	        $in_data[$n-1] = $vals[$n];
  	    }
  	}
    }
    print OUTFILE "$i";


#    print "$in_data[$0] * ( $i - $last_in_ts ) / $in_interval + $last_in_data[0] * ( $in_ts - $i ) / $in_interval ) \n";

    for ($n = 0; $n<=$#in_data; ++$n) {
	print OUTFILE "\t".($in_data[$n] * ($i - $last_in_ts)/$in_interval + $last_in_data[$n] * ($in_ts - $i)/$in_interval);
    }
    print OUTFILE "\n";
	
}

