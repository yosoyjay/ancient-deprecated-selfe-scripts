#!/usr/bin/perl

use CORIE;

if ($#ARGV < 5) {
    print "usage: $0 <fin> <fout> <upper left x> <ul y> <lr x> <lr y> <optional: max height>\n";
    exit 1;
}

my $fin = $ARGV[0];
my $fout = $ARGV[1];
my $ulx = $ARGV[2];
my $uly = $ARGV[3];
my $lrx = $ARGV[4];
my $lry = $ARGV[5];

my $max_ht = 5; #99999999999;
if ($#ARGV==6) {
    my $max_ht = $ARGV[6];
}

open(IN, "$fin") or die "could not open $fin for reading\n";
open(OUT, ">$fout") or die "could not open $fout for writing\n";

my $i;
while (!eof (IN)) {
    $line = <IN>;
    @vals = split(/\s+/, $line);
#    if ( $vals[0] >= $ulx && $vals[0] <= $lrx && $vals[1] <= $uly && $vals[1] >= $lry ) {
    if ( $vals[0] >= $ulx && $vals[0] <= $lrx && $vals[1] <= $uly && $vals[1] >= $lry && $vals[2] < $max_ht) {
        print OUT $line;
    }
}
close(OUT);
close(IN);
