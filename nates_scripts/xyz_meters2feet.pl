#!/usr/bin/perl

use CORIE;

if ($#ARGV < 2) {
    print "usage: $0 <fin> <fout> <0|1: 0 to go m to ft, 1 to go ft to m>\n";
    exit 1;
}

my $fin = $ARGV[0];
my $fout = $ARGV[1];
my $factor = 1/0.3048;
if ($ARGV[2] == 1) {
    $factor = 0.3048;
}

open(IN, "$fin") or die "could not open $fin for reading\n";
open(OUT, ">$fout") or die "could not open $fout for writing\n";

my $i=0;
while (!eof (IN)) {
    $line = <IN>;
    @vals = split(/\s+/, $line);
    print OUT ($vals[0]*$factor)." ".($vals[1]*$factor)." $vals[2]\n";
    ++$i;
}
close(OUT);
close(IN);

