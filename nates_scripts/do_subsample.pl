#!/usr/bin/perl

use CORIE;

if ($#ARGV < 3) {
    print "usage: $0 <fin> <fout> <subsample factor> <subsample index>\n";
    exit 1;
}

my $fin = $ARGV[0];
my $fout = $ARGV[1];
my $subsample_factor = $ARGV[2];
my $index = $ARGV[3];

open(IN, "$fin") or die "could not open $fin for reading\n";
open(OUT, ">$fout") or die "could not open $fout for writing\n";

#my @input = <IN>;
#close(IN);
#my $i;
#for($i=0; $i<=$#input; ++$i) {
#    if (($i+1) % $subsample_factor == $index) {
#        print OUT $input[$i];
#    }
#}
#close(OUT);

my $i=0;
while (!eof (IN)) {
    $line = <IN>;
    @vals = split(/\s+/, $line);
    if (($i+1) % $subsample_factor == $index) { 
        print OUT $line;
    }
    ++$i;
}
close(OUT);
close(IN);

