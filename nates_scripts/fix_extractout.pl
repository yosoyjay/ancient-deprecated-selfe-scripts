#! /usr/bin/perl -w

##########
########## fixes a (hopefully) one time problem where a run had to be restarted but the internal time was not fixed
##########

use strict;

if (@ARGV != 2) { 
    print "$0 <fin> <fout>\n";
    exit(1);
}
my $fin=$ARGV[0]; 
my $fout=$ARGV[1];

open(IN, $fin) or die "could not open $fin for reading\n";
my @input = <IN>;
close(IN);

open(OUT, ">$fout") or die "could not open $fout for writing\n";

my ($i, $j);
my $interval = 900/86400;

for ($i=0; $i<=$#input; ++$i) {
    my @vals = split(/\s+/, $input[$i]);
    print OUT "".(($i+1)*$interval)."";
    for ($j=1; $j<=$#vals; ++$j) {
	print OUT " ".$vals[$j];
    }
    print OUT "\n";
}

close(OUT);



