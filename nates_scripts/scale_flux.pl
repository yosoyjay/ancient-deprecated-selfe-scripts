#!/usr/bin/perl

use strict;

if (@ARGV != 5) {
    die "usage: $0 <infile> <outfile> <adjusted column> <factor> <factor_min>\n";
}
# 

open(INFILE, "$ARGV[0]") or die "could not open $ARGV[0] for reading\n";
my @input = <INFILE>;
close(INFILE);
open(OUTFILE, ">$ARGV[1]");

my $adj_col = $ARGV[2];
my $factor = $ARGV[3];
my $factor_min = $ARGV[4];
my @vals;

my ($i, $j);

for ($i=0; $i<=$#input; ++$i) { 
    @vals = split(/\s+/, $input[$i]);
#    print "$#vals    @vals\n";
    if ($vals[$adj_col] < -$factor_min) {   #flux values are negative
        $vals[$adj_col] = -$factor_min+($vals[$adj_col]+$factor_min)*$factor;
    }
    for ($j=0; $j<=$#vals; ++$j) {
	print OUTFILE "$vals[$j]\t";
    }
    print OUTFILE "\n";
	
}
close(OUTFILE);
