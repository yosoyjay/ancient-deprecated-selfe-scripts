#!/usr/bin/perl

use strict;
use POSIX qw(ceil floor);

if ($#ARGV!=4) {
    print "$0 <base dir> <output dir> <year> <julian start day> <julian end day>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $out_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];

my $fout = "hab_opp_$year"."_$sday"."_$eday";

my ($i, $j, $k);
my @vals;
my $variables;
my @values;

for ($i=$sday; $i<=$eday; ++$i) {
    my $dir=sprintf("%d-%03d", $year, $i);
    print "$dir\n";

    open(INPUT, "$base_dir"."/$dir/process/regional_habitat.ascii") or die "could not open $base_dir"."/$dir/process/regional_habitat.ascii\n";
    my @input = <INPUT>;
    close(INPUT);

    if ($i==$sday) {
	chomp ($input[1]);
	$variables = $input[1];
    }
    chomp($input[2]);
    $values[$i-$sday] = $input[2]; #split(/\s+/, $input[2]);
}

open(DAT, ">$out_dir$fout".".dat") or die "could not open $fout".".dat for writing\n";
open(CSV, ">$out_dir$fout".".csv") or die "could not open $fout".".csv for writing\n";

print(DAT $variables."\n");

for ($i=0; $i<=$#values; ++$i) {
    print (DAT "$values[$i]\n");
    @vals = split(/\s+/, $values[$i]);
    for ($j=1; $j<=$#vals; ++$j) {
	print(CSV "$vals[$j]");
	if ($j<$#vals) {
	    print(CSV ", ");
	}
	else {
	    print(CSV "\n");
	}
    }
}

close(DAT);
close(CSV);
