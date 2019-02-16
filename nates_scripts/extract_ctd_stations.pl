#!/usr/bin/perl

use strict;
#use CORIE;

my $debug=0;

if ($#ARGV!=1) {
    print "$0 <target_dir> <fout>\n";
    exit 0;
}

my $code = "/home/users/yinglong/SELFE/read_output7b_group_z2";

my @command;
my $source_dir = $ARGV[0];
my $fout = $ARGV[1];

#my @stas = ("van", "sth", "lon", "wau", "ska", "ast", "ham");
#my @stas = ("bon", "wil", "mor", "van", "cls", "sth", "lon", "bvr", "wau", "ska", "ast", "ham");
my $i;
my $j;

my @outstrings = ();

for ($i=0; $i<=$#stas; ++$i) {
    unlink("extract.out");

    open (APROG, "| $code") or die "couldn't do the fancy thing\n";
    print APROG "elev.61\n";
    print APROG "9\n";
    print APROG "$x[$i] $y[$i]\n";
    print APROG "0\n";
    close (APROG);

    sleep(10);
    
    open (IN, "extract.out") or die "could not open extract.out\n";
    my @input = <IN>;
    for ($j=0; $j<=$#input; ++$j) {
	$input[$j] =~ /(\S+)\s+(\S+)/;
	if ($i==0) {
	    $outstrings[$j] = $1." ".$2;
	} else {
	    $outstrings[$j] = $outstrings[$j]." ".$2;
	}
    }
    close (IN);

    print "finished $stas[$i]\n";
}

open (OUT, ">$fout") or die "could not open $fout for writing\n";
for ($i=0; $i<=$#outstrings; ++$i) {
    print OUT $outstrings[$i]."\n";
}
close(OUT);
