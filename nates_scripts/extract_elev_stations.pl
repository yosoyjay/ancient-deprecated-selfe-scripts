#!/usr/bin/perl

use strict;
#use CORIE;

my $debug=0;

if ($#ARGV!=1) {
    print "$0 <target_dir> <fout>\n";
    exit 0;
}

my $code = "/home/users/yinglong/SELFE/read_output7b";

my @command;
my $source_dir = $ARGV[0];
my $fout = $ARGV[1];

# my @xy = ((4.38456e+005, 2.20653e+005), (4.31539e+005, 2.46803e+005), (4.19968e+005, 2.73806e+005), (3.84362e+005, 2.82938e+005), (3.82265e+005, 2.92944e+005), (3.57432e+005, 2.87557e+005), (3.43733e+005, 2.87572e+005));
#my @x = (4.38456e+005, 4.31618e+005, 4.19968e+005, 3.84362e+005, 3.82265e+005, 3.57432e+005, 3.43733e+005);
#my @y = (2.20653e+005, 2.46634e+005, 2.73806e+005, 2.82938e+005, 2.92944e+005, 2.87557e+005, 2.87572e+005);
#my @x = (494759, 442098, 440341, 4.38456e+005, 433282, 4.31618e+005, 4.19968e+005, 402530, 3.84362e+005, 3.82265e+005, 3.57432e+005, 3.43733e+005);
#my @y = (218816, 188191, 208059, 2.20653e+005, 221660, 2.46634e+005, 2.73806e+005, 282530, 2.82938e+005, 2.92944e+005, 2.87557e+005, 2.87572e+005);

### the following reflect pre-shift (i.e. off by [80,20]) station points
# my @x = (494759, 444326, 440341, 4.38456e+005, 432200, 4.31618e+005, 4.19968e+005, 402530, 3.84362e+005, 3.82265e+005, 3.57432e+005, 3.43733e+005);
# my @y = (218816, 190078, 208059, 2.20653e+005, 221660, 2.46634e+005, 2.73806e+005, 282832, 2.82938e+005, 2.92944e+005, 2.87557e+005, 2.87572e+005);

### corrected points:
my @x = (494759-80, 444326-80, 440341-80, 4.38456e+005-80, 432200-80, 4.31618e+005-80, 4.19968e+005-80, 402530-80, 3.84362e+005-80, 3.82265e+005, 3.57432e+005-80, 3.43827e+005);
my @y = (218916-20, 190078-20, 208059-20, 2.20653e+005-20, 221660-20, 2.46634e+005-20, 2.73806e+005-20, 282832-20, 2.82938e+005-20, 2.92944e+005, 2.87557e+005-20, 2.87572e+005);

#my @stas = ("van", "sth", "lon", "wau", "ska", "ast", "ham");
my @stas = ("bon", "wil", "mor", "van", "cls", "sth", "lon", "bvr", "wau", "ska", "ast", "ham");
my $i;
my $j;

my @outstrings = ();

for ($i=0; $i<=$#stas; ++$i) {
    unlink("extract.out");

    open (APROG, "| $code") or die "couldn't do the fancy thing\n";
    print APROG "elev.61\n";
    print APROG "11\n";
    print APROG "$x[$i] $y[$i]\n";
    print APROG "3\n";
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
