#!/usr/bin/perl

use strict;

my $debug=0;

if ($#ARGV!=6) {
    print "$0 <base_dir> <run name> <year> <julian start day> <julian end day> <executable> <hotstart prefix>\n";
}

my @command;
my $home_dir = $ARGV[0]."/".$ARGV[1]."/";
my $year = $ARGV[2];

my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $ex = $ARGV[5];
my $hs_pre = $ARGV[6];

my @dirs;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {

    chdir("$home_dir$dirs[$i]/run");

    @command=("./$ex");
    print "\n./$ex from $home_dir$dirs[$i]/run\n";
    system(@command);
    
    chdir("$home_dir");
   
    if ($i<$#dirs) {
	if (-e "$home_dir$dirs[$i+1]/run/hotstart.in") {
	    unlink("$home_dir$dirs[$i+1]/run/hotstart.in");
	}

	@command=("ln", "-sf", "$home_dir$dirs[$i]"."/run/$hs_pre"."_hotstart", "$home_dir$dirs[$i+1]/run/hotstart.in");
	print "@command\n";
	system(@command);
    }
}
