#!/usr/bin/perl

use strict;
use CORIE;
use POSIX qw(ceil floor);

my $debug=0;

if ($#ARGV!=6) {
    print "usage: $0 <base_dir> <run name> <year> <julian start day> <julian end day> <ref dir> <time step>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0]."/".$ARGV[1]."/";
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $ref_dir = $ARGV[5];
my $tstep = $ARGV[6];

my @dirs;
my ($m1, $m2, $d1, $d2);

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    
    if (!-e "$home_dir$dirs[$i]") {
	mkdir("$home_dir$dirs[$i]");
    }
    if (!-e "$home_dir$dirs[$i]/run") {
	mkdir("$home_dir$dirs[$i]/run");
    }

    ($m1, $d1) = getmonthday($year, $sday+$i);
    ($m2, $d2) = getmonthday($year, $sday+$i+1);

    chdir("$home_dir$dirs[$i]/run");
#    runcommand("/home/workspace/ccalmr/mazulauf/amb10xx/netcdf/cvs_stuff/forecasts/bin/atmos_nc/scripts/make_sflux_links.csh","3 $year $m1 $d1 $year $m2 $d2", 0);
    
    runcommand("cp", "$ref_dir/copy/* $home_dir$dirs[$i]/run/", 0);
    setup_bctides("$ref_dir/copy/bctides.in", "$home_dir$dirs[$i]/run/", $year, $m1, $d1);
}



##############################################

sub setup_bctides() {
# currently this thing only does the date: i.e. no tides.  See /home/users/hyde/scripts/setup_run.pl for code 
#  that does tides for the old param.in style

    my $pfile = $_[0];
    my $out_dir = $_[1];
    my $y = $_[2];
    my $m = $_[3];
    my $d = $_[4];
    my $jj;

    open(BCTIDES, "$pfile") or die "could not open $pfile\n";
    my @input = <BCTIDES>;
    close(BCTIDES);

    open(BCTIDES_OUT, ">$out_dir/bctides.in") or die "could not open $out_dir/bctides.in for writing\n";
    my @vals;
    
    my $temp=sprintf("%02d/%02d/%d 00:00:00 PST", $m, $d, $y);
    print BCTIDES_OUT "$temp\n";
    for ($jj=1; $jj<=$#input; ++$jj) {
	print BCTIDES_OUT "$input[$jj]";
    }
}

