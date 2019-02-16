#!/usr/bin/perl

use strict;
use CORIE;
use POSIX qw(ceil floor);

if ($#ARGV!=4) {
    print "$0 <base_dir> <target dir> <year> <julian start day> <julian end day>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $local_day;
my ($i, $j, $k);
my $dir;

my @link_files = ("elev.61", "hvel.64", "temp.63", "salt.63", "zcor.63");

for ($i=$sday; $i<=$eday; ++$i) {
    my $from_dir=sprintf("%d-%03d", $year, $i);

    for ($j=-1; $j<=1; ++$j) {
	my $to_dir = sprintf("%d-%03d", $year, $i+$j);
	if (!-e "$target_dir/$to_dir") {
	    mkdir("$target_dir/$to_dir");
	}
	if (!-e "$target_dir/$to_dir/run") {
	    mkdir("$target_dir/$to_dir/run");
	}

	runcommand("rm", "$target_dir/$to_dir/run/hgrid.gr3 $target_dir/$to_dir/run/1_*", 1);
	runcommand("rm", "$target_dir/$to_dir/run/hgrid.gr3 $target_dir/$to_dir/run/1_*", 0);
	runcommand("ln", "-s $base_dir/$from_dir/run/hgrid.gr3 $target_dir/$to_dir/run/hgrid.gr3", 1);
	runcommand("ln", "-s $base_dir/$from_dir/run/hgrid.gr3 $target_dir/$to_dir/run/hgrid.gr3", 0);	
	
	$local_day = $j+2;

	for ($k=0; $k<=$#link_files; ++$k) {
	    runcommand("ln", "-s $base_dir/$from_dir/run/$local_day"."_$link_files[$j] $target_dir/$to_dir/run/1_$link_files[$j]", 1);
	    runcommand("ln", "-s $base_dir/$from_dir/run/$local_day"."_$link_files[$j] $target_dir/$to_dir/run/1_$link_files[$j]", 0);
	}
    }
}

###############################################

sub get_daily_flux() {
    my $fname = $_[0];
    my $day = $_[1];
    
    open(INPUT, "$fname") or die "could not open $fname for reading\n";
    my @input = <INPUT>;

    my ($ii, $jj);
    my @vals;
    my @ret;

    @vals = split(/\s+/, $input[0]);
    my $dt = $vals[0];
    my $npts = 86400/$dt;
    my $offset = ($day-1)*$npts;

    for ($ii=0; $ii<$npts; ++$ii) {
	@vals = split(/\s+/, $input[$ii+$offset]);
	$ret[$ii][0]=$dt*($ii+1);
	for ($jj=1; $jj<=$#vals; ++$jj) {
	    $ret[$ii][$jj] = $vals[$jj];
	}
    }

    return @ret;
}
