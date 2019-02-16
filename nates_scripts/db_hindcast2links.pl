#!/usr/bin/perl

use strict;
use CORIE;
use POSIX qw(ceil floor);

if ($#ARGV!=5) {
    print "$0 <source dir> <target dir> <source_day_1> <source_last_day> <year> <julian start day>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $base_start = $ARGV[2];
my $base_end = $ARGV[3];
my $target_year = $ARGV[4];
my $target_day_start = $ARGV[5];

my ($i, $j, $k);

my @link_files = ("elev.61", "hvel.64", "temp.63", "salt.63", "zcor.63", "wind.62"); #, "trcr_1.63", "trcr_2.63");

for ($i=$base_start; $i<=$base_end; ++$i) {

    my $target_day = $target_day_start+$i-$base_start;
    my $dir=sprintf("%d-%03d", $target_year, $target_day);

    print "doing $base_dir / day $i\n";

    if (!-e "$target_dir/$dir") {
	mkdir("$target_dir/$dir");
    }
    if (!-e "$target_dir/$dir/run") {
	mkdir("$target_dir/$dir/run");
    }   

    runcommand("ln", "-sf $base_dir/hgrid.gr3 $target_dir/$dir/hgrid.gr3", 1);
    runcommand("ln", "-sf $base_dir/hgrid.gr3 $target_dir/$dir/hgrid.gr3", 0);

    for ($j=0; $j<=$#link_files; ++$j) {
	runcommand("ln", "-sf $base_dir/outputs/$i"."_$link_files[$j] $target_dir/$dir/run/1_$link_files[$j]", 1);
	runcommand("ln", "-sf $base_dir/outputs/$i"."_$link_files[$j] $target_dir/$dir/run/1_$link_files[$j]", 0);
    }
}

###############################################

# sub get_daily_flux() {
#     my $fname = $_[0];
#     my $day = $_[1];
#     
#     open(INPUT, "$fname") or die "could not open $fname for reading\n";
#     my @input = <INPUT>;
# 
#     my ($ii, $jj);
#     my @vals;
#     my @ret;
# 
#     @vals = split(/\s+/, $input[0]);
#     my $dt = $vals[0];
#     my $npts = 86400/$dt;
#     my $offset = ($day-1)*$npts;
# 
#     for ($ii=0; $ii<$npts; ++$ii) {
# 	@vals = split(/\s+/, $input[$ii+$offset]);
# 	$ret[$ii][0]=$dt*($ii+1);
# 	for ($jj=1; $jj<=$#vals; ++$jj) {
# 	    $ret[$ii][$jj] = $vals[$jj];
# 	}
#     }
# 
#     return @ret;
# }
