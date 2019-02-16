#!/usr/bin/perl

use strict;
use CORIE;
use POSIX qw(ceil floor);

if ($#ARGV!=5) {
    print "$0 <source dir> <target dir> <year> <julian start day> <julian end day> <source db>\n";
    exit 1;
}

my $base_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $db = $ARGV[5];

my ($i, $j, $k);

my @link_files = ("elev.61", "hvel.64", "temp.63", "salt.63", "zcor.63", "wind.62"); #, "trcr_1.63", "trcr_2.63");

for ($i=$sday; $i<=$eday; ++$i) {
    my $dir=sprintf("%d-%03d", $year, $i);

    my $w = floor(($i-1)/7)+1;
    my $local_day = ($i-1) % 7 + 1;

    my $db_dir = sprintf("%d-%02d-%s", $year, $w, $db);
    print "doing $db_dir\n";

    if (-e "$base_dir/$db_dir/run") {

	runcommand("ln", "-s $base_dir/$db_dir/run/hgrid.gr3 $target_dir/hgrid.gr3", 1);
	runcommand("ln", "-s $base_dir/$db_dir/run/hgrid.gr3 $target_dir/hgrid.gr3", 0);

	for ($j=0; $j<=$#link_files; ++$j) {
#	runcommand("ln", "-s $base_dir/$db_dir/run/$local_day"."_$link_files[$j] $target_dir/$dir/run/1_$link_files[$j]", 1);
	    runcommand("ln", "-s $base_dir/$db_dir/run/$local_day"."_$link_files[$j] $target_dir/".($i-$sday+1)."_$link_files[$j]", 0);
	}
    } else {
	die "$base_dir/$db_dir/run does not exist"
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
