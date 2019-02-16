#!/usr/bin/perl

use strict;
use CORIE;

my $base_dir = "/home/workspace/ccalmr38/hindcasts/";

my @runs = ("yaqalsea_hindcast_2", "coos_hindcast_1", "willapa_calib4");
my @shapes_dir = ("yaqalsea_files/", "coos_hindcast_1/shapes/", "willapa_files/");
my $dbg = 1;
my $key;

my ($i, $j);

for ($i=0; $i<=$#runs; ++$i) {
    my %year_days;

    opendir(DIR, "$base_dir$runs[$i]/");
    my @names = readdir(DIR);
    closedir(DIR);

    for($j=0; $j<=$#names; ++$j) {
#	print "$names[$j]\n";

	if (-e "$base_dir/$runs[$i]/$names[$j]/run/2_elev.61") {
	    if (-e "$base_dir/$runs[$i]/$names[$j]/process/cross_section_data.dat"){}
	    else {
		my ($y, $d) = split(/-/, $names[$j]);

		print "y: $y, d: $d\n";

		my $key;
		my $found=0;
		foreach $key (keys %year_days) {
		    if ($key==$y) {
			$found=1;
		    }
		}
		if (!$found) {
		    @{ $year_days{$y} }= ($d);
		}
		else {
		    my @temp = $year_days{$y};
		    $year_days{$y}[@temp] = $d;
		}
	    }
	}
    }
    my @tkey = keys %year_days;
    print "n keys: ".scalar(@tkey)."\n";

    foreach $key (keys %year_days) {
	print "key: $key, year_days: $year_days{$key}\n";
	foreach my $d (@{ $year_days{$key} }) {
	    print "\t$d\n";
	}
	my @temp = sort(@{$year_days{$key}});
	runcommand("python", "/home/hyde/ourestuaries/pythonlib/cmop/do_daily_extract.py $runs[$i] $base_dir base_dir/$shapes_dir[$i] $key $temp[0] ".($temp[$#temp]-$temp[0]+1)." 0", $dbg);
    }
}