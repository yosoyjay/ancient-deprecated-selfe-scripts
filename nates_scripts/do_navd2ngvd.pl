#!/usr/bin/perl

## The frefpts file has a selection of points converted form NAVD88 to NGVD27 using the on-line version of VERTCON.
## Vertical translation is done be finding the closest node from this file.
## An example is /home/users/hyde/models/columbia/grids/navd2ngvdPtsOSPN.txt

use CORIE;

if ($#ARGV != 2) {
    print "usage: perl do_ngvd2navd.pl <input file> <output file> <ft_meter_factor -- 0.3048 or 3.2808>\n";
    exit 1;
}

$fin = $ARGV[0];
$fout = $ARGV[1];
$factor = $ARGV[2];

$frefpts = "/home/users/hyde/models/columbia/grids/navd2ngvdPtsOSPN.txt";

open (REF, $frefpts) or die "could not open $frefpts\n";
@lines = <REF>;
close(REF);
for ($i=1; $i<=$#lines; ++$i) {
    @vals = split(/\s+/, $lines[$i]);
    $x[$i-1] = $vals[1];
    $y[$i-1] = $vals[2];
    $zdiff[$i-1] = $vals[3];    
}

open(IN, $fin) or die "could not open $fin";
open(OUT, ">$fout") or die "could not open $fout\n";

while (!eof (IN)) {
    $line = <IN>;
    @vals = split(/,/, $line);
    
    if ($#vals !=2) {
	@vals = split(/\s+/, $line);
    }

    $curr_diff = 9999999;
    $curr_index = 9999999;
    for ($i=0; $i<=$#x; ++$i) {
	$dist = distance($vals[0], $x[$i], $vals[1], $y[$i]);
	if ($dist < $curr_diff) {
	    $curr_diff = $dist;
	    $curr_index = $i;
	}
    }

    if ($#vals==2 && $vals[1] =~ /[0-9]/) {
	#### a minus sign if positive up, plus if positive down
	print OUT "$vals[0], $vals[1], ".($vals[2]*$factor-$zdiff[$curr_index])."\n";
    }
}

close(IN);
close(OUT);



############## subroutines ###############

sub distance {
    $xdiff = $_[0]-$_[2];
    $ydiff = $_[1]-$_[3];

    return sqrt($xdiff*$xdiff+$ydiff*$ydiff);
}


