#!/usr/bin/perl

use CORIE;

if ($#ARGV < 3) {
    print "usage: rasterAscii2xyz.pl <fin> <fout> <subsample factor> <max elev>\n";
    exit 1;
}

$fin = $ARGV[0];
$fout = $ARGV[1];
$subsample_factor = $ARGV[2];
$max_elev = $ARGV[3];

open(IN, "$fin") or die "could not open $fin\n";;
open(OUT, ">$fout") or die "could not open $fout for writing\n";;

$line = <IN>;
@vals = split(/\s+/, $line);
$ncols = $vals[1];

$line = <IN>;
@vals = split(/\s+/, $line);
$nrows = $vals[1];

$line = <IN>;
@vals = split(/\s+/, $line);
$xllcorner = $vals[1];

$line = <IN>;
@vals = split(/\s+/, $line);
$yllcorner = $vals[1];

$line = <IN>;
@vals = split(/\s+/, $line);
$cellsize = $vals[1];

$line = <IN>;
@vals = split(/\s+/, $line);
$nodata = $vals[1];

print "cellsize: $cellsize\n";

$iline = 1;
while (!eof (IN)) {

#    print "$iline of $nrows\n";

    $line = <IN>;
    @vals = split(/\s+/, $line);

#    $yloc = $yllcorner + ($nrows-$iline)*$cellsize + $cellsize/2 + 20;
    $yloc = $yllcorner + ($nrows-$iline)*$cellsize + $cellsize/2;

    if ($iline % $subsample_factor==0) {

	for ($j=1; $j<=scalar(@vals); ++$j) {
#	$irand = int(rand($subsample_factor));
	
#	    $xloc = $xllcorner + $j*$cellsize + $cellsize/2 + 80;
	$xloc = $xllcorner + $j*$cellsize + $cellsize/2;

#	if ($vals[$j-1] != $nodata  &&  $irand % $subsample_factor == 0) {
	    if ($vals[$j-1] != $nodata  &&  $j % $subsample_factor == 0 && $vals[$j-1] <= $max_elev) {
		print OUT "$xloc $yloc ".($vals[$j-1])."\n";
	    }
	}

    }

    $iline = $iline+1;
}
