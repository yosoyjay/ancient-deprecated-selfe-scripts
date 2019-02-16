#!/usr/bin/perl

use strict;
use CORIE;

if ($#ARGV!=1) {
    print "<source hgrid> <out_dir>\n";
    exit 0;
}


my $hgrid = $ARGV[0];
my $out_dir = $ARGV[1];

open(HGRID, $hgrid) or die "could not open $hgrid for reading\n";
my @in = <HGRID>;
close(HGRID);

my ($ne, $nn) = split(/\s+/, $in[1]);
my $i;
my @vals;

open(INTERPOL, ">$out_dir/interpol.gr3") or die "could not open $out_dir/interpol.gr3 for writing\n";
print(INTERPOL "interpol\n");
print(INTERPOL $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    if ($vals[3] > 100) {
	print INTERPOL "$vals[0] $vals[1] $vals[2] 1\n";
    } else {
	print INTERPOL "$vals[0] $vals[1] $vals[2] 2\n";
    }
}

for ($i=0; $i<$ne; ++$i) {
    print INTERPOL $in[$nn+2+$i];
}

close(INTERPOL);

############################################################
############# based on /home/users/yinglong/Scripts/gen_nudge.f90

## center pt:
my $x0=333906;
my $y0=240000;  #293427;

## Minor and major axes
my $rx = 170000;
my $ry = 370000;  ## outer_rx/inner_rx, same for ry, must by > 1
my $rat = 1.7;
my $rmax = 90/1.2/86400;

open(NUDGE, ">$out_dir/nudge.gr3") or die "could not open $out_dir/new_nudge.gr3 for writing\n";
print(NUDGE "nudge\n");
print(NUDGE $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);

    my $rr = ($vals[1]-$x0)**2/($rx**2)+($vals[2]-$y0)**2/($ry**2);

    my $tnu = $rmax*($rr-1)/($rat*$rat-1);

#    my $ret = $rmax;
#    if ($rmax>$tnu) {
#	$ret = $tnu;
#    }
#    if ($ret<0) {
#	$ret = 0;
#    }

    $tnu = max(0, min($rmax, $tnu));
    print NUDGE "$vals[0] $vals[1] $vals[2] $tnu\n";
}

for ($i=0; $i<$ne; ++$i) {
    print NUDGE $in[$nn+2+$i];
}

close(NUDGE);

runcommand("cp", "$out_dir/nudge.gr3 $out_dir/s_nudge.gr3", 0);
runcommand("cp", "$out_dir/nudge.gr3 $out_dir/t_nudge.gr3", 0);

############################################## constant .gr3s #####################################

my $xlsc = 0.5;
open(XLSC, ">$out_dir/xlsc.gr3") or die "could not open $out_dir/xlsc.gr3 for writing\n";
print(XLSC "xlsc = $xlsc\n");
print(XLSC $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print XLSC "$vals[0] $vals[1] $vals[2] $xlsc\n";
}

for ($i=0; $i<$ne; ++$i) {
    print XLSC $in[$nn+2+$i];
}
close(XLSC);

my $watertype = 7;
open(WATERTYPE, ">$out_dir/watertype.gr3") or die "could not open $out_dir/watertype.gr3 for writing\n";
print(WATERTYPE "watertype = $watertype\n");
print(WATERTYPE $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print WATERTYPE "$vals[0] $vals[1] $vals[2] $watertype\n";
}

for ($i=0; $i<$ne; ++$i) {
    print WATERTYPE $in[$nn+2+$i];
}
close(WATERTYPE);

my $albedo = 0.06;
open(ALBEDO, ">$out_dir/albedo.gr3") or die "could not open $out_dir/albedo.gr3 for writing\n";
print(ALBEDO "albedo = $albedo\n");
print(ALBEDO $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print ALBEDO "$vals[0] $vals[1] $vals[2] $albedo\n";
}

for ($i=0; $i<$ne; ++$i) {
    print ALBEDO $in[$nn+2+$i];
}
close(ALBEDO);

my $diffmax = 1.0;
open(DIFFMAX, ">$out_dir/diffmax.gr3") or die "could not open $out_dir/diffmax.gr3 for writing\n";
print(DIFFMAX "diffmax = $diffmax\n");
print(DIFFMAX $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print DIFFMAX "$vals[0] $vals[1] $vals[2] $diffmax\n";
}

for ($i=0; $i<$ne; ++$i) {
    print DIFFMAX $in[$nn+2+$i];
}
close(DIFFMAX);

my $diffmin = 0.000001;
open(DIFFMIN, ">$out_dir/diffmin.gr3") or die "could not open $out_dir/diffmin.gr3 for writing\n";
print(DIFFMIN "diffmin = $diffmin\n");
print(DIFFMIN $in[1]);

for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);
    print DIFFMIN "$vals[0] $vals[1] $vals[2] $diffmin\n";
}

for ($i=0; $i<$ne; ++$i) {
    print DIFFMIN $in[$nn+2+$i];
}
close(DIFFMIN);

####################################################################
################## do some of the more complex ones... #############
###### hgrid.ll ######

open(TMP, ">$out_dir/tmp.txt") or die "could not open $out_dir/tmp.txt for writing\n";
for ($i=2; $i<=$nn+1; ++$i) {
    print(TMP $in[$i]);
}
close(TMP);

unlink("ll.tmp");

open(SPCS_IN, "| /home/users/hyde/bin/spcs.pl");
print SPCS_IN "tmp.txt\n";
print SPCS_IN "ll.tmp\n";
print SPCS_IN "2\n";
print SPCS_IN "1\n";
print SPCS_IN "8\n";
print SPCS_IN "1\n";
print SPCS_IN "3601\n";
print SPCS_IN "wo\n";
print SPCS_IN "2\n";
close SPCS_IN;

open(TMP, "ll.tmp") or die "could not open ll.tmp, an intermediate file in hgrid.ll generation\n";
my @ll = <TMP>;
close(TMP);
open(HGRIDLL, ">$out_dir/hgrid.ll") or die "could not opne hgrid.ll for writing";
print(HGRIDLL "hgrid.ll\n$ne $nn\n");
for ($i=0; $i<=$nn+1; ++$i) {
    print(HGRIDLL $ll[$i]);
}
for ($i=0; $i<$ne; ++$i) {
    print HGRIDLL $in[$nn+2+$i];
}
close(HGRIDLL);

##### windrot_geo2proj.gr3 #####
runcommand("/home/users/cseaton/bin/rotate_wind/rotate_wind_spcs2ll","-input $out_dir/hgrid.ll -output $out_dir/windrot_geo2proj.gr3 -ll2spcs", 0);

##### drag.gr3 #####

my @estuary_region = (321643, 302813, 531393, 177637);
my @upstream_region = (403060, 246529, 510757, 182296);
##### drag.gr3 #####
my $basedrag = 0.002;
my $intermediatedrag = 0.001;
my $upstreamdrag = 0.0005;
my $slope = -0.00000008244;
my $intercept = 0.0307;

open(DRAG, ">$out_dir/drag.gr3");
print DRAG "drag.gr3\n$ne $nn\n";
for ($i=0; $i<$nn; ++$i) {
    my @vals = split(/\s+/, $in[$i+2]);

    my $dragnu;

    if ($vals[1] > $estuary_region[0] && $vals[1] < $estuary_region[2] && $vals[2] < $estuary_region[1] && $vals[2] > $estuary_region[3]) {
	if ($vals[1] > $upstream_region[0] && $vals[1] < $upstream_region[2] && $vals[2] < $upstream_region[1] && $vals[2] > $upstream_region[3]) {
	    $dragnu = $upstreamdrag;
	    

	} else {
	    $dragnu = min(max($slope*$vals[1] + $intercept, $intermediatedrag), $basedrag);
	}
    } else {
	$dragnu = $basedrag;
    }


    print DRAG "$vals[0] $vals[1] $vals[2] $dragnu\n";
}
for ($i=0; $i<$ne; ++$i) {
    print DRAG $in[$nn+2+$i];
}

close(DRAG);







