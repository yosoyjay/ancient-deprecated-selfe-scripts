#!/usr/bin/perl
#
# Get and parse real time elevation from NOAA and USGS.
#
# Usage: getexdata.pl [db host name] [number of days from present to retrieve]
#
# There are current 2 agencies with data on web pages that are used by this
# application, NOAA and USGS. NOAA uses PST and MLLW in meters, USGS uses
# PST/PDT and
# stage height in feet. Both are converted to NGVD29 using conversions
# located in the forecast/sites table.
#
# Opens database forecasts, reads table sites and writes to table
# elevation.
#

#use strict;
use CORIE;

# Clean up this query
my $siteid;
my $refdepth;
my ( $endyr, $endmo );

my ( $se, $mi, $hr, $da, $mo, $yr, $wday, $yday, $isdst ) = localtime( time - $ndays * 86400 );

$mo += 1;
$yr += 1900;
my $currentyr = $yr;
my @stations  = ( 'tpoin', 'neah', 'tokepoint', 'westport', 'lapush', 'garibaldi', 'southbeach' );
#my @stations  = ( 'tpoin', 'neah', 'tokepoint' );
my @ids       = ( '9439040', '9443090', '9440910', '9441102', '9442396', '9437540', '9435380' );
my @starts    = ( 1999, 1999, 1999, 1999, 1999, 1999, 1999 );
my @datums    = ( 4, 4, 4, 4, 4, 4, 4 );
my @dshifts   = ( 1.028, 1.136, 1.101, 1.105, 1.095, 0, 1.009 );

my @ndays    = ( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );

for ( $i = 0 ; $i < @stations ; $i++ ) {
#for ( $i = 0 ; $i < 1 ; $i++ ) {
    my $station = $stations[$i];
    my $id      = $ids[$i];
    my $datum   = $datums[$i];
    my $start   = $starts[$i];
    my $cor   = $dshifts[$i];

# http://tidesandcurrents.noaa.gov/data_listing.shtml?bdate=20050101&edate=20050131&wl_sensor_hist=W1&relative=&datum=$datum&unit=0&shift=s&stn=9440910+Toke+Point%2C+WA&type=Historic+Tide+Data&format=View+Data

    my $baseurl = "http://tidesandcurrents.noaa.gov/data_listing.shtml?bdate=STARTDATE&edate=ENDDATE&wl_sensor_hist=W1&relative=&datum=$datum&unit=0&shift=s&stn=$id&type=Historic+Tide+Data&format=View+Data&listing=1";
    my $outfile = '-';
    my $cday;

    for ( my $yy = $starts[0] ; $yy <= $starts[0] ; $yy++ ) {
        #for ( my $mm = 2 ; $mm <= 3 ; $mm++ ) {
        for ( my $mm = 1 ; $mm <= 12 ; $mm++ ) {
			open( OUT, ">$station" . "_" . $mm ."_1999.txt" );

            $starttime = sprintf( "%4d%02d%02d", $yy, $mm, 1 );
            $endtime   = sprintf( "%4d%02d%02d", $yy, $mm, $ndays[$mm - 1] );
            $url       = $baseurl;
            $url =~ s/STARTDATE/$starttime/;
            $url =~ s/ENDDATE/$endtime/;
            print "Getting $url\n";
            my @elev = `wget -o /dev/null -O $outfile \"$url\"`;
            foreach $line (@elev) {
                chomp $line;
                if ( $line =~ /Station\s+Date\s+Time\s+Pred 6$/ ) {
                    print "file contains no observations.\n";
                    next;
                }
                if ( $line =~ /^\d/ ) {
                    $line =~ s///;
                    $line =~ s/\s{2,4}/ /g;
                    ( $siteid, $d1, $d2, $pred, $val ) = split ( /\s/, $line );
                    if ( $d1 !~ /^\d/ || $val !~ /^-{0,1}\d/ ) {
                        next;
                    }

                    $val  = sprintf( "%.3f", $val * 1.0 - $cor );
                    $pred = sprintf( "%.3f", $pred * 1.0 - $cor );
                    $yr = substr $d1, 0, 4;
                    $mo = substr $d1, 4, 2;
                    $da = substr $d1, 6, 2;
                    ( $hr, $mi, $se ) = split ( /:/, $d2 );
                    $yr *= 1.0;
                    $mo *= 1.0;
                    $da *= 1.0;
                    $hr *= 1.0;
                    $mi *= 1.0;
                    $se *= 1.0;
                    $cday = "$se $mi $hr $da $mo $yr";
                    $cday = getcorieday($se, $mi, $hr, $da, $mo, $yr);
                    print OUT "$mo-$da-$yr $hr:$mi:$se, $val\n";
                }    # if data line
            }    # line
		}    # month
		close(OUT);
    }    # year
}    #station
