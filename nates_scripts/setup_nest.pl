#!/usr/bin/perl

use strict;
#use CORIE;
use POSIX qw(ceil floor);

my $debug=0;

if ($#ARGV!=6) {
    print "$0 <base_dir> <bounding model dir> <year> <julian start day> <julian end day> <ref dir> <time step>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0]."/";
my $bound_model = $ARGV[1];
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $ref_dir = $ARGV[5];
my $tstep = $ARGV[6];

my @dirs;
my ($m1, $m2, $d1, $d2);

#constants or inputs by another form...
my $nest_db = 16;  # -1 for forecast format
my $model_nvrt = 6;
my $nope_line = "2 42 0";
my $new_dt = $tstep;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$home_dir$dirs[$#dirs]\n";

    if (!-e "$home_dir$dirs[$#dirs]") {
	mkdir("$home_dir$dirs[$#dirs]");
    }
    if (!-e "$home_dir$dirs[$#dirs]/run") {
	mkdir("$home_dir$dirs[$#dirs]/run");
    }
}

for ($i=0; $i<=$#dirs; ++$i) {
    chdir("$ref_dir");
#    runcommand("rm","fort.*",0);
    setup_nest("$bound_model", "$home_dir$dirs[$i]/run/", "$ref_dir/", $year, $sday+$i, $nest_db, $model_nvrt, $nope_line, $new_dt);
}

##############################################

sub setup_nest() {
    my $source_dir = $_[0];
    my $target_dir = $_[1]; 
    my $ref_dir = $_[2];
    my $y = $_[3];
    my $d = $_[4];
    my $db = $_[5];  # -1 for forecast format
    my $nvrt = $_[6];
    my $nope_line = $_[7];
    my $new_dt = $_[8];

    my $exe = "interpolate_variables_selfe4_1day";
    my $timeint = "timeint_3Dth_2";

    my $model_dir;
    my $local_day;
    
    if ($db == -1) {
	$model_dir = sprintf("%d-%03d", $year, $d);
	$local_day = 1;
    } else {
	my $w = floor(($d-1)/7)+1;
	$local_day = ($d-1) % 7 + 1;

	$model_dir = sprintf("%d-%02d-%02d", $year, $w, $db);
    }

    print "$y-$d == $model_dir, $local_day\n";

    runcommand("rm", "$ref_dir/*.61", $debug);
    runcommand("rm", "$ref_dir/*.63", $debug);
    runcommand("rm", "$ref_dir/*.64", $debug);

    runcommand("ln", "-s $source_dir$model_dir/run/$local_day"."_elev.61 $ref_dir/$local_day"."_elev.61", $debug);
    runcommand("ln", "-s $source_dir$model_dir/run/$local_day"."_salt.63 $ref_dir/$local_day"."_salt.63", $debug);
    runcommand("ln", "-s $source_dir$model_dir/run/$local_day"."_temp.63 $ref_dir/$local_day"."_temp.63", $debug);
    runcommand("ln", "-s $source_dir$model_dir/run/$local_day"."_hvel.64 $ref_dir/$local_day"."_hvel.64", $debug);
    

    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
    print INTERP "1 $local_day\n";
    close INTERP;

    runcommand("$ref_dir/$exe", "", $debug);
    runcommand("cp", "fort.16 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "1 1\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("$ref_dir/$timeint", "", $debug);
    runcommand("cp", "th.new $target_dir/elev3D.th", $debug);

    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
    print INTERP "2 $local_day\n";
    close INTERP;
    runcommand("$ref_dir/$exe", "", $debug);
    runcommand("cp", "fort.17 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "1 $nvrt\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("$ref_dir/$timeint", "", $debug);
    runcommand("cp", "th.new $target_dir/temp3D.th", $debug);
    
#    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
#    print INTERP "2 $local_day\n";
#    close INTERP;
#    runcommand("$ref_dir/interp_selfe3_1_day", "", 0);
    runcommand("cp", "fort.18 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "1 $nvrt\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("$ref_dir/$timeint", "", $debug);
    runcommand("cp", "th.new $target_dir/salt3D.th", $debug);

    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
    print INTERP "3 $local_day\n";
    close INTERP;
    runcommand("$ref_dir/$exe", "", $debug);
    runcommand("cp", "fort.17 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "1 ".(2*$nvrt)."\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("$ref_dir/$timeint", "", $debug);
    runcommand("cp", "th.new $target_dir/uv3D.th", $debug);
}


sub runcommand {
    my ( $e, $a, $debug ) = @_;

    if ($debug) { # in case $debug is undef
    } else {
       $debug = 0;
    }

    # invoke the real system call here
    my $err = 0;
    if ($debug) {
        print "$e $a\n";
    } else {
        $err = system("$e $a");
    }

    return $err;
}
