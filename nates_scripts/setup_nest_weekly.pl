#!/usr/bin/perl

use strict;
#use CORIE;
use POSIX qw(ceil floor);

my $debug=0;

if ($#ARGV!=7) {
    print "$0 <base_dir> <new_db> <bounding model dir> <year> <julian start day> <julian end day> <ref dir> <time step>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0]."/";
my $new_db = $ARGV[1];
my $bound_model = $ARGV[2];
my $year = $ARGV[3];
my $sday = $ARGV[4];
my $eday = $ARGV[5];
my $ref_dir = $ARGV[6];
my $tstep = $ARGV[7];

my @dirs;
my ($m1, $m2, $d1, $d2);

#constants or inputs by another form...
my $nest_db = 14;  # -1 for forecast format
my $model_nvrt = 26;
my $nope_line = "1 32";
my $new_dt = $tstep;

my $i;
for ($i=$sday; $i<=$eday; $i=$i+7) {
    my $week=floor($i/7)+1;
    $dirs[$#dirs+1]=sprintf("%d-%02d-%s", $year, $week, $new_db);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    chdir("$ref_dir");
#    runcommand("rm","fort.*",0);
    setup_nest("$bound_model", "$home_dir$dirs[$i]/run/", "$ref_dir/", $year, floor(($sday+$i*7)/7)+1, $nest_db, $model_nvrt, $nope_line, $new_dt);
}

##############################################

sub setup_nest() {
    my $source_dir = $_[0];
    my $target_dir = $_[1]; 
    my $ref_dir = $_[2];
    my $y = $_[3];
    my $w = $_[4];
    my $db = $_[5];  # -1 for forecast format
    my $nvrt = $_[6];
    my $nope_line = $_[7];
    my $new_dt = $_[8];

    my $debug = 0;

    my $model_dir;
    my $local_day;
    my $ii;

    $model_dir = sprintf("%d-%02d-%02d", $year, $w, $db);

#    print "$y-$d == $model_dir, $local_day\n";

    runcommand("rm", "$ref_dir/*.61", $debug);
    runcommand("rm", "$ref_dir/*.63", $debug);
    runcommand("rm", "$ref_dir/*.64", $debug);

    for ($ii=1; $ii<8; ++$ii) {
        runcommand("ln", "-s $source_dir$model_dir/run/$ii"."_elev.61 $ref_dir/$ii"."_elev.61", $debug);
        runcommand("ln", "-s $source_dir$model_dir/run/$ii"."_salt.63 $ref_dir/$ii"."_salt.63", $debug);
        runcommand("ln", "-s $source_dir$model_dir/run/$ii"."_temp.63 $ref_dir/$ii"."_temp.63", $debug);
        runcommand("ln", "-s $source_dir$model_dir/run/$ii"."_hvel.64 $ref_dir/$ii"."_hvel.64", $debug);
    }    

    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
    print INTERP "1 7\n";
    close INTERP;

    runcommand("$ref_dir/interpolate_variables_selfe3", "", $debug);
    runcommand("cp", "fort.16 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "7 1\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("./timeint_3Dth", "", $debug);
    runcommand("cp", "th.new $target_dir/elev3D.th", $debug);

    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
    print INTERP "2 7\n";
    close INTERP;
    runcommand("$ref_dir/interpolate_variables_selfe3", "", $debug);
    runcommand("cp", "fort.17 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "7 $nvrt\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("./timeint_3Dth", "", $debug);
    runcommand("cp", "th.new $target_dir/temp3D.th", $debug);
    
    runcommand("cp", "fort.18 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "7 $nvrt\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("./timeint_3Dth", "", $debug);
    runcommand("cp", "th.new $target_dir/salt3D.th", $debug);

    open(INTERP, ">$ref_dir/interpolate_variables_selfe.in") or die "could not open $ref_dir/interpolate_variables_selfe.in for writing\n";
    print INTERP "3 7\n";
    close INTERP;
    runcommand("$ref_dir/interpolate_variables_selfe3", "", $debug);
    runcommand("cp", "fort.17 th.old", $debug);
    open(TIMEINT, ">$ref_dir/timeint.in") or die "could not open $ref_dir/timeint.in for writing\n";
    print TIMEINT "7 ".(2*$nvrt)."\n$nope_line\n$new_dt\n";
    close TIMEINT;
    runcommand("./timeint_3Dth", "", $debug);
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
