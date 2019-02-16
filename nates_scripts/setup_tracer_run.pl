#!/usr/bin/perl

use strict;
use CORIE;

my $debug=0;

if ($#ARGV!=4) {
    print "$0 <source_dir> <target_dir> <codes_dir> <hotstart prefix> <n processors>\n";
    exit 0;
}

print "Note, this version assumes you're using forecast format.  <source_dir> is '2009-134', not '2009-134/run'.  It assumes the <source_dir> is the next day so that it can simply copy input files and then start from the target day.  This is because the data in run directory, say, 2009-134 actually starts at 2009-133.  Thus, the input data in 2009-135 actually starts on 2009-134.  Confused yet?\nDo not run this program in the background!!! It requires some input from you!\n";

my @command;
my $source_dir = $ARGV[0];
my $target_dir = $ARGV[1];
my $codes_dir = $ARGV[2];
my $hstart = $ARGV[3];
my $nprocs = $ARGV[4];

my $xy = '307319 321385';
my $secs = '49500'; 

# don't change these
my @bctidesin1 = ("  1 !itr_met", "  3 !nope", "  0 !itrtype", "  0 !itrtype", "  0 !itrtype");
# change as you will
my @bctidesin2 = ("  10.e6 !dye_mass in kg", "  0. !release depth in m (positive)", "  100. !dye_radius in meter", "  $xy !dye_x,dye_y", "  $secs  !dye_release in sec");

my $selfe = "pelfe_2.0hb_sirius";

my ($m1, $d1, $m2, $d2);

my $i;

#set up directory and 
if (!-e "$target_dir") {
    mkdir("$target_dir");
}
if (!-e "$target_dir/run") {
    mkdir("$target_dir/run");
}
if (!-e "$target_dir/run/outputs") {
    mkdir("$target_dir/run/outputs");
}
runcommand("cp","$source_dir/run/*.gr3 $target_dir/run/", 0);
runcommand("cp","$source_dir/run/*.th $target_dir/run/", 0); 
runcommand("cp","$source_dir/run/*.in $target_dir/run/", 0); 
runcommand("rm","$target_dir/run/hotstart.in", 0);  #making sure there's no confusion...
runcommand("ln","-s $source_dir/run/sflux/ $target_dir/run/", 0);
runcommand("cp","$codes_dir/$selfe $target_dir/run/", 0);
runcommand("cp","$codes_dir/runhind.qsub $target_dir/run/", 0);
runcommand("cp","$source_dir/run/*.ll $target_dir/run/", 0);


runcommand("ln","-s $source_dir/run/outputs/*hotstart* $target_dir/run/outputs/", 0);
runcommand("ln","-s $source_dir/run/outputs/*local* $target_dir/run/outputs/", 0);
chdir("$target_dir/run/outputs");
runcommand("$codes_dir/combine_hotstart2", "", 0);
# print STDIN "$nprocs 0\n";
# print STDIN "$hstart\n";
runcommand("cp","hotstart.in ..", 0);
runcommand("rm","*", 0);

# adjust param.in and bctides.in

open(PIN, "$target_dir/run/param.in") or die "could not open $target_dir/run/param.in\n";
my @pin = <PIN>;
close(PIN);
runcommand("mv", "$target_dir/run/param.in $target_dir/run/param.in.old", 0);
open(POUT, ">$target_dir/run/param.in") or die "could not open $target_dir/run/param.in for writing\n";
for ($i=0; $i<=$#pin; ++$i) {

    if ($pin[$i] =~ /ntracers\s*=/) {
	print POUT "  ntracers = 1\n";
    } elsif ($pin[$i] =~ /rnday\s*=/) {
	print POUT "  rnday = 2 !total run time in days\n";
    } elsif ($pin[$i] =~ /testout\s*=/) {
	print POUT "  trcr_1.63 = 1\n";
	print POUT "$pin[$i]";
    } else {
	print POUT "$pin[$i]";
    }

}
close(POUT);

open(BIN, "$target_dir/run/bctides.in") or die "could not open $target_dir/run/bctides.in\n";
my @bin = <BIN>;
close(BIN);
runcommand("mv", "$target_dir/run/bctides.in $target_dir/run/bctides.in.old", 0);
open(BOUT, ">$target_dir/run/bctides.in") or die "could not open $target_dir/run/bctides.in for writing\n";
for ($i=0; $i<=$#bin; ++$i) {
    print BOUT "$bin[$i]";
}
for ($i=0; $i<=$#bctidesin1; ++$i) {
    print BOUT "$bctidesin1[$i]\n";
}
for ($i=0; $i<=$#bctidesin2; ++$i) {
    print BOUT "$bctidesin2[$i]\n";
}
close(BOUT);

# create hotstart.in in the correct format
