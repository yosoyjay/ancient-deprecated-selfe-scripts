#! /usr/bin/perl -w

#Drives combine_output*.f90 to combine all variables.
#Runs inside dir 'outputs'.

if(@ARGV != 4) 
  { 
    print "$0 <start_day # w/o leading 0s> <end_day #> <hotstart prefix> <n processors>\n";
    exit(1);
  }
$start_day=$ARGV[0]; 
$end_day=$ARGV[1];
$hstart=$ARGV[2];
$nprocs=$ARGV[3];

# Echo inputs
print "$0 $start_day $end_day\n";

$code = "/home/users/hyde/scripts/fortran/combine_output4_sirius";
$hstart_code = "/home/users/hyde/scripts/fortran/combine_hotstart2_sirius";
@vars=('elev.61','pres.61','airt.61','shum.61','srad.61','flsu.61','fllu.61','radu.61','radd.61','flux.61','evap.61','prcp.61','wind.62','wist.62','dahv.62','vert.63','temp.63','salt.63','conc.63','tdff.63','vdff.63','kine.63','mixl.63','zcor.63','qnon.63','hvel.64');

foreach $var (@vars)
{
  if(-e "$start_day\_0000_$var")
  {
    open(IN,">combine_output.in"); 
    print IN "$var\n $start_day $end_day\n";
    close(IN);
    system "$code";
    print "done combining $var...\n";
  } #if -e
} #foreach

open(IN, ">combine_hotstart.in");
print IN "$nprocs 0\n $hstart\n";
close(IN);
system "$hstart_code";
print "done combining $hstart"."_hotstart\n";
