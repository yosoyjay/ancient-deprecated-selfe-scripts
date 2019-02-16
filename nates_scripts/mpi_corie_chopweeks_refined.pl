#! /usr/bin/perl -w

# Divide the results from MPI SELFE for CORIE into weeks for processing
# Runs inside parent dir (i.e., one level above outputs/) on any system. 
# The time series must start from week 'start week' 
# but may not end at the end of a week. 
# The output individual weeks are inside parent dir.
@vars=('elev.61','pres.61','airt.61','shum.61','srad.61','flsu.61','fllu.61','radu.61','radd.61','flux.61','evap.61','prcp.61','wind.62','wist.62','dahv.62','vert.63','temp.63','salt.63','conc.63','tdff.63','vdff.63','kine.63','mixl.63','zcor.63','qnon.63','hvel.64');

if(@ARGV != 4) 
{
  print "Usage: $0 <year> <starting week # w/o leading 0> <# of days in total> <run id> \n";
  exit(1);
}

$year=$ARGV[0]; 
$sweek=$ARGV[1]; 
$ndays=$ARGV[2]; 
$runid=$ARGV[3];

$nweeks=int(($ndays-1)/7)+1; #no. of weeks
for($week=$sweek; $week<$sweek+$nweeks; $week++)
{
  $week1=sprintf('%02d',$week);
#  system "rm -rf $year-$week1-$runid/";
  if(!-e "$year-$week1-$runid")
  {
    system "mkdir $year-$week1-$runid/";
    system "mkdir $year-$week1-$runid/run";
    system "cd $year-$week1-$runid/run/; ln -sf ../../hgrid.* .";
    system "cd $year-$week1-$runid/run/; ln -sf ../../param.in .";
    system "cd $year-$week1-$runid/run/; ln -sf ../../vgrid.in .";
  }
}

foreach $var (@vars)
{
#  print "$var\n";
  if(-e "./outputs/1_$var")
  {
    for($day=1; $day<=$ndays; $day++)
    {
      $week=int(($day-1)/7)+$sweek;
      $week1=sprintf('%02d',$week);
      if(!-e "$year-$week1-$runid") {die "Dir not created: $year-$week1-$runid\n";}
      $new_day=$day-($week-$sweek)*7;
      system "cd $year-$week1-$runid/run/; ln -sf ../../outputs/$day\_$var $new_day\_$var";
    }
  } #if -e
} #foreach

