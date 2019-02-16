#! /usr/bin/perl -w
#
# $Revision: 1.2 $
# $Log: do_hind_process.pl,v $
# Revision 1.2  2006/03/14 15:53:28  pturner
# Changes.
#
# Revision 1.1.1.1  2004/11/15 17:28:52  pturner
# Selfe processing.
#
# Revision 1.5  2004/10/26 13:47:36  pturner
# Fixups.
#
# Revision 1.4  2004/10/11 23:18:15  pturner
# Adding baseDir argument.
#
# Revision 1.3  2004/10/08 04:40:34  pturner
# More base dir changes.
#
# Revision 1.2  2004/10/08 00:09:17  gliu
# check in for sigma grid
#
# Revision 1.1.1.1  2004/09/22 17:00:21  gliu
# Hindcast processing version 06
#
#
#

use strict;
use Cwd;
use File::Basename;
use CORIE;

#my $cwd = getcwd();
#my ($name, $path, $suffix) = fileparse($0, ("pl"));
#my $name = basename($0);
#my $path = dirname($0);
#if ($path eq "." || $path eq "./") {
#    $path = $cwd;
#} 
#print "$0: $name, $path\n";
#exit(1);

if (@ARGV != 5)
{ 
  print "USAGE: $0 <base dir> <input dir> <year> <week> <o/n>\n";  
  print "     input dir: contains elcirc input and output files\n";
  print "     base dir: directory which has this script\n";
  print "     year: the actual year of the run\n";
  print "     week: the week number of the run\n";
  print "     o/n: official run or not\n";

  exit(1);
}

my $debug = 0;
my $err = 0; 

my $inputDir = $ARGV[1];
my $outputDir = $ARGV[1];
if (!-e "$outputDir/log")
{ $err = runcommand("/bin/mkdir", "$outputDir/log", $debug); }

my $xml = "$inputDir/log/output.xml";
setXMLOutputFile($xml);

my $baseDir = $ARGV[0];

my $runYear = $ARGV[2];
my $runWeek = $ARGV[3];

my $official = $ARGV[4];
my $official_flg;
if ($official eq "o")
{ $official_flg = 1; }
else
{ $official_flg = 0; }

my $i; my @tempset=();

my $startDay; my $endDay; 
  
$startDay = ($runWeek-1)*7+1;
$startDay = getcoriedayfromyearday(0,0,0,$startDay,$runYear);

if ($runWeek == 53)
{
  if ($runYear % 4 == 0) { $endDay = $startDay + 2; }
  else { $endDay = $startDay + 1; }
}
else
{ $endDay = $startDay + 7; } 

#$err = runcommand("/bin/rm", "-rf $outputDir/images", $debug);
#$err = runcommand("/bin/rm", "-rf $outputDir/process", $debug);
#$err = runcommand("/bin/rm", "-rf $outputDir/data_files", $debug);

# to prepare the directories for output
if (!-e "$outputDir")
{ $err = runcommand("/bin/mkdir", "$outputDir", $debug); }

if (!-e "$outputDir/process")
{ $err = runcommand("/bin/mkdir", "$outputDir/process", $debug); }

if (!-e "$outputDir/images")
{ $err = runcommand("/bin/mkdir", "$outputDir/images", $debug); }

$err = runcommand("$baseDir/riverforcing_hind.pl", "$baseDir $inputDir", $debug);
$err = runcommand("$baseDir/do_intrusionlength_hind.pl", "$baseDir $inputDir", $debug);

if (!-e "$outputDir/images/stats")
{ $err = runcommand("/bin/mkdir", "$outputDir/images/stats", $debug); }

# extract station model files into process directory

$err = runcommand("$baseDir/Extract/master_extract.pl", "$baseDir /home/workspace/ccalmr/hindcasts/reference/new_stations.sta $inputDir $outputDir/process", $debug);

# do skill assessment after extracting stations

$err = runcommand("$baseDir/SkillAssessment/do_sa_hind.pl", "$inputDir", $debug);

# generate station images

$err = runcommand("$baseDir/StationImages/master_stationImages.pl", "$baseDir $inputDir $runYear $runWeek 1", $debug);

# generate iso images

$err = runcommand("$baseDir/IsoImages/master_isoImages.pl", "$baseDir $inputDir 1", $debug); 

# generate cruise images

my $cdays  = '';
for (my $c =$startDay; $c<$endDay; $c++){
   $cdays.=" -cday $c";
}

my @outputpath;
my $outrun;
my $outrunyear;
my $outrunweek;
my $outrunv;
my $outrunvv;
@outputpath = split(/\//,$outputDir);
$outrun = pop(@outputpath);
my $outputbase = join("/",@outputpath);
($outrunyear, $outrunweek, $outrunv) = split(/\-/,$outrun);
$outrunvv = 'vv' . $outrunv;
$err = runcommand("$baseDir/CruiseImages/virtual_cruise_today_sigmaz.pl", "-rundir $outputbase $cdays -noadcp -notsg -h $outrunv -bin $baseDir/CruiseImages/", $debug);


# Temporary for dev. pturner 2005-04-15
#exit(0);

if ($official_flg == 1)
{
  my @outputpath;
  my $outrun;
  my $outrunyear;
  my $outrunweek;
  my $outrunv;
  my $outrunvv;

  @outputpath = split(/\//,$outputDir);
  $outrun = pop(@outputpath);
  ($outrunyear, $outrunweek, $outrunv) = split(/\-/,$outrun);
  $outrunvv = 'vv' . $outrunv;

#  chdir '/disk/ccalmr/hindcasts';
  chdir '/home/workspace/ccalmr/hindcasts';

  $err = runcommand("/bin/ln", "-s $outputDir $outrun", $debug);
  $err = runcommand("/bin/ln", "-s $outputDir/images images/$outrun", $debug);
  $err = runcommand("/bin/ln", "-s $outrun $outrunyear-$outrunweek-$outrunvv", $debug);

  if ($outrun =~  /\d{4}-\d{2}-(\w+)/ && $outrun =~ /\d$/)
  {
    chdir '/home/workspace/ccalmr/hindcasts/database';
    $err = runcommand("/bin/ln", "-s ../$outrun $outrun", $debug); 
  }
  else 
  {
    chdir '/home/workspace/ccalmr/hindcasts/calibration';
    $err = runcommand("/bin/ln", "-s ../$outrun $outrun", $debug);
    $err = runcommand("/bin/ln", "-s ../$outrunyear-$outrunweek-$outrunvv $outrunyear-$outrunweek-$outrunvv", $debug);
  }
}
