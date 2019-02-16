#!/usr/bin/perl

use Cwd;
use CMOP::CMOP;

######## set command line options ######################################

# presently set up so that the minimum number of arguments is 7
# (see below)

if ( @ARGV < 7) {
  my $buf = <<EOF;
  This utility should be executed from within the run directory.
  
  It creates the sflux subdirectory, as well as the links to data
  and the sflux_inputs.txt file.
  
  you must supply the following command line arguments:
  
  (1) data_source_code
  (2) year_start
  (3) month_start
  (4) day_start
  (5) year_stop
  (6) month_stop
  (7) day_stop
  (8) utc_start  (optional argument)
  (9) start_hour (optional argument)
  
  Example:
  % make_sflux_links.pl 1 2007 1 1 2007 1 7
  
  Numeric values for data_source_code:
  1 - NAM
  2 - GFS
  3 - NARR
  4 - NNRP
EOF
  print $buf;
  exit(1);
}

my $data_source_code = $ARGV[0];
my $year_start       = $ARGV[1];
my $month_start      = $ARGV[2];
my $day_start        = $ARGV[3];
my $year_stop        = $ARGV[4];
my $month_stop       = $ARGV[5];
my $day_stop         = $ARGV[6];

# set optional arguments
my $utc_start = '8.0';
my $start_hour = '0.0';

if ( @ARGV == 8 ) {
  $utc_start  = $ARGV[7];
}

if ( @ARGV == 9 ) {
  my $start_hour = $ARGV[8];
}

# output options to screen
print " data_source_code = $data_source_code\n";
print " year_start       = $year_start\n";
print " month_start      = $month_start\n";
print " day_start        = $day_start\n";
print " year_stop        = $year_stop\n";
print " month_stop       = $month_stop\n";
print " day_stop         = $day_stop\n";
print " utc_start        = $utc_start\n";
print " start_hour       = $start_hour\n";

if ($year_start < 2004) {
    $data_source_code = 3; # use narr for early stuff.
}

######## calculate starting and stopping Julian dates ##################

print "Note: including additional day at end (PST/UTC difference)\n";

my $cjd_start = getcorieday(0, 0, 0, $day_start, $month_start, $year_start);
my $cjd_stop  = getcorieday(0, 0, 0, $day_stop, $month_stop, $year_stop) + 1;

print "cjd_start = $cjd_start\n";
print "cjd_stop  = $cjd_stop\n";

######## specify the models to be used as data sources, locations, etc #

my $src_1 = "";
my $src_2 = "";
my $domain_1 = ".";
my $domain_2 = ".";
my $data_grandparent = '/home/workspace/ccalmr37/sflux_data/';

$src_1 = "nam";
$domain_1 = ".local.";

if ($data_source_code == 1) {
  $src_1 = "nam";
  $domain_1 = ".local.";
} elsif ($data_source_code == 2) {
  $src_1 = "gfs";

  my $cjd_switch    = getcorieday(0, 0, 0, 20, 10, 2005);
  my $cjd_switch_p1 = getcorieday(0, 0, 0, 21, 10, 2005);

  if ($cjd_start < $cjd_switch) {
    if ($cjd_stop <= $cjd_switch) {
      $domain_1 = ".local.";
    } elsif ($cjd_stop == $cjd_switch_p1) {
      $domain_1 = ".local.";
      print "Removing additional day, due to GFS local/global switch!\n";
      $cjd_stop = $cjd_switch
    } else {
      print "echo Run extends across GFS local/global switch!\n";
      exit(1);
    }
  }
} elsif ($data_source_code == 3) {
  $src_1 = "narr";
  $domain_1 = ".";

} elsif ($data_source_code == 4) {
  $src_1 = "nnrp"

} else {
  print "undefined data_source_code. . .\n";
  exit(1);
}

print "Model Sources\n";
print "src_1 = $src_1\n";
print "src_2 = $src_2\n";


######## create the destination directory (replace if it exists) #######

# create name of destination directory
my $local_dir = getcwd();
my $dest_dir  = "$local_dir/sflux/";
my $old_dest_dir  = "$local_dir/sflux.old/";

# fail if directory already exists
if ( -e $dest_dir ) {
  print "destination directory already exists!\n";
  print "renaming it!\n";
  system("/bin/rm -rf $old_dest_dir");
  system("/bin/mv $dest_dir $old_dest_dir");
}

# create destination directory
mkdir ($dest_dir);

print "created destination directory: $dest_dir\n";

######## create the sflux_inputs file ##################################
my $inputs_file = 'sflux_inputs.txt';

print "Creating sflux_inputs file at:\n";
print "  $dest_dir/$inputs_file\n";

my $buf = <<EOF;
&sflux_inputs
start_year  = $year_start,
start_month = $month_start,
start_day   = $day_start,
start_hour  = $start_hour,
utc_start   = $utc_start,
/
EOF

open(OUTP, ">$dest_dir/$inputs_file") or die "Unable to open $dest_dir/$inputs_file\n";;
print OUTP $buf;
close(OUTP);

######## loop over file types and dates, creating links ################

# loop over possible file types
#foreach my $file_type ("air", "rad", "prc") {
foreach my $file_type ("air", "rad", "prc") {
  print "-----------------------------------------------------------------\n";
  print "Creating links for file_type: $file_type\n";

# initialize counter for files of this type
  my $type_counter = 0;

# top of loop over Julian days
  my $cjd = $cjd_start;
  for ($cjd = $cjd_start; $cjd<=$cjd_stop; $cjd++) {
    $cjd_switch_p1 = getdatefromcorieday($cjd);
    my @day = getdatefromcorieday($cjd);
    my $day   = $day[3];
    my $month = $day[4];
    my $year = $day[5];
    my $date_string = sprintf("%4d_%02d_%02d", $year, $month, $day);
#    print "echo Julian date and date_string : $cjd $date_string [@day]\n";

# determine where data will come from for this day (and go there)
    my $data_parent = "$data_grandparent$src_1/";
    my $data_dir    = sprintf("$data_parent%4d_%02d/", $year, $month);
    chdir($dest_dir);

# loop over files that satisfy required name (this skips missing files)
    my $lsarg = "$src_1"."_$file_type$domain_1$date_string.nc";
#    print "[$data_dir] [$lsarg]\n";
    foreach my $in_file (`ls $data_dir/$lsarg`) {
      chomp($in_file);
      $type_counter += 1;
      $type_counter_label = sprintf("%03d", $type_counter);
      my $link_name = "sflux_".$file_type."_1.".$type_counter_label.".nc";
      #print "creating link from:\n";
      #print "$data_dir/$in_file\n";
      #print "to:\n";
      #print "$dest_dir/$link_name\n";
#      print ("/bin/ln -s $in_file $dest_dir/$link_name\n");
# rather than symboloc link, copy the file - a sirus cluster not NFS mounts on nodes problem.
      system("/bin/cp -p $in_file $dest_dir/$link_name");
      #system("/bin/ln -s $in_file $dest_dir/$link_name");
      #system("/bin/ln -s $data_dir/$in_file $dest_dir/$link_name");
# end of foreach loop that selects files
    }
  }
# end of loop over file types
}

