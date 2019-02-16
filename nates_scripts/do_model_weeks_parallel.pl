#!/usr/bin/perl

use strict;

my $debug=0;

if ($#ARGV!=6) {
    print "$0 <base_dir> <run name> <year> <julian start day> <julian end day> <hotstart prefix> <nprocs>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0]."/".$ARGV[1]."/";
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
#my $ex = $ARGV[5];
my $hs_pre = $ARGV[5];
my $nprocs = $ARGV[6];
my $qsub = "qsub -pe orte $nprocs run.qsub";

my $lines = "#!/bin/bash\n#\$ -cwd\n#\$ -j y\n#\$ -S /bin/bash\n";
my $final_line = "/usr/mpi/gcc/openmpi-1.2.5/bin/mpirun --mca btl_tcp_if_include eth0 -np \$NSLOTS ";

my @dirs;

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {

    #### need to copy all the sflux links locally ####
    chdir("$home_dir$dirs[$i]/run/sflux/");
    opendir(DIR, ".");
    my $file;
    while ( $file = readdir(DIR)) {
	next if $file eq "." or $file eq "..";
	system("cp", "$file", "temp.tmp");
	system("mv temp.tmp $file");
    }
    closedir(DIR);

    #### make sure the outputs directory exists ####
    print "$home_dir$dirs[$i]/run/outputs\n";
    if (! -e "$home_dir$dirs[$i]/run/outputs") {
	mkdir("$home_dir$dirs[$i]/run/outputs");
    }

    chdir("$home_dir$dirs[$i]/run");
    open(FOUT, ">run.qsub") or die "could not open run.qsub for writing\n";
    print FOUT "$lines";
    print FOUT "$final_line$home_dir$dirs[$i]/run/pelfe\n";
    close(FOUT);

    @command=("$qsub");
    print "\n$qsub from $home_dir$dirs[$i]/run\n";
    system(@command);

    chdir("$home_dir$dirs[$i]/run/outputs/");    
    while (! -e "$hs_pre"."_0000_hotstart") {
        print "going to sleep...\n";
	sleep(1800);
    }
    sleep(60); # make sure the executable is done writing
    
    chdir("$home_dir$dirs[$i]/run/outputs/");
    system(("perl", "/home/users/hyde/scripts/combine_output.pl", "1", "1", "$hs_pre", "$nprocs"));
    system(("cp", "1_elev.61", "1_zcor.63", "1_temp.63", "1_salt.63", "1_hvel.64", ".."));
    system(("cp", "hotstart.in", "../$hs_pre"."_hotstart"));

    chdir("$home_dir$dirs[$i]/run/sflux/");    
    system(("rm", "*.nc"));
 
    chdir("$home_dir"); 

    if ($i<$#dirs) { 
	if (-e "$home_dir$dirs[$i+1]/run/hotstart.in") {
	    unlink("$home_dir$dirs[$i+1]/run/hotstart.in");
	}

	@command=("ln", "-sf", "$home_dir$dirs[$i]"."/run/$hs_pre"."_hotstart", "$home_dir$dirs[$i+1]/run/hotstart.in");
	print "@command\n";
	system(@command);
    }
}
