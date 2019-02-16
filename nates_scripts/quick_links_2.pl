#!/usr/bin/perl

use strict;
use POSIX qw(ceil floor);

if ($#ARGV!=1) {
    print "$0 <home dir> <target dir>\n";
    exit 1;
}

my $home_dir = $ARGV[0];
my $target_dir = $ARGV[1];


opendir(DIR,"$home_dir") or die "could not open directory $home_dir\n";
my @days = readdir(DIR) or die "unable to read $home_dir\n";
closedir(DIR);
my $day;

foreach $day (@days) {
    next if ($day eq ".");
    next if ($day eq "..");

    if (-d "$home_dir/$day") {
	if (!-e "$target_dir/$day") {
	    mkdir("$target_dir/$day");
	}	    
	if (-d "$home_dir/$day/process") {
	    opendir(DIR, "$home_dir/$day/process") or die "could not open directory $home_dir/$day/process\n";
	    my @names = readdir(DIR) or die "unable to read $home_dir/$day/process\n";
	    closedir(DIR);
	    my $name;
	    
	    if (!-e "$target_dir/$day/process") {
		mkdir("$target_dir/$day/process");
	    }	    

	    foreach $name (@names) {
		next if ($name eq ".");
		next if ($name eq "..");
		if (-l "$target_dir/$day/process/$name") {
		    runcommand("rm", "$target_dir/$day/process/$name", 0);		    
		}
		if (-f "$target_dir/$day/process/$name") {
		    runcommand("mv", "$target_dir/$day/process/$name $target_dir/$day/process/$name".".bkup", 0);
		}
		runcommand("ln", "-s $home_dir/$day/process/$name $target_dir/$day/process/$name", 0);		
	    }
	}
    }
}


########################################################

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
