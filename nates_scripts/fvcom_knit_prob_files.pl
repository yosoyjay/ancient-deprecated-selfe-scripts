#!/usr/bin/perl
# 
# Knit together timeseries in fvcom output files
# 

use strict;
use CORIE;
use POSIX qw(ceil floor);

my $debug=0;
my ($i, $j);

if ($#ARGV!=4) {
    print "usage: $0 <input dir> <run_dir_prefix (eg day)> <start #> <end #> <output_dir>\n";
    exit 1;
}

my $indir = $ARGV[0];
my $rname = $ARGV[1];
my $start = $ARGV[2];
my $end = $ARGV[3];
my $out_dir = $ARGV[4];

my $file;

for ($i=$start; $i<=$end; ++$i) {
    print "$indir/$rname$i/output/timeseries/*.dat\n";

    chdir("$indir/$rname$i/output/timeseries/");

    my @dir = <*.dat>;
    print "$#dir\n";

    for ($j=0; $j<=$#dir; ++$j) {
	my @ts;
	my $k;
	if (-e "$out_dir/$dir[$j]") {
	    open(OUT, ">>$out_dir/$dir[$j]") or die "could not open $out_dir/$dir[$j] for writing\n";
	    @ts = read_fvcom_ts("$indir/$rname$i/output/timeseries/$dir[$j]");
	    for ($k=2; $k<=$#ts; ++$k) {
		print OUT $ts[$k];
	    }
	    close(OUT);
	}
	else {
	    open(OUT, ">$out_dir/$dir[$j]") or die "could not open $out_dir/$dir[$j] for writing\n";		
	    @ts = read_fvcom_ts("$indir/$rname$i/output/timeseries/$dir[$j]");
	    for ($k=0; $k<=$#ts; ++$k) {
		print OUT $ts[$k];
	    }
	    close(OUT);		
	}
    }
    
    closedir(DIR);
}


##############################################################################

sub read_fvcom_ts {
    my $fname = $_[0];
    open (TS, $fname) or die "could not open $fname for reading\n";

    my @lines = <TS>;
    $lines[0]=~/(\d+)/;
    my $num = $1;
    my @ret;
    $ret[0] = $num;

    $ret[1] = $lines[6];

    my $ii;
    for ($ii=9; $ii<=$#lines; ++$ii) {
	$ret[$#ret+1] = $lines[$ii];
    }

    return @ret;
}
