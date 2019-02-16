#!/usr/bin/perl

use strict;
use CORIE;
use POSIX qw(ceil floor);

my $debug=0;

if ($#ARGV!=6) {
    print "usage: $0 <base_dir> <run name> <year> <julian start day> <julian end day> <ref dir> <time step>\n";
    exit 1;
}

my @command;
my $home_dir = $ARGV[0]."/".$ARGV[1]."/";
my $year = $ARGV[2];
my $sday = $ARGV[3];
my $eday = $ARGV[4];
my $ref_dir = $ARGV[5];
my $tstep = $ARGV[6];

my @dirs;
my ($m1, $m2, $d1, $d2);

my $i;
for ($i=$sday; $i<=$eday; ++$i) {
    $dirs[$#dirs+1]=sprintf("%d-%03d", $year, $i);
    print "$dirs[$#dirs]\n";
}

for ($i=0; $i<=$#dirs; ++$i) {
    
    if (!-e "$home_dir$dirs[$i]") {
	mkdir("$home_dir$dirs[$i]");
    }
    if (!-e "$home_dir$dirs[$i]/run") {
	mkdir("$home_dir$dirs[$i]/run");
    }

    ($m1, $d1) = getmonthday($year, $sday+$i);
    ($m2, $d2) = getmonthday($year, $sday+$i+1);

    chdir("$home_dir$dirs[$i]/run");
#    runcommand("/home/workspace/ccalmr/mazulauf/amb10xx/netcdf/cvs_stuff/forecasts/bin/atmos_nc/scripts/make_sflux_links.csh","3 $year $m1 $d1 $year $m2 $d2", 0);
    
    runcommand("cp", "$ref_dir/copy/* $home_dir$dirs[$i]/run/", 0);

#    runcommand("cp", "$ref_dir/copy/tid_e $home_dir$dirs[$i]/run/", 0);
#    runcommand("cp", "$ref_dir/copy/fort.8 $home_dir$dirs[$i]/run/", 0);
#    setup_paramin("$ref_dir/copy/param.in", "$home_dir$dirs[$i]/run/", "$ref_dir/setup_paramin_args.txt", $year, $m1, $d1, $tstep, "42");
#    setup_nudge("$home_dir/nudge", "$home_dir$dirs[$i]/run/", "$home_dir/temp/", $year, $m1, $d1);
}



##############################################

sub setup_paramin() {
    my $pfile = $_[0];
    my $out_dir = $_[1];
    my $args_file = $_[2];
    my $y = $_[3];
    my $m = $_[4];
    my $d = $_[5];
    my $ts = $_[6];
    my $nbnds = $_[7];
    my $jj;

    my $mth_txt = " 4 4 4 4 ! mouth";
    # constants for tide.out values:
 #   my $tO1_line = 12;
 #   my $tK1_line = 20;
 #   my $tQ1_line = 10;
 #   my $tP1_line = 18;
 #   my $tK2_line = 42;
 #   my $tN2_line = 31;
 #   my $tM2_line = 35;
 #   my $tS2_line = 40;
    
    # constants for param.in values.  Adjust as appropriate (primarily output_line, based on # of open boundaries and open boundary nodes)
    ######### for reference only.  changed to read setup_paramin_args.txt for these values
    #    my $datetime_line = 1;
    #    my $dt_line       = 15;        # 13; 
    #    my $nws_line      = 23;        # 20; 
    #    my $output_line   = 619;     #39;     # 618;    
    #    my $O1_line       = 32;      #-1;      # 29; 
    #    my $K1_line       = 34;      #-1;      # 31; 
    #    my $Q1_line       = 36;      #-1;      # 33; 
    #    my $P1_line       = 38;      #-1;      # 35; 
    #    my $K2_line       = 40;      #-1;      # 37; 
    #    my $N2_line       = 42;      #-1;      # 39; 
    #    my $M2_line       = 44;      #-1;      # 41; 
    #    my $S2_line       = 46;      #-1;      # 43; 

    open(PARAMINARGS, "$args_file") or die "could not open $args_file\n";
    my @input = <PARAMINARGS>;
    close(PARAMINARGS);
    my $datetime_line = (split(/\t/, $input[0]))[0];
    my $dt_line       = (split(/\t/, $input[1]))[0];
    my $nws_line      = (split(/\t/, $input[2]))[0];
    my $output_line   = (split(/\t/, $input[3]))[0];
    my $O1_line       = (split(/\t/, $input[4]))[0];
    my $K1_line       = (split(/\t/, $input[5]))[0];
    my $Q1_line       = (split(/\t/, $input[6]))[0];
    my $P1_line       = (split(/\t/, $input[7]))[0];
    my $K2_line       = (split(/\t/, $input[8]))[0];
    my $N2_line       = (split(/\t/, $input[9]))[0];
    my $M2_line       = (split(/\t/, $input[10]))[0];
    my $S2_line       = (split(/\t/, $input[11]))[0];
    my $boundary_line = (split(/\t/, $input[12]))[0];

    my $tcoluline1;
    my $tcoluline2 = "     8615  SVI tidal model      GMT 461512400 !Do not change this line";

    chdir("$out_dir");

    open(PARAMIN, "$pfile") or die "could not open $pfile\n";
    open(PARAMOUT, ">$out_dir/param.in") or die "could not open $out_dir/param.in for writing\n";
    my @vals;

#    open(TIDECOLU, ">tide_colu.com") or die "could not open $out_dir/tide_colu.com for writing\n";
#
#    my $yrthhu = floor($y/100);
#    my $yrteon = $y - $yrthhu*100;
#    
#    $tcoluline1 = sprintf("  08%02d%02d%02d%02d", $d, $m, $yrteon, $yrthhu);
#    print TIDECOLU "$tcoluline1\n$tcoluline2\n";
#    close(TIDECOLU);
#    
#    runcommand("./tid_e", "< tide_colu.com > tide.out", 0);
#    
#    open(TIDEOUT, "$out_dir/tide.out") or die "could not open $out_dir/tide.out\n";
#    my @tout = <TIDEOUT>;
#    
#    @{$vals[0]} = split(/\s+/, $tout[$tO1_line]);
#    @{$vals[1]} = split(/\s+/, $tout[$tK1_line]);
#    @{$vals[2]} = split(/\s+/, $tout[$tQ1_line]);
#    @{$vals[3]} = split(/\s+/, $tout[$tP1_line]);
#    @{$vals[4]} = split(/\s+/, $tout[$tK2_line]);
#    @{$vals[5]} = split(/\s+/, $tout[$tN2_line]);
#    @{$vals[6]} = split(/\s+/, $tout[$tM2_line]);
#    @{$vals[7]} = split(/\s+/, $tout[$tS2_line]);
    
    my @pin = <PARAMIN>;

    for ($jj=0; $jj<=$#pin; ++$jj) {
	if ($jj==$dt_line) {
	    print PARAMOUT "$ts   dt\n";
	}
	elsif ($jj==$datetime_line) {
	    my $temp=sprintf("%02d/%02d/%d 00:00:00 PST", $m, $d, $y);
#	    print "$temp\n";
	    print PARAMOUT "$temp\n";
	}
	elsif ($jj==$nws_line) {
	    my @vars = split(/\s+/, $pin[$jj]);
	    print PARAMOUT "$vars[0] $ts    nws\n";
	}
	elsif ($jj==$output_line) {
	    print PARAMOUT (900/$ts)." ".(900/$ts*96)."\n";
	}
	elsif ($jj==$boundary_line) {
	    print PARAMOUT "$nbnds$mth_txt\n";
	}
	elsif ($jj==$O1_line) {
	    my @tmp = split(/\s+/, $pin[$O1_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[0][2] $vals[0][3]\n";
	}
	elsif ($jj==$K1_line) {
	    my @tmp = split(/\s+/, $pin[$K1_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[1][2] $vals[1][3]\n";
	}
	elsif ($jj==$Q1_line) {
	    my @tmp = split(/\s+/, $pin[$Q1_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[2][2] $vals[2][3]\n";
	}
	elsif ($jj==$P1_line) {
	    my @tmp = split(/\s+/, $pin[$P1_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[3][2] $vals[3][3]\n";
	}
	elsif ($jj==$K2_line) {
	    my @tmp = split(/\s+/, $pin[$K2_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[4][2] $vals[4][3]\n";
	}
	elsif ($jj==$N2_line) {
	    my @tmp = split(/\s+/, $pin[$N2_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[5][2] $vals[5][3]\n";
	}
	elsif ($jj==$M2_line) {
	    my @tmp = split(/\s+/, $pin[$M2_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[6][2] $vals[6][3]\n";
	}
	elsif ($jj==$S2_line) {
	    my @tmp = split(/\s+/, $pin[$S2_line]);
	    if (scalar(@tmp) > 3) {
		$tmp[0] = $tmp[1];
	    }
	    print PARAMOUT " $tmp[0] $vals[7][2] $vals[7][3]\n";
	}
	else {
	    print PARAMOUT $pin[$jj];
	}	    
    } 
    close(PARAMOUT);
}


# target directory must have hgrid.ll hgrid.gr3 vgrid.in estuary.gr3 already in place
# base dir must have the exectuable (hard-coded in the function)

sub setup_nudge() {
    my $base_dir = $_[0];
    my $out_dir = $_[1];
    my $temp_dir = $_[2];
    my $y = $_[3];
    my $m = $_[4];
    my $d = $_[5];

    my $ex = "readncom7e_nate";   
    my $datein_l2 = "0";
    my $datein_l3 = "SELFE";
    my $datein_l4 = "10. 33.5";
    my $datein_l5 = "0";

    runcommand("cp", "$out_dir/hgrid.* $out_dir/vgrid.in $out_dir/estuary.gr3 $temp_dir", 0);
    runcommand("cp", "$base_dir/$ex $temp_dir", 0);
    runcommand("ln", "-s /usr/local/netcdf/include/TYPESIZES.mod $temp_dir/TYPESIZES.mod", 0);
    runcommand("ln", "-s /usr/local/netcdf/include/NETCDF.mod $temp_dir/NETCDF.mod", 0);

    chdir("$temp_dir");

    open(DATEIN, ">date.in") or die "could not open $temp_dir/date.in for writing\n";

    print DATEIN "$y $m $d 1\n";
    print DATEIN "$datein_l2\n";
    print DATEIN "$datein_l3\n";
    print DATEIN "$datein_l4\n";
    print DATEIN "$datein_l5\n";

    runcommand("$ex", "", 0);
    
    runcommand("cp", "salt_nu.in temp_nu.in $out_dir", 0);
}
