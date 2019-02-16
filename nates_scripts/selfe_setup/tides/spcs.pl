#!/usr/bin/perl -w

print "Enter input file name to be converted: ";
chop($filein = <STDIN>);

print "Enter output file name: ";
chop($fileout = <STDIN>);
 
print "What type of data is to be converted?\n  [1] long/lat\n  [2] state plane, nad27\n  [3] state plane, nad83\nEnter data type: ";
chop($intype = <STDIN>);

print "What type of output do you want?\n  [1] long/lat\n  [2] state plane, nad27\n  [3] state plane, nad83\nEnter data type: ";
chop($outtype = <STDIN>);

if ($intype == 1 || $outtype == 1) {
    print "Enter a digit [0-8] for the precision 'proj' uses for output: ";
    chop($prec = <STDIN>);
    die "Out of range precision...\n" unless ($prec < 9 && $prec >=0);
}

print "Which style is your input file?\n  [1] integer x-value y-value z-value\n          or\n  [2] integer x-value y-value\nEnter style number: ";
chop($style = <STDIN>);

print "Enter state plane coordinate system zone number: ";
chop($zone = <STDIN>);
die "Invalid spcs zone number...\n" unless ($zone =~ /\d+/);

print "Enter NADCON correction region: ";
chop($nadcor = <STDIN>);

print "What are the units for the cartesian system?\n  [1] feet\n  [2] meters\nEnter units type: ";
chop($units = <STDIN>);
if ($units == 1) {
    $unit = "us-ft";
} elsif ($units == 2) {
    $unit = "m";
} else {
    die "Invalid unit number...\n";
}


pipe(PROJR,PROJW);
if (fork) {
    close(PROJR);
    open(INFILE, "<" . $filein) ||
        die "Cannot open $filein for input\n";
#    $numctr=1;
    while (<INFILE>) {
        if ($style == 1) {
            ($num,$x,$y,$z) = split(/\s+/,$_,4);
            chop($z);
            print PROJW "$x $y $z\n";
        } elsif ($style == 2) {
            ($num,$x,$y) = split(/\s+/,$_,3);
            chop($y);
            print PROJW "$x $y\n";
        } else {
            die "Invalid style number...\n";
        }
    }
    exit;
} else {
    close(PROJW);
    open(STDIN, "<&PROJR");

    $errstring = "Not an allowed conversion...\nFor input # to ouput #, you may \ngo from 1 to 2,3\n   from 2 to 1,3\n    or from 3 to 1,2\n";

    if ($intype == 1) {
        if ($outtype == 2) {
            open(PROJ, "proj +init=nad27:$zone +units=$unit -W${prec}|") ||
                die "Cannot open pipe to proj...\n";
        } elsif ($outtype == 3) {
            open(PROJ, "proj +init=nad83:$zone +units=$unit -W${prec}|") ||
                die "Cannot open pipe to proj...\n";
        } else {
            die $errstring;
        }
    } elsif ($intype == 2) {
        if ($outtype == 1) {
            open(PROJ, "proj -I +init=nad27:$zone +units=$unit -W${prec}|") ||
                die "Cannot open pipe to proj...\n";
        } elsif ($outtype == 3) {
            if ($unit eq "us-ft") {
                open(PROJ, "nad2nad -i 27,spcs=$zone,feet -o 83,spcs=$zone,feet -r $nadcor|") ||
                    die "Cannot open pipe to proj...\n";
            } elsif ($unit eq "m") {
                open(PROJ, "nad2nad -i 27,spcs=$zone -o 83,spcs=$zone -r $nadcor|") ||
                    die "Cannot open pipe to proj...\n";
            }
        } else {
            die $errstring;
        }
    } elsif ($intype == 3) {
        if ($outtype == 1) {
            open(PROJ, "proj -I +init=nad83:$zone +units=$unit -W${prec}|") ||
                die "Cannot open pipe to proj...\n";
        } elsif ($outtype == 2) {
            if ($unit eq "us-ft") {
                open(PROJ, "nad2nad -i 83,spcs=$zone,feet -o 27,spcs=$zone,feet -r $nadcor|") ||
                    die "Cannot open pipe to proj...\n";
            } elsif ($unit eq "m") {
                open(PROJ, "nad2nad -i 83,spcs=$zone -o 27,spcs=$zone -r $nadcor|") ||
                    die "Cannot open pipe to proj...\n";
            }
        } else {
            die $errstring;
        }
    } else {
        die "Invalid input file type number...\n";
    }

    open(OUTFILE, ">" . $fileout) ||
        die "Cannot open $fileout for output\n";
    $ctr=0;
    $newctr=0;
    while(<PROJ>) {
        $newctr++;
        if ($outtype == 1) {
            if ($style == 1) {
                ($longt,$latt,$depth) = ($_ =~ /([\w\.'"]+)\s+([\w\.'"]+)\s+([\-\d\.]+)/);
            } else {
                ($longt,$latt) = ($_ =~ /([\w\.'"]+)\s+([\w\.'"]+)/);
            }
            ($longdeg,$longmin,$longsec,$longsgn) = ($longt =~ /(\d+)d(\d+)\'([\d.]+)\"(\w)/);
            ($latdeg,$latmin,$latsec,$latsgn) = ($latt =~ /(\d+)d(\d+)\'([\d.]+)\"(\w)/);
            $newlong = $longdeg+($longmin+($longsec/60))/60;
            $newlat = $latdeg+($latmin+($latsec/60))/60;
            if ($longsgn eq "W") {
                $newlong *= (-1);
            }
            if ($latsgn eq "S") {
                $newlat *= (-1);
            }
            if ($style == 1) {
                if ($newctr%1 == 0) {
                    printf(OUTFILE "%d %.${prec}f %.${prec}f %.3f\n", $ctr+1, $newlong, $newlat, $depth);
                }
            } else {
                if ($newctr%1 == 0) {
                    printf(OUTFILE "%d %.${prec}f %.${prec}f\n", $ctr+1, $newlong, $newlat);
                }
            }
        } else {
            if ($style == 1) {
                ($xout,$yout,$zout) = split(/\s+/,$_,3);
                chop($zout);
                if ($newctr%1 == 0) {
                    printf(OUTFILE "%d %.8f %.8f %.3f\n", $ctr+1, $xout, $yout, $zout);
                }
            } else {
                ($xout,$yout) = split(/\s+/,$_,2);
                chop($yout);
                if ($newctr%1 == 0) {
                    printf(OUTFILE "%d %.8f %.8f\n", $ctr+1, $xout, $yout);
                }
            }
        }
        $ctr++;
    }
}
print "Done\n";

