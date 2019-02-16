#!/usr/bin/perl

use CORIE;

if (@ARGV < 3) {
    print "usage $0 <base_dir> <file_base> <out_file>\n";
}

$input_dir = $ARGV[0];
$file_base = $ARGV[1];
$fout = $ARGV[2];


open(OUT, ">$fout") or die "could not open $fout for writing\n";;
opendir(DIR, $input_dir);

my @mtimes = ();
my %fnames = {};

while($file = readdir(DIR)) {
    if ($file =~ /$file_base/) {
	@fstat = stat($input_dir."/".$file);
#	print $file." ".$fstat[7]." \n";
	$mtimes[$#mtimes+1] = $fstat[9];
	$fnames{$fstat[9]} = $file;
    }
}

#print @mtimes;

@times = sort {$a <=> $b} @mtimes;

for ($i=1; $i<=$#times; ++$i) {
    print OUT $fnames{$times[$i]}." - ".$fnames{$times[$i-1]}."   ".(($times[$i] - $times[$i-1]) / 60 )." mins \n";
#    print OUT (($times[$i] - $times[$i-1]) / 3600 )."\n";
}

print OUT $fnames{$times[$#times]}." - ".$fnames{$times[0]}."   ".(($times[$#times] - $times[0]) / 3600 )." hours\n";

print "# hours per model day == ".(scalar(@times) / (($times[$#times] - $times[0]) / 3600 ))."\n";

close(OUT);
