#! /usr/bin/perl -w 

#26/06/15	Fiona J. Whelan
#This program will format the output from UPARSE, results.txt, into seqs_otus.txt for use with qiime.
#
#Usage: editUPARSEForQiime.pl results.txt

#Check usage
if (($#ARGV+1) != 2) {
	print "\nUsage: editUPARSEForQiime.pl results.txt gg\n";
	print "Will produce seqs_<gg>_otus.txt\n";
	exit;
}
$seqs = $ARGV[0];
$gg = $ARGV[1];
open (SEQS, "<", $seqs) or die "$!";
open (SEQS_OTUS, ">", "seqs_".$gg."_otus.txt") or die "$!";
#load each sequence into memory in the form of a hash
%reference = ();
@seqlines = <SEQS>;
$a = 0;
$woline = $seqlines[$a];
$a++;
while($a < $#seqlines) {
	$woline =~/(.*?);size=.*?;\t(.*?)\t.*?\t\*\t(OTU_\d*)/;
	$key = $3;
	$value = $1;
	if ($2 ne "chimera") {
		if (exists $reference{$key}) {
			$existingValue = $reference{$key};
			$existingValue = $existingValue."\t".$value;
			$reference{ $key } = $existingValue;
		} else {
			$reference{ $key } = $value;
		}
	}
	$a++;
	$woline = $seqlines[$a];
}
close SEQS;
#go through all $keys in numerical order until one doesn't exist
$a=1;
$key = "OTU_".$a;
while(exists $reference{$key}) {
	print SEQS_OTUS $key."\t".$reference{$key}."\n";
	$a++;
	$key = "OTU_".$a;
}
close SEQS_OTUS;
