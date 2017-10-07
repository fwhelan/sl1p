#! /usr/bin/perl -w

#08/05/13	Fiona J. Whelan
#This progam will remove all
#	Root
#	Root;Other
#	Unassignable
#from an inputted otu table to produce a version lacking all unassigned.
#
#Usage: removeRoot.pl infile outfile\

#Check usage
if (($#ARGV+1) != 2) {
	print "\nUsage: removeRoot.pl infile outfile\nDo not include file endings: script assumes you have a infile.biom and infile.txt\n";
	exit;
}
$inf = $ARGV[0];
$otf = $ARGV[1];
open (IN, "<", $inf.".txt") or die "$!";
open (TMP, ">", "noRoot_filtered_otus.txt") or die "$!";
@lines = <IN>;
$i = 0;
$woline = $lines[$i];
$i++;
while($i <= $#lines+1) {
	if ($woline=~/Root$/ or $woline=~/Unassignable$/ or $woline=~/Root;Other$/ or $woline=~/No blast hit$/ or $woline=~/Unclassified/) {
		$woline=~/^(.*?)\t/;
		print TMP $1."\n";
		#print $1."\n";

	}
	$woline = $lines[$i];
	$i++;
}
close IN;
close TMP;
#replace any spaces in the pwd path with \/s
$inf=~s/\s/\\ /g;
$otf=~s/\s/\\ /g;
`filter_otus_from_otu_table.py -i $inf.biom -o $otf.biom -e noRoot_filtered_otus.txt`;
`rm noRoot_filtered_otus.txt`;
