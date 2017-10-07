#!/usr/bin/perl -w
#This script will reformat the rep_set.fna.clust output file from running
#	<<AbundantOTU -i ../splits_dir/seqs.fna -o rep_set.fna>>
#into a file that qiime can use as input into make_otu_table.py (namely, of the same format as seq_otus.txt
#
#Usage: abundantOTU_to_qiime.py infile outfile

#Check usage
if (($#ARGV+1) != 2) {
	print "\nUsage: abundantOTU_to_qiime.py infile outfile\n";
	exit;
}
$inf = $ARGV[0];
$otf = $ARGV[1];
open(IN, "<", $inf) or die "$!";
open(OUT, ">", $otf) or die "$!";
@lines = <IN>;
$i = 0;
$woline = $lines[$i];
$i++;
while($i < $#lines){
	#Each new consensus cluster is marked with a line beginning with an #Consensus x
	if($woline!~/^#Consensus/) {
		print "Error, finding header line on line:\n";
		print $woline;
		exit;
	}
	$woline=~/^#Consensus\s(\d+).*/;
	$cons = $1;
	print OUT $cons;
	print $cons."\n";
	$woline=$lines[$i];
	$i++;	
	do{
		#$woline=~/\s+seq\s+\d+\s+(\w+?_\d+)\s.*/;
		$woline=~/.*?(\S+?_\d+)\s.*/;
		$seqid = $1;
		print OUT "\t".$seqid;
		$woline=$lines[$i];
		$i++;	
	}
	while($woline!~/^#Consensus/ && $i < $#lines);
	print OUT "\n";
}
close IN;
close OUT;
