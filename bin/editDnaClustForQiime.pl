#! /usr/bin/perl -w 

#10/05/13	Fiona J. Whelan
#This program will format the output from dnaclust into rep_set.fna and seqs_otus.txt for use with qiime.
#
#Usage: editDnaClustForQiime.pl dnaclustoutputfile

#Check usage
if (($#ARGV+1) != 3) {
	print "\nUsage: editDnaClustForQiime.pl seqs.fna dnaclustoutputfile gg\n";
	print "Will produce rep_set.fna and seqs_gg_otus.txt\n";
	exit;
}
$seqs = $ARGV[0];
$inf = $ARGV[1];
$gg = $ARGV[2];
open (IN, "<", $inf) or die "$!";
open (SEQS, "<", $seqs) or die "$!";
open (SEQS_OTUS, ">", "seqs_".$gg."_otus.txt") or die "$!";
open (REP_SET, ">", "rep_set.fna") or die "$!"; 
#load each sequence into memory in the form of a hash
%reference = ();
@seqlines = <SEQS>;
$a = 0;
$woline = $seqlines[$a];
$a++;
while($a < $#seqlines) {
	$woline =~/>(.*?)\s.*/;
	$key = $1;
	$value = $seqlines[$a];
	$a++;
	$reference{ $key } = $value;
	$woline = $seqlines[$a];
	$a++;
}
close SEQS;
#go through IN, create rep_set and seqs_otus
@lines = <IN>;
$i = 0;
$woline = $lines[$i];
while($i < $#lines+1) {
	$woline =~/^(.*?)\t/;
	$rep = $1;
	$woline =~/(.*)\t\n/;
	print SEQS_OTUS $i."\t".$1."\n";
	print REP_SET ">".$i." ".$rep."\n";
	print REP_SET $reference{$rep};
	$i++;
	$woline = $lines[$i];
}
close IN;
close SEQS_OTUS;
close REP_SET;
