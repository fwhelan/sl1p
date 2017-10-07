#! /usr/bin/perl -w 

#04/06/13	Fiona J. Whelan
#This program will format the output from CLUSTOM into rep_set.fna and seqs_otus.txt for use with qiime.
#
#Usage: editCLUSTOMForQiime.pl _ClusterFolder _repFile

#Check usage
if (($#ARGV+1) != 2) {
	print "\nUsage: editDnaClustForQiime.pl _ClusterFolder _repFile\n";
	print "Will produce rep_set.fna and seqs_otus.txt\n";
	exit;
}
$clustFolder = $ARGV[0];
$rep = $ARGV[1];
open (REP_IN, "<", $rep) or die "$!";
open (SEQS_OTUS, ">", "seqs_otus.txt") or die "$!";
open (REP_SET, ">", "rep_set.fna") or die "$!"; 

#go through _repFile; alter for rep_set.fna
@lines = <REP_IN>;
$z = 0;
$woline = $lines[$z];
while ($z < $#lines) {
	$woline=~/^\>(.*?)\tcluster=C(\d*)/;
	print $1." ".$2."\n";
	print REP_SET ">".$2." ".$1."\n";
	$z++;
	$woline = $lines[$z];
	$woline=~/(\w*)/;
	print REP_SET $1."\n";
	$z++;
	$woline = $lines[$z];
}
close REP_IN;
close REP_SET;

#go through cluster folder; format for seqs_otu.txt
chdir($clustFolder);
system("ls *.cluster > fofn_clusters");
open(FOFN, "<", "fofn_clusters") or die "$!\n";
@lines = <FOFN>;
$y = 0;
$woline = $lines[$y];
while ($y < $#lines+1) {
	$woline=~/(.*?)\.cluster/;
	print $1."\n";
	print SEQS_OTUS $1;
	chomp($woline);
	open(CLUST, "<", $woline) or die "$!\n";
	@clustfiles = <CLUST>;
	$c = 0;
	$woclustline = $clustfiles[$c];
	while ($c < $#clustfiles+1) {
		$woclustline=~/^\>(\w*).*\n/;
		print SEQS_OTUS "\t".$1;
		$c = $c + 2; #skip seq line, only need name
		$woclustline = $clustfiles[$c];
	}
	close CLUST;
	print SEQS_OTUS "\n";
	$y++;
	$woline = $lines[$y];
}
close FOFN;
close SEQS_OTUS;
