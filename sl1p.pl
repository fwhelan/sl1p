#! /usr/bin/perl -w
use Switch;
use strict; use warnings;
#Author: Fiona J Whelan <whelanfj@mcmaster.ca>
#Last Modified Date: April 4 2017
#License: GPL
my $sl1p_ver = "4.3";
my $bin = $ENV{'SL1P_BIN'};
my $db = $ENV{'SL1P_DB'};
if ((! defined $bin) && (! defined $db)) {
	print "SL1P_BIN and SL1P_DB ENV variables are not initialized. Please see INSTALL.txt for details.\n";
	exit;
} elsif (! defined $bin) {
	print "SL1P_BIN ENV variable is not initialized. Please see INSTALL.txt for details.\n";
	exit;
} elsif (! defined $db) {
	print "SL1P_DB ENV variable is not initialized. Please see INSTALL.txt for details.\n";
	exit;
}
#Will automate sequence "processing": raw sequences off the sequencer through to QIIME analysis after asking a series of questions as to the clustering and taxonomic assignment methods etc. the user would like to use. Defaults are suggested.
#
#Usage: sl1p.pl <# of runs> <fofn for each run> <projectname>
#		e.g. sl1p.pl 2 fofn1.txt fofn2.txt Fiona_Proj

my $usage =  "\nUsage: sl1p.pl <# of runs> <fofn for each run> <projectname> [-r region] [-s seq info file] [-b barcode location] [-q quality filter] [-k keep filtered output] [-d taxonomy database] [-p OTU picking algorithm] [-c clustering threshold] [-l linkage] [-m chimera checking] [-t taxonomic assignment method] [-x timed results] [-u CPU threads] [-f force overwrite] \n";

my $help = "\nWill automate sequence \"processing\": raw v3 sequences off the sequencer through to QIIME analysis after asking a series of questions as to the clustering and taxonomic assignment methods etc. the user would like to use. Defaults are suggested.\n\nAll .fastq files must be in the directory that you call the script from.  Multiple runs must be split into multiple fofn files.\n
	-r which region of the 16S rRNA gene sequenced; if other, -s must be defined [v3/v34/v4/other] (default: v3)
	-s a sequence information file including primers, and min and max sequence length (default: none)
	-b location of the sequence barcode [fwd/rev/mixed] (default: fwd)
	-q the quality threshold for filtering (default: 30)
	-k keep .fastq output from quality filtering process [y/n] (default:n)
	-d the greengenes database that you would like to use for taxonomic assignment [gg2011/gg2013/silva111/all] (default: gg2011)
	-p OTU picking algorithm [abundantotu/uclust/cdhit/dnaclust/uclust-ref/uclust-ref-strict/mothur/blast/uparse/all] (default: abundantotu)
	-c OTU picking clustering threshold [0-100] (default: 97)
	-l the linkage method to use with mothur [nearest/average/furthest/all] (default: average)
	-m whether the user would like to employ chimera checking [y/n] (default: y)
	-t the taxonomic assignment method [blast/rdp-training/all] (default: rdp-training)
	-e certainty cutoff for taxonomic assignment as designed by Qiime [0-1] (default: 0.5)
	-x do you want timed OTU picking and overall results [y/n\ (default: n)
	-u how many CPU threads you would like to envoke [0-9] (default: 1)
	-f force overwrite .fna, map, picked_otu, and split_dir files from previous runs [y/n] (default: n)\n";

#Check if help flag used
my @vars = splice @ARGV;
my $search_for = "-h";
my ($h_index)= grep { $vars[$_] eq $search_for } 0..$#vars;
if (defined($h_index)) {
	print "\nsl1p.pl v$sl1p_ver Fiona J Whelan <whelanfj\@mcmaster.ca>\n";
        print "$usage";
	print "$help";
        exit;
}

#Check usage
if (($#vars+1) < 3) {
        print "\nsl1p.pl v$sl1p_ver Fiona J Whelan <whelanfj\@mcmaster.ca>\n";
	print "\nError: Mismatch in number of input arguments. There must be minimum 3 arguments as shown below. \n";
        print $usage;
	print "\nUse sl1p.pl -h for more details\n";
        print "\nExiting....\n";
        exit;
}

#Save var input
my $numruns = $vars[0];
if ($numruns !~/^\d+$/) {
	print "The # of runs inputted ($numruns) does not appear to be a number.\n\n";
	print $usage;
	print "\nUse sl1p.pl -h for more details\n";
        print "\nExiting....\n";
        exit;
}
my @fofns = [];
for (my $i = 1; $i <= $numruns; $i++) {
        $fofns[$i-1] = $vars[$i];
        #ensure that file exists
        unless (-e $fofns[$i-1]) {
                print "Error: fofn file doesn't exist: $fofns[$i-1]\n$!\n";
        	print $usage;
		print "Use sl1p.pl -h for more details\n";
        	print "\nExiting....\n";
                exit;
        }
	#ensure that file has >=2 lines
	my $cmd = `wc -l $fofns[$i-1]`;
	$cmd=~/^(\d+).*$/; $cmd=$1;
	if ($cmd < 2) {
		print "Error: fofn file $fofns[$i-1] doesn't have >=2 lines. Paired-end reads are expected.\n";
		print $usage;
		print "\nPlease see the readme for more details\n";
		print "\nExiting...\n";
		exit;
	} 
}
#Check to be sure project name exists
my $proj;
if ($numruns+1 < $#vars+1) {
	$proj = $vars[$numruns+1];
	if ($proj=~/^-/) {
		print "Please provide a project name.\n";
  	        print $usage;
		print "Use sl1p.pl -h for more details\n";
  	    	print "\nExiting....\n";
  		exit;
	}
} else {
	print "Please provide a project name.\n";
	print $usage;
	print "Use sl1p.pl -h for more details\n";
	print "\nExiting....\n";
	exit;
}

#Setup log file
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
(my $date, my $min, my $hour, my $mday, my $mon) = localtime();
$date = $months[$mon]."".$mday."_".$hour."h".$min."m";
my $log = "log_sl1p_".$date.".log";
open LOG, ">", $log or die "$!\n";
#Include full path in $log so that I don't have to trace what sub-folder I'm in.
my $pwd = `pwd`;
chomp($pwd);
$log = $pwd."/".$log;
print LOG "sl1p V$sl1p_ver\n";
print LOG "---------------\n";
print LOG "<usage: sl1p.pl ". join(" ", @vars)."\n";#.join(" ", @opts).">\n";
#Setup err file
my $err = "log_sl1p_".$date.".err";
open ERR, ">", $err or die "$!\n";
#Include full path in $err so that I don't have to trace what sub-folder I'm in.
$err = $pwd."/".$err;
#Save full path to a variable for time so I don't have to trace what sub-folder I'm in.
my $time_pwd = $pwd;
#Setup seq info file
open SEQ, ">", "log_seq_info.txt" or die "$!\n";

###############################################################
#Check for user inputs
###############################################################
print LOG "###############################################################\n";
print LOG "#Check for user inputs\n";
print LOG "###############################################################\n";
my $region = "v3";
my $seqinfofile = "";
my $barloc = "fwd";
my $mult_bc = 0;
my $qual = 30;
my $keep = "n";
(my $fwd_primer, my $rev_primer, my $fwd_revcomp_primer, my $rev_revcomp_primer, my $overlap, my $splitlib_min, my $splitlib_max) = &getPrimers($barloc, $region);
my @prim; my @fwd_primer; my @rev_primer; my @fwd_revcomp_primer; my @rev_revcomp_primer;
my $gg = "gg2011";
my $clust = "abundantotu";
my $clustThres = "97";
my $linkage = "average";
my $taxon = "rdp-training";
my $taxon_cutoff = "0.5";
my $time = "n";
my $force = 0;
my $threads = 1;
my $chimeraChecking = "y";

if ($#vars+1 > $numruns+2) {
	my @opts = splice @vars, $numruns+2;
	#-r region
	$search_for = "-r";
	my ($r_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($r_index)) {
		if ($r_index+1 < $#opts+1) {
			$region = $opts[$r_index+1];
			if($region ne "v3" && $region ne "v34" && $region ne "v4" && $region ne "other") {
        	        	print "-r which region of the 16S rRNA gene sequenced [v3/v34i/v4/other] (default: v3)\n";
				print "Not a valid response, try again.\n";
				exit;
			} else {
				if ($region eq "v3" || $region eq "v34" || $region eq "v4") {
					($fwd_primer, $rev_primer, $fwd_revcomp_primer, $rev_revcomp_primer, $overlap, $splitlib_min, $splitlib_max) = &getPrimers($barloc, $region);
					print LOG "sequence barcode: $barloc\n";
				} else {
					#for -r other, -s must be provided
					$search_for = "-s";
					my ($s_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
					if (defined($s_index)) {
						if ($s_index+1 < $#opts+1) {
							$seqinfofile = $opts[$s_index+1];
							($fwd_primer, $rev_primer, $fwd_revcomp_primer, $rev_revcomp_primer, $overlap, $splitlib_min, $splitlib_max, $barloc) = &getPrimersFile ($seqinfofile);
						} else {
							print "You provided -s with no option.\n";
							print "-s a sequence information file including primers, and min/max sequence length (default: none)\n";
							print "Not a valid response, try again.\n";
							exit;

						}
					} else {
						print "You must specify -s when using -r other.\n";
						exit;
					} 
				}
			}
		} else {
			print "You provided -r with no option.\n";
			print "-r which region of the 16S rRNA gene sequenced [v3/v34] (default: v3)\n";
                        print "Not a valid response, try again.\n";
                        exit;
		}
	}
	#-b barcode location
	$search_for = "-b";
        my ($b_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($b_index)) {
		if ($b_index+1 < $#opts+1) {
			$barloc = $opts[$b_index+1];
			if ($barloc ne "fwd" && $barloc ne "rev" && $barloc ne "mixed") {
				print "-b location of the sequence barcode [fwd/rev/mixed] (default: fwd)\n";
        	        	print "Not a valid response, try again.\n";
				exit;
			}
        		if ($barloc eq "rev" || $barloc eq "fwd") {
        			($fwd_primer, $rev_primer, $fwd_revcomp_primer, $rev_revcomp_primer, $overlap, $splitlib_min, $splitlib_max) = &getPrimers($barloc, $region);
        	        	print LOG "sequence barcode: $barloc\n";
        		} else { #mixed
        	        	$mult_bc = 1;
        	        	print "So. You've mixed barcodes, eh. Tisk, tisk. Alright then, you'll have to tell me, for each fofn file, whether the barcode is on the fwd or rev primer.\n";
        	        	for ($a=0; $a <=$#fofns; $a++) {
        	        		print "For $fofns[$a], [fwd/rev]?\n";
        	        	        while(1) {
        	                		$prim[$a] = <STDIN>;
        	                        	chomp($prim[$a]);
        	                        	if ($prim[$a] ne "fwd" && $prim[$a] ne "rev") {
        	                        		print "Not valid, try again.\n";
        	                        	} else {
        	                                	if ($prim[$a] eq "fwd" || $prim[$a] eq "rev") {
        	                        			($fwd_primer[$a], $rev_primer[$a], $fwd_revcomp_primer[$a], $rev_revcomp_primer[$a], $overlap, $splitlib_min, $splitlib_max) = &getPrimers($prim[$a], $region);
        	                                	        print LOG "sequence barcode for $fofns[$a]: $prim[$a]\n";
        	                                	}
        	                        	last;
        	                        	}
        	                	}
        	        	}
        		}
		} else {
			print "You provided -b with no option.\n";
			print "-b location of the sequence barcode [fwd/rev/mixed] (default: fwd)\n";
                        print "Not a valid response, try again.\n";
                        exit;
		}
	}
	#-q quality filter
	$search_for = "-q";
	my($q_index) = grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($q_index)) {
		if ($q_index+1 < $#opts+1) {
			$qual = $opts[$q_index+1];
			if ($qual < 0 || $qual > 40) {
				print "-q the quality threshold for filtering (default: 30)\n";
				print "Not a valid response, try again.\n";
				exit;
			}
		} else {
			print "You provided -q with no option.\n";
		        print "-q the quality threshold for filtering (default: 30)\n";
			print "Not a valid response, try again.\n";
			exit;
		}
	}
	#-k keep fastq
	$search_for = "-k";
	my($k_index) = grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($k_index)){
		if($k_index+1 < $#opts+1) {
			$keep = $opts[$k_index+1];
			if ($keep ne "y" && $keep ne "n") {
				print "-k keep the .fastq output of the quality filtering process [y/n] (default: n)\n";
				print "Not a valid response, try again.\n";
				exit;
			}
		}
	}
	#-d database
	$search_for = "-d";
        my($d_index) = grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($d_index)) {
		if ($d_index+1 < $#opts+1) {
			$gg = $opts[$d_index+1];
			if($gg ne "gg2011" && $gg ne "gg2013" && $gg ne "silva111" && $gg ne "all") {
				print "-d the greengenes database that you would like to use for taxonomic assignment [gg2011/gg2013/silva111/all] (default: gg2011)\n";
                		print "Not a valid response, try again.\n";
				exit;
			}
		} else {
			print "You provided -d with no option.\n";
			print "-d the greengenes database that you would like to use for taxonomic assignment [gg2011/gg2013/silva111/all] (default: gg2011)\n";
                	print "Not a valid response, try again.\n";
                	exit;
		}
	}
	#-p OTU picking algorithm
	$search_for = "-p";
        my($p_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($p_index)) {
		if ($p_index+1 < $#opts+1) {
			$clust = $opts[$p_index+1];
			#if ($clust ne "abundantotu" && $clust ne "uclust" && $clust ne "cdhit" && $clust ne "dnaclust" && $clust ne "uclust-ref" && $clust ne "mothur" && $clust ne "blast" && $clust ne "uparse" && $clust ne "all") 
			if ($clust ne "abundantotu" && $clust ne "uclust" && $clust ne "cdhit" && $clust ne "dnaclust" && $clust ne "uclust-ref" && $clust ne "uclust-ref-strict" && $clust ne "mothur" && $clust ne "blast" && $clust ne "uparse" && $clust ne "swarm" && $clust ne "usearch61" && $clust ne "all")
			{
				print "-p OTU picking algorithm [abundantotu/uclust/cdhit/dnaclust/uclust-ref/uclust-ref-strict/mothur/blast/uparse/all] (default: abundantotu) \n";
				print "Not a valid response, try again.\n";
				exit;
			}
		} else {
			print "You provided -p with no option.\n";
			print "-p OTU picking algorithm [abundantotu/uclust/cdhit/dnaclust/uclust-ref/uclust-ref-strict/mothur/blast/uparse/all] (default: abundantotu) \n";
	        	print "Not a valid response, try again.\n";
			exit;
		}
	}
	#-c OTU picking clustering threshold
	$search_for = "-c";
	my($c_index) = grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($c_index)) {
		if ($c_index+1 < $#opts+1) {
			$clustThres = $opts[$c_index+1];
			if (!($clustThres >= 0 && $clustThres <= 100)) {
				print "-c OTU picking clustering threshold [0-100] (default: 97) \n";
				print "Not a valid response, try again.\n";
				exit;
			} else { #check to make sure there is a rep set and taxonomy for the associated clustering threshold
				if ($gg eq "gg2011") {
					if (! -e "$db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta") {
						print "There is no corresponding greengenes database references file for the clustering threshold you provided.\n";
						print "Please try again.\n";
						print "-c OTU picking clustering threshold [0-100] (default: 97)\n";
						exit;	
					}
				} elsif ($gg eq "gg2013") {
					my @files = ();
					$files[0] = "$db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta";
					$files[1] = "$db/gg_13_8_otus/taxonomy/".$clustThres."_otu_taxonomy.txt";
					foreach (@files) {
						if (! -e $_) {
							print "There is no corresponding greengenes database reference file for the clustering threshold you provide.\n";
							print "Please try again.\n";
							print "-c OTU picking clustering threshold [0-100] (default: 97)\n";
							exit;
						}
					}
				} elsif ($gg eq "silva111") {
					my @files = ();
					$files[0] = "$db/Silva_111_post/taxonomy/".$clustThres."_Silva_111_taxa_map_RDP_6_levels.txt";
					$files[1] = "$db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta";
					foreach (@files) {
						if (! -e $_) {
							print "There is no corresponding Silva database reference file for the clustering threshold you provide.\n";
							print "Please try again.\n";
							print "-c OTU picking clustering threshold [0-100] (default: 97)\n";
							exit;
						}
					}
				} else { #all
					my @files = ();
					$files[0] = "$db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta";
					$files[1] = "$db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta";
					$files[2] = "$db/gg_13_8_otus/taxonomy/".$clustThres."_otu_taxonomy.txt";
					$files[3] = "$db/Silva_111_post/taxonomy/".$clustThres."_Silva_111_taxa_map_RDP_6_levels.txt";
					$files[4] = "$db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta";
					foreach (@files) {
						if (! -e $_) {
							print "There is at least one missing corresponding database reference file for the clustering threshold you provide.\n";
							print "Please try again.\n";
							print "-c OTU picking clustering threshold [0-100] (default: 97)\n";
							exit;
						}
					}
				}
			}
		} else {
			print "You provided -c with no option.\n";
			print "-c OTU picking clustering threshold [0-100] (default: 97) \n";
			print "Not a valid response, try again.\n";
			exit;
		}
	}
	#-m chimera checking
	$search_for = "-m";
	my($m_index) = grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($m_index)) {
		if ($m_index+1 < $#opts+1) {
			$chimeraChecking = $opts[$m_index+1];
			if ($chimeraChecking ne "y" && $chimeraChecking ne "n") {
				print "-m whether the user would like to employ chimera checking [y/n] (default: y)\n";
				print "Not a valid response, try again.\n";
				exit; 
			}
		} else {
			print "You provided -m with no option.\n";
			print "-m whether the user would like to employ chimera checking [y/n] (default: y)\n";
			print "Not a valid response, try again.\n";
			exit;
		}
	}
	#-l linkage for mothur
	$search_for = "-l";
        my($l_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($l_index)) {
		if ($l_index+1 < $#opts+1) {
			$linkage = $opts[$l_index+1];
			if($linkage ne "nearest" && $linkage ne "average" && $linkage ne "furthest" && $linkage ne "all") {
                		print "-l the linkage method to use with mothur [nearest/average/furthest/all] (default: average)\n";
                	        print "Not a valid response, try again.\n";
				exit;
                	}
		} else {
			print "You provided -l with no option.\n";
			print "-l the linkage method to use with mothur [nearest/average/furthest] (default: average)\n";
                        print "Not a valid response, try again.\n";
			exit;
		}
	}	
	#-t taxonomic assignment method
	$search_for = "-t";
        my($t_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($t_index)) {
		if ($t_index+1 < $#opts+1) {
			$taxon = $opts[$t_index+1];
			if($taxon ne "blast" && $taxon ne "rdp-training" && $taxon ne "all") {
				print "-t the taxonomic assignment method [blast/rdp-training/all] (default: rdp-training)\n";
                        	print "Not a valid response, try again.\n";
				exit;
			}
		} else {
			print "You provided -t with no option.\n";
			print "-t the taxonomic assignment method [blast/rdp-training/all] (default: rdp-training)\n";
                        print "Not a valid response, try again.\n";
			exit;
		}
	}
	#-e taxonomic assignment classification cutoff
	$search_for = "-e";
	my($e_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($e_index)) {
		if ($e_index+1 < $#opts+1) {
			$taxon_cutoff = $opts[$e_index+1];
			if($taxon_cutoff < 0 || $taxon_cutoff > 1) {
				print "-e certainty cutoff for taxonomic assignment as designed by Qiime [0-1] (default: 0.5)\n";
				print "Not a valid response, try again.\n";
				exit;
			}
		} else {
			print "You provided -e with no option.\n";
			print "-e certainty cutoff for taxonomic assignment as designed by Qiime [0-1] (default: 0.5)\n";
			print "Not a valid response, try again.\n";
			exit;
		}
	}
	#-x use timed code
	$search_for = "-x";
        my($x_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
        if (defined($x_index)) {
		if ($x_index+1 < $#opts+1) {
                	$time = $opts[$x_index+1];
			if($time ne "y" && $time ne "n") {
				#print "-x do you want timed results for each step of the analysis [y/n] (default: n) \nWARNING: the timed code is slightly outdated- use with caution and report any adverse behaviour to whelanfj\@mcmaster.ca\n";
				print "-x do you want timed OTU picking and overall results [y/n] (default: n)\n";
				print "Not a valid response, try again.\n";
				exit;
			}
		} else {
			print "You provided -x with no option.\n";
			#print "-x do you want timed results for each step of the analysis [y/n] (default: n) \nWARNING: the timed code is slightly outdated- use with caution and report any adverse behaviour to whelanfj\@mcmaster.ca\n";
			print "-x do you want timed OTU picking and overall results [y/n] (default: n)\n";
                        print "Not a valid response, try again.\n";
                        exit;

		}
        }
	#-f force overwrite
	$search_for = "-f";
        my($f_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($f_index)) {
		$force = 1;
	}
	#-u CPU threads
	$search_for = "-u";
	my($u_index)= grep { $opts[$_] eq $search_for } 0..$#opts;
	if (defined($u_index)) {
		if ($u_index+1 < $#opts+1) {
			$threads = $opts[$u_index+1];;
			if ($threads !~/^\d+/) {
				print "-u how many CPU threads you would like to envoke [0-9] (default: 1)\n";
                        	print "Not a valid response, try again.\n";
                        	exit;
			}
		} else {
			print "You provided -u with no option.\n";
			print "-u how many CPU threads you would like to envoke [0-9] (default: 1)\n";
			print "Not a valid response, try again.\n";
			exit;
		}
		
	}
}
###############################################################
#Check for possible overwrites
###############################################################
print LOG "###############################################################\n";
print LOG "#Check for possible overwrites\n";
print LOG "###############################################################\n";
$pwd = `pwd`;
my $cmd;
chomp($pwd);
my $otu_folder;
if (($clust eq "all") && ($gg eq "all")) {
	#abundantotu
	$otu_folder = "picked_otus_abundantotu_gg2011";
	chomp($otu_folder);
	if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_abundantotu_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	$otu_folder = "picked_otus_abundantotu_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_abundantotu_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_abundantotu_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_abundantotu_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#blast
	$otu_folder = "picked_otus_blast_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_blast_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_blast_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_blast_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_blast_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_blast_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#cdhit
	$otu_folder = "picked_otus_cdhit_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_cdhit_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_cdhit_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_cdhit_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_cdhit_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_cdhit_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#uclust
	$otu_folder = "picked_otus_uclust_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_uclust_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_uclust_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#uclust-ref
	$otu_folder = "picked_otus_uclust-ref_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust-ref_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_uclust-ref_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust-ref_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_uclust-ref_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust-ref_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#uclust-ref-strict
	$otu_folder = "picked_otus_uclust-ref-strict_gg2011";
	chomp($otu_folder);
	if (-d $otu_folder) {
		if ($force) {
			$cmd = "rm -rf picked_otus_uclust-ref-strict_gg2011";
			system($cmd);
		} else {
			print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
			exit;
		}
	}
	$otu_folder = "picked_otus_uclust-ref-strict_gg2013";
	chomp($otu_folder);
	if (-d $otu_folder) {
		if ($force) {
			$cmd = "rm -rf picked_otus_uclust-ref-strict_gg2013";
			system($cmd);
		} else {
			print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
			exit;
		}
	}
	$otu_folder = "picked_otus_uclust-ref-strict_silva111";
	chomp($otu_folder);
	if (-d $otu_folder) {
		if ($force) {
			$cmd = "rm -rf picked_otus_uclust-ref-strict_silva111";
			system($cmd);
		} else {
			print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
			exit;
		}
	}
	#dnaclust
	$otu_folder = "picked_otus_dnaclust_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_dnaclust_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_dnaclust_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_dnaclust_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_dnaclust_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_dnaclust_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#mothur_average
	$otu_folder = "picked_otus_mothur_average_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_mothur_average_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_mothur_average_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_mothur_average_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_mothur_average_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_mothur_average_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	#uparse
	$otu_folder = "picked_otus_uparse_gg2011";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uparse_gg2011";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_uparse_gg2013";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uparse_gg2013";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }
	$otu_folder = "picked_otus_uparse_silva111";
	chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uparse_silva111";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }
        }	
} elsif($clust eq "all") {
	#abundantotu
	$otu_folder = "picked_otus_abundantotu_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_abundantotu_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#blast
	$otu_folder = "picked_otus_blast_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_blast_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#cdhit
	$otu_folder = "picked_otus_cdhit_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_cdhit_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#uclust
	$otu_folder = "picked_otus_uclust_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#uclust-ref
	$otu_folder = "picked_otus_uclust-ref_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uclust-ref_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#uclust-ref-strict
	$otu_folder = "picked_otus_uclust-ref-strict_$gg";
	chomp($otu_folder);
	if (-d $otu_folder) {
		if ($force) {
			$cmd = "rm -rf picked_otus_uclust-ref-strict_$gg";
			system($cmd);
		} else {
			print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
			exit;
		}
	}
	#dnaclust
	$otu_folder = "picked_otus_dnaclust_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_dnaclust_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#mothur_average
	$otu_folder = "picked_otus_mothur_average_$gg";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_mothur_average_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }
	#mothur_nearest
	#$otu_folder = "picked_otus_mothur_nearest_*/";
        #chomp($otu_folder);
        #if (-d $otu_folder) {
        #        if ($force) {
        #                $cmd = "rm -rf picked_otus_mothur_nearest_*";
        #                system($cmd);
        #        } else {
        #                print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
        #                exit;
        #        }    
        #}
	#mothur_furthest
	#$otu_folder = "picked_otus_mothur_furthest_*/";
        #chomp($otu_folder);
        #if (-d $otu_folder) {
        #        if ($force) {
        #                $cmd = "rm -rf picked_otus_mothur_furthest_*";
        #                system($cmd);
        #        } else {
        #                print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
        #                exit;
        #        }    
        #}
	#uparse
	$otu_folder = "picked_otus_uparse_*/";
        chomp($otu_folder);
        if (-d $otu_folder) {
                if ($force) {
                        $cmd = "rm -rf picked_otus_uparse_$gg";
                        system($cmd);
                } else {
                        print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
                        exit;
                }    
        }	
} else {
	$otu_folder = "picked_otus_".$clust."_$gg";
	chomp($otu_folder);
	if (-d $otu_folder) {
		if ($force) {
			$cmd = "rm -rf picked_otus_".$clust."_$gg";
			system($cmd);
		} else {
			print "Folder <picked_otus_> already exists. Add -f to overwrite.\n";
			exit;
		}
	}
}
my $splits_folder = "splits_dir/";
chomp($splits_folder);
if (-d $splits_folder) {
	if ($force) {
		$cmd = "rm -rf splits_dir";
		system($cmd);
	} else {
		print "Folder <splits_dir> already exists. Add -f to overwrite.\n";
		exit;
	}
}
if (-e "$proj.fna") {
	if ($force) {
		system("rm $proj.fna");
	} else {
		print "File $proj.fna already exists. Add -f to overwrite.\n";
		exit;
	}
}
if (-e "map_$proj.txt") {
	if ($force) {
		system("rm map_$proj.txt");
	} else {
		print "File map_$proj.txt already exists. Add -f to overwrite.\n";
		exit;
	}
}
###############################################################
#Pre-check for software installations/version compatibility
###############################################################
print LOG "###############################################################\n";
print LOG "#Pre-check for software installations/version compatibility\n";
print LOG "###############################################################\n";
#region; #barloc; #mult_bc; 
#$qual;
my $exc;
$exc = `$bin/sickle/sickle --version`;
if ($exc!~/sickle version/) {
	print "sickle is not installed; please run `perl SETUP.pl` to install.\n";
	exit;
}
#$gg;
if (($gg eq "gg2011") || ($gg eq "all")) {
	if (! -d "$db/gg_otus_4feb2011/rep_set/") {
		print "Greengenes 2011 is not installed; please see instructions in INSTALL.txt or choose a different database.\n";
		exit;
	}
} elsif (($gg eq "gg2013") || ($gg eq "all")) {
	if (! -d "$db/gg_13_8_otus/rep_set/") {
		print "Greengenes 2013 is not installed; please see instructions in INSTALL.txt or choose a different database.\n";
		exit;
	}
} elsif (($gg eq "silva111") || ($gg eq "all")) {
	if (! -d "$db/Silva_111_post/rep_set/") {
		print "Silva 111 is not installed; please see instructions in INSTALL.txt or choose a different database.\n";
		exit;
	} 
}
#$clust; #$clustThres; #$linkage;
if (($clust eq "abundantotu") || ($clust eq "all")) {
	if (! -e "$bin/AbundantOTU+0.93b/bin/AbundantOTU+") {
		print "AbundantOTU+0.93b is not installed; please run `perl SETUP.pl` to install or choose a different OTU picking algorithm.\n";
		exit;
	}
} elsif (($clust eq "uclust") || ($clust eq "cdhit") || ($clust eq "uclust-ref") || ($clust eq "uclust-ref-strict")  || ($clust eq "blast") || ($clust eq "all")) {
	$exc = `print_qiime_config.py`;
	if ($exc =~/QIIME library version:\s*Not installed\./) {
		print "Qiime is not installed; please see instructions in INSTALL.txt\n";
		exit;
	}
} elsif (($clust eq "dnaclust") || ($clust eq "all")) {
	if (! -e "$bin/dnaclust_linux_release3/dnaclust") {
		print "dnaclust is not installed; please run `perl SETUP.pl` to install or choose a different OTU picking algorithm.\n";
		exit;
	}
} elsif (($clust eq "mothur") || ($clust eq "all")) {
	$exc = `print_qiime_config.py -tf`;
	if ($exc =~/FAIL: test_mothur_supported_version/) {
		print "mothur is not installed; please see instructions in INSTALL.txt or choose a different OTU picking algorithm.\n";
		exit;
	}
	#Remove unnecessary logfile created by -tf
	`rm mothur.*.logfile`;
} elsif (($clust eq "uparse") || ($clust eq "all")) {
	if (! -e "$bin/usearch") {
		print "usearch is not installed; please see instructions in INSTALL.txt or choose a different OTU picking algorithm.\n";
		exit;
	}
}
#$taxon;#$time; #$force; #$threads;
#$chimeraChecking;
if ($chimeraChecking eq "y") {
	$exc = `print_qiime_config.py -tf`;
	if ($exc =~/FAIL:  test_usearch_supported_version/) {
		print "usearch is not installed and visible to Qiime; please see instructions in INSTALL.txt or choose -m n.\n";
		exit;
	}
	#Remove unncessary logfile created by -tf
	`rm mothur.*.logfile`;
}
###############################################################
#Record QIIME Information to Log File
###############################################################
print LOG "###############################################################\n";
print LOG "#Record QIIME Information\n";
print LOG "###############################################################\n";
my $info = `print_qiime_config.py`;
print LOG $info."\n";

###############################################################
#Ensure Illumina Naming Scheme
#Clip with cutadapt
###############################################################
print"###############################################################\n";
print"Clip with cutadapt\n";
print"###############################################################\n";
print LOG "###############################################################\n";
print LOG "#Clip with cutadapt\n";
print LOG "###############################################################\n";
my $illumina;
for($a = 0; $a < $#fofns+1; $a++) {
	open (FOFN, "<", $fofns[$a]) or die "$!";
	my @files = <FOFN>;
	my $filecount = 0;	
	#mkdir
	unless (-e "pandaseq_logs$a") {
		$cmd = "mkdir pandaseq_logs$a/";
		system($cmd);
	}
	while($filecount < $#files+1) {
		my $in = $files[$filecount];
		chomp($in);
		#Check if file $in exists; exit gracefully if it doesn't
		if (!(-e $in)) {
			print "File $in in user input fofn doesn't not exist.\n";
			print "Exiting...\n";
			exit;
		}
		$filecount += 1;
		#Ensure Illumina File Naming Scheme
		if ($in=~/((.*?_[ACGT]*?_L\d{3}_R\d_\d{3})\.fastq(.gz|))$/) { #old scheme
			#print $in."\n";
			#$in = $1;
			#print $in."\n";
			#$out = "pandaseq_logs$a/$in";
			$illumina = "old";
		} elsif ($in=~/((.*?_S\d+_L\d{3}_R\d_\d{3})\.fastq(.gz|))$/) { #new scheme
			$illumina = "new";
		} else {
			print "Error: $in is not a properly formatted Illumina .fastq file- aborting.\n";
			exit;
		}
		$in=~/.*\/(.*?_.*?_L\d{3}_R\d_\d{3}\.fastq(.gz|))$/;
		my $path = $1;
		print $path."\n";
		my $out = "pandaseq_logs$a/$path";
		#cutadapt
		if ($in=~/.*_R1_.*/) {
			if ($mult_bc) {
				print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $rev_revcomp_primer[$a] $in 2>&1 | tee -a $err\n";
				`$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $rev_revcomp_primer[$a] $in 2>&1 | tee -a $err`;
			} else {
				print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $rev_revcomp_primer $in 2>&1 | tee -a $err\n";
				`$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $rev_revcomp_primer $in 2>&1 | tee -a $err`;
			}
		} elsif ($in=~/.*_R2_.*/) {
			if ($mult_bc) {
				print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $fwd_revcomp_primer[$a] $in 2>&1 | tee -a $err\n";
				`$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $fwd_revcomp_primer[$a] $in 2>&1 | tee -a $err`;
			} else {
				print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $fwd_revcomp_primer $in 2>&1 | tee -a $err\n";
				`$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 17 -m 0 -o $out -a $fwd_revcomp_primer $in 2>&1 | tee -a $err`;
			}
		} else {
			print "Error: $in is not a R1 or R2 file- aborting.\n";
			exit;
		}
	}
	close FOFN;
}
###############################################################
#Align with PANDAseq
###############################################################
print"###############################################################\n";
print"Align with PANDAseq\n";
print"###############################################################\n";
print LOG "###############################################################\n";
print LOG "#Align with PANDAseq\n";
print LOG "###############################################################\n";
my @sample_names = (); my @orig_reads = (); my @post_pandaseq = (); my @post_cutadapt = (); my @post_sickle = ();
my $samplef; my $sampler; my $barcodef; my $barcoder;
my $fwd; my $rev;
my $barco; my $head;
for($b = 0; $b < $#fofns+1; $b++) {
	my $dir = "pandaseq_logs$b";
	chdir($dir);
	my $mapf = "map_".$proj."".$b.".txt";
	my $fast = $proj."".$b.".fa";
	open (FOFN, "<", "../".$fofns[$b]) or die "$!";
	open (MAPF, ">", $mapf) or die "$!";
	#Set up mapfile header
	print MAPF "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n";
	#Set up global fasta file; make sure it's empty
	if (my ($grabbing) = glob($fast)) {
		my $cmd0 = "rm $fast";
		system($cmd0);
	}
	$cmd = "touch $fast";
	system($cmd);
	#Run pandaseq on each set of sequences; generate mapFile
	my @data = <FOFN>;
	my $i = 0;
	while ($i < $#data+1) {
		$samplef = $sampler = $barcodef = $barcoder = "";
		$head = "";
		$fwd = $data[$i];
		chomp($fwd);
		$rev = $data[$i+1];
		chomp($rev);
		$i += 2;
		if ($illumina eq "old") {
			$fwd=~/(([^\/]+?)_([ACGT]+)_L001_R1_001\.fastq(.gz|))/;
			$fwd = $1 and $samplef = $2 and $barcodef = $3;
			$rev=~/(([^\/]+?)_([ACGT]+)_L001_R2_001\.fastq(.gz|))/;
			$rev = $1 and $sampler = $2 and $barcoder = $3;
		} elsif ($illumina eq "new") {
			$fwd=~/(([^\/]+?)_S\d+_L001_R1_001\.fastq(.gz|))/;
			$fwd = $1 and $samplef = $2;
			if ($3=~/gz/) {
				$head = `gzip -cd $fwd | head -n 1`;
			} else {
				$head = `head -n 1 $fwd`;
			}
			$head=~/\@.+:\d+:.+:\d+:\d+:\d+:\d+ 1:N:0:([ACGT]+)\+[ACGT]+/;
			$barcodef = $1;
			$rev=~/(([^\/]+?)_S\d+_L001_R2_001\.fastq(.gz|))/;
			$rev = $1 and $sampler = $2;
			if ($3=~/gz/) {
				$head = `gzip -cd $rev | head -n 1`;
			} else {
				$head = `head -n 1 $rev`;
			}
			$head=~/\@.+:\d+:.+:\d+:\d+:\d+:\d+ 2:N:0:([ACGT]+)\+[ACGT]+/;
			$barcoder = $1;
		}
		print $fwd."\n";
		print "$samplef\t$barcodef\n";
		#Check
		if ($samplef ne $sampler || $barcodef ne $barcoder || $samplef eq "" || $sampler eq "" || $barcodef eq "" || $barcoder eq "") {
			print "Regex for fwd and rev do not match:\n";
			print "samplef=".$samplef."\tsampler=".$sampler."\n";
			print "barcodef=".$barcodef."\tbarcoder=".$barcoder."\n";
			exit;
		}
		#Make mapFile- check SampleID for illegal characters
		$samplef=~s/-|\+|\s/\./g;
		print $samplef."\n";
		print MAPF $samplef."\t".$barcodef."\tCCTACGGGAGGCAGCAG\t$proj\n";
		$barco = length $barcodef;
		#Run pandaseq cmd
		$fwd=~/([^\/]+?)_L001_R1_001\.fastq(.gz|)/;
		my $fastq = $1.".fastq";
		my $pandalog = $1.".log.bz2";
		$cmd = "";
		#multiple barcode locations?
		if ($mult_bc) {
			$cmd = "(pandaseq -A simple_bayesian -f ".$fwd." -r ".$rev." -p $fwd_primer[$b] -q $rev_primer[$b] -N -t 0.70 -T 1 -F > ".$fastq.") 2>&1| bzip2 -c > ".$pandalog;
			#add -o 10 to specify a minimum overlap of 10 bases
		} elsif ($region eq "other") {
			$cmd = "(pandaseq -A simple_bayesian -f ".$fwd." -r ".$rev." -N -t 0.70 -T 1 -F > ".$fastq.") 2>&1| bzip2 -c > ".$pandalog;
		} else {
			$cmd = "(pandaseq -A simple_bayesian -f ".$fwd." -r ".$rev." -p $fwd_primer -q $rev_primer -N -t 0.70 -T 1 -F > ".$fastq.") 2>&1| bzip2 -c > ".$pandalog;
		}
		
		my $cmd2 = "bzgrep STAT ".$pandalog." | tail -n 9";
		print $cmd."\n";
		print LOG $cmd."\n";
		system($cmd);
		print $cmd2."\n";
		print LOG $cmd2."\n";
		system($cmd2);
		#Get original read information from wc -l
		my $lines;
		if ($fwd=~/.*\.fastq\.gz/) {
			$lines = `zcat $fwd | wc -l`;
		} else {
			$lines = `wc -l $fwd`;
		}
		$lines =~/(\d+).*/; $lines = $1;
		$lines = $lines/4;
		push @orig_reads, $lines;
		#Save sample name, original, and post pandaseq number of reads
		$samplef=~s/_/./g;
		push @sample_names, $samplef;
		$cmd = `bzgrep \"STAT\\sREADS\" $pandalog  | tail -n 1`;
		$cmd =~/.*\t(\d*)$/;
		#push @orig_reads, $1;
		$cmd = `bzgrep \"STAT\\sOK\" $pandalog | tail -n 1`;
		$cmd =~/.*\t(\d*)$/;
		push @post_pandaseq, $1;
		#Delete any reads that still have Illumina sequencing primers in them.
		#delete reads with the binding site for the Illumina sequencing primers, V3_F primer; -O 5 came from research
        	print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 5 -m 0 -a acactctttccctacacgacgctcttccgatct --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a $err\n";
        	system("$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 5 -m 0 -a acactctttccctacacgacgctcttccgatct --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a ../temp.txt");
        	#$cmd = `grep "Processed reads" ../temp.txt | tail -n 1`;
        	$cmd = `grep "Total reads processed" ../temp.txt | tail -n 1`;
		$cmd =~/.*?([\d,]+).*?/;
        	my $all = $1;
		$all=~s/,//g;
        	$cmd = `grep "Reads with adapters" ../temp.txt | tail -n 1`;
        	$cmd =~/.*?([\d,]+).*?\(.*?\)/;
        	my $lost = $1;
		$lost=~s/,//g;
        	my $left = $all - $lost;
		print LOG "cat ../temp.txt >> $err\n";
        	`cat ../temp.txt >> $err`;
        	`rm ../temp.txt`;
        	print LOG "mv tmp.fa $fastq\n";
        	`mv tmp.fa $fastq`;
        	#delete reads with the REVCOMP of the binding site for the Illumina sequencing primers, V3_F primer; -O 8 came from research
        	print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 8 -m 0 -a agatcggaagagcgtcgtgtagggaaagagtgt --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a $err\n";
        	system("$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 8 -m 0 -a agatcggaagagcgtcgtgtagggaaagagtgt --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a ../temp.txt");
		$cmd = `grep "Total reads processed" ../temp.txt | tail -n 1`;
                $cmd =~/.*?([\d,]+).*?/;
                $all = $1;
		$all=~s/,//g;
                $cmd = `grep "Reads with adapters" ../temp.txt | tail -n 1`;
		$cmd =~/.*?([\d,]+).*?\(.*?\)/;
        	$lost = $1;
		$lost=~s/,//g;
        	#$left = $left + ($all - $lost);
		$left = $left - $lost;
        	print LOG "cat ../temp.txt >> $err\n";
        	`cat ../temp.txt >> $err`;
        	`rm ../temp.txt`;
        	print LOG "mv tmp.fa $fastq\n";
        	`mv tmp.fa $fastq`;
        	#delete reads with the binding site for the Illumina sequencing primers, V3_R primer; -O 11 came from research
        	print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 11 -m 0 -a gtgactggagttcagacgtgtgctcttccgatct --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a $err\n";
        	system("$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 11 -m 0 -a gtgactggagttcagacgtgtgctcttccgatct --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a ../temp.txt");
        	$cmd = `grep "Total reads processed" ../temp.txt | tail -n 1`;
                $cmd =~/.*?([\d,]+).*?/;
                $all = $1;
		$all=~s/,//g;
                $cmd = `grep "Reads with adapters" ../temp.txt | tail -n 1`;
		$cmd =~/.*?([\d,]+).*?\(.*?\)/;
        	$lost = $1;
		$lost=~s/,//g;
        	#$left = $left + ($all - $lost);
		$left = $left - $lost;
        	print LOG "cat ../temp.txt >> $err\n";
        	`cat ../temp.txt >> $err`;
        	`rm ../temp.txt`;
        	print LOG "mv tmp.fa $fastq\n";
        	`mv tmp.fa $fastq`;
        	#delete reads with the REVCOMP of the binding site for the Illumina sequencing primers, V3_R primer; -O 8 came from research
        	print LOG "$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 8 -m 0 -a agatcggaagagcacacgtctgaactccagtcac --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a $err\n";
        	system("$bin/cutadapt-1.8.1/bin/cutadapt -f fastq -e 0 -n 1 -O 8 -m 0 -a agatcggaagagcacacgtctgaactccagtcac --discard-trimmed -o tmp.fa $fastq 2>&1 | tee -a ../temp.txt");
        	$cmd = `grep "Total reads processed" ../temp.txt | tail -n 1`;
                $cmd =~/.*?([\d,]+).*?/;
                $all = $1;
		$all=~s/,//g;
                $cmd = `grep "Reads with adapters" ../temp.txt | tail -n 1`;
		$cmd =~/.*?([\d,]+).*?\(.*?\)/;
        	$lost = $1;
		$lost=~s/,//g;
        	#$left = $left + ($all - $lost);
		$left = $left - $lost;
		print LOG "cat ../temp.txt >> $err\n";
		`cat ../temp.txt >> $err`;
		`rm ../temp.txt`;
		print LOG "mv tmp.fa $fastq\n";
		`mv tmp.fa $fastq`;
        	push @post_cutadapt, $left;
		#Remove low quality paired-end reads post pandaseq
		print LOG "$bin/sickle/sickle se -q $qual -l $splitlib_min -f $fastq -t sanger -o tmp.fa 2>&1 | tee -a ../temp.txt\n";
		system("$bin/sickle/sickle se -q $qual -l $splitlib_min -f $fastq -t sanger -o tmp.fa 2>&1 | tee -a ../temp.txt");
		print "$bin/sickle/sickle se -q $qual -l $splitlib_min -f $fastq -t sanger -o tmp.fa 2>&1 | tee -a ../temp.txt\n";
		`mv tmp.fa $fastq`;
		#Calculate the number of sequences post sickle
		$left = `grep "FastQ records kept:" ../temp.txt`;
		$left=~/FastQ records kept:\s(\d*)/;
		push @post_sickle, $1;
		#Convert fastq to fasta
		`grep -A 1 "^\@M" $fastq > tmp.fa`;
		`perl -i -pe 's/\@/>/g' tmp.fa`;
		`perl -i -pe 's/--\n//g' tmp.fa`;
		#Cat results to global fasta file
		$cmd = "cat tmp.fa >> $fast";
		print LOG $cmd."\n";
		system($cmd);
		`rm tmp.fa`;
		#Cat log results to $err
		print LOG "cat ../temp.txt >> $err\n";
		`cat ../temp.txt >> $err`;
		`rm ../temp.txt`;
	}
	close FOFN;
	close MAPF;
	#rid of underscores in mapfile
	`perl -pi -e 's/_/\./g' $mapf`;
	chdir('..');
}
#print to seq info file
print SEQ "sample name\t".join("\t", @sample_names)."\n";
print SEQ "original\t".join("\t", @orig_reads)."\n";
print SEQ "post_pandaseq\t".join("\t", @post_pandaseq)."\n";
print SEQ "post_cutadapt\t".join("\t", @post_cutadapt)."\n";
print SEQ "post_sickle\t".join("\t", @post_sickle)."\n";
###############################################################
###############################################################
print"###############################################################\n";
print "Flip any Forward Barcoded Reads\n";
print"###############################################################\n";
print LOG"###############################################################\n";
print LOG "#Flip any Forward Barcoded Reads\n";
print LOG"###############################################################\n";
my $fi; my $dir; my $inf; my $line; my $revcomp;
if ($barloc eq "rev") {
	#perfect, don't have to do anything
} elsif ($barloc eq "fwd") {
	for ($fi = 0; $fi < $#fofns+1; $fi++) {
		$dir = "pandaseq_logs$fi";
		chdir($dir);
		$inf = $proj."".$fi.".fa";
		#flip reads
		open (INF, "<", $inf) or die "$!";
		open (TMP, ">", "temp.fa") or die "$!";
		foreach $line (<INF>) {
			if ($line=~/>.*/) {
				print TMP $line;
			} else {
				chomp($line);
				$revcomp = reverse($line);
				$revcomp=~tr/ACGTacgt/TGCAtgca/;
				print TMP $revcomp."\n";	
			}
		}
		close INF;
		close TMP;
		`mv temp.fa $inf`;
		chdir('..');
	}
} else { #mixed {
	for ($fi = 0; $fi < $#fofns+1; $fi++) {
		if ($prim[$fi] eq "fwd") {
			$dir = "pandaseq_logs$fi";
			chdir($dir);
			$inf = $proj."".$fi.".fa";
			#flip reads
			open (INF, "<", $inf) or die "$!";
			open (TMP, ">", "temp.fa") or die "$!";
			foreach $line (<INF>) {
				if ($line=~/>.*/) {
					print TMP $line;
				} else {
					chomp($line);
					$revcomp = reverse($line);
					$revcomp=~tr/ACGTacgt/TGCAtgca/;
					print TMP $revcomp."\n";
				}
			}
			close INF;
			close TMP;
			`mv temp.fa $inf`;
			chdir('..');
		}
	}
}
###############################################################
#Hack Barcodes for input into QIIME
###############################################################
print"###############################################################\n";
print"Hack Barcodes for input into QIIME\n";
print"###############################################################\n";
print LOG "###############################################################\n";
print LOG "#Hack Barcodes for input into QIIME\n";
print LOG "###############################################################\n";
#Add barcodes?
my $c; my $add; my $map; my $otf; my @data; my $i; my $barcode; my $m; my @mapdata; my $newline;
for($c = 0; $c < $#fofns+1; $c++) {
	$add = "";
	if ($numruns > 64) {
		print "sl1p currently can only handle 64 fofn files.\n";
		print "Exiting...\n";
		exit;
	}
	if ($numruns > 16) {# && $numruns <= 64) {
		if (int($c/16) == 0) { $add = $add."A"; }
		if (int($c/16) == 1) { $add = $add."C"; }
		if (int($c/16) == 2) { $add = $add."G"; }
		if (int($c/16) == 3) { $add = $add."T"; }
	}
	if ($numruns > 4) {# && $numruns <= 16) {
		if (int($c/4%4) == 0) { $add = $add."A"; }
		if (int($c/4)%4 == 1) { $add = $add."C"; }
		if (int($c/4)%4 == 2) { $add = $add."G"; }
		if (int($c/4)%4 == 3) { $add = $add."T"; }
	}
	if ($numruns > 1) {# && $numruns <= 4) {
		if ($c%4 == 0) { $add = $add."A"; }
		if ($c%4 == 1) { $add = $add."C"; }
		if ($c%4 == 2) { $add = $add."G"; }
		if ($c%4 == 3) { $add = $add."T"; }
	}
	$dir = "pandaseq_logs$c";
	chdir($dir);
	$map = "map_".$proj."".$c.".txt";
	$inf = $proj."".$c.".fa";
	$otf = $proj."".$c.".fna";
	open (IN, "<", $inf) or die "$!";
	open (OUT, ">", $otf) or die "$!";
	@data = <IN>;
	$i = 0;
	while ($i < $#data+1) {
		$line = $data[$i];
		#$sub = $barco * -1 + 1;
		if ($illumina eq "old") {
			$barcode = substr $line, -7;
		} elsif ($illumina eq "new") {
			$barcode = substr $line, -14, 6;
		}
		chomp($barcode);
		if ($add ne "") {
			$barcode = $add.$barcode;
		}
		chomp($line);
		#$sub = $barco * -1;
		if ($illumina eq "old") {
			$line = substr $line, 0, -6;
		} elsif ($illumina eq "new") {
			$line = substr $line, 0, -13;
		}
		print OUT $line.$barcode."\n";
		$i += 1;
		$line = $data[$i];
		print OUT $barcode.$line;
		$i += 1;
	}
	close IN;
	close OUT;
	#update mapfile
	if ($add ne "") {
		open (MAP, "<", $map) or die "$!";
		@mapdata = <MAP>;
		$m = 1; #skip header line
		while ($m <= $#mapdata) {
			$line = $mapdata[$m];
			$line =~/(.*?)\t([ACGT]*?)\t(.*)/;
			$newline = $1."\t".$add."".$2."\t".$3."\n";
			system("perl -pi -e 's/$line/$newline/' $map");
			$m += 1;
		}
	}
	close MAP;
	chdir('..');
}
###############################################################
#Consolidate run files
###############################################################
print"###############################################################\n";
print"Consolidate run files\n";
print"###############################################################\n";
print LOG "###############################################################\n";
print LOG "#Consolidate run files\n";
print LOG "###############################################################\n";
mkdir ("pandaseq_logs");
my $d; my $mapf; my $fast; my $cmd2;
for($d = 0; $d < $#fofns+1; $d++) {
	$dir = "pandaseq_logs$d";
	chdir($dir);
	$mapf = "map_".$proj."".$d.".txt";
	$fast = $proj."".$d.".fna";
	$cmd = "cp $fast ../$fast";
	$cmd2 = "cp $mapf ../$mapf";
	system($cmd);
	system($cmd2);
	#If map file isn't the first, remove header
	if ($d != 0) {
		$cmd = "perl -i -pe 's/#SampleID\\tBarcodeSequence\\tLinkerPrimerSequence\\tDescription\\n//g' ../$mapf";
		system($cmd);
	}
	$cmd = "rm $fast";
	$cmd2 = "rm $mapf";
	system($cmd);
	system($cmd2);
	#$cmd = "rm *\.fastq";
	#system($cmd);
	$cmd = "rm *\.fa";
	system($cmd);
	#Move log files into the main pandaseq_logs folder
	$cmd = "mv *.log.bz2 ../pandaseq_logs";
	system($cmd);
	#If -k flag used, keep fastq; else remove
	if ($keep eq "y") {
		$cmd = "mv *\.fastq ../pandaseq_logs";
	} else {
		$cmd = "rm *\.fastq";
	}
	system($cmd);
	chdir('..');
}
#Remove all numbered, sub pandaseq_log folders
my $cmd0; my $cmd1;
$cmd = "rm -r pandaseq_logs[0-9]*";
system($cmd);
if ( -e "map_".$proj.".txt") {
	$cmd0 = "rm map_".$proj.".txt";
	system($cmd0);
}
if ( -e $proj.".fna") {
	$cmd1 = "rm ".$proj.".fna";
	system($cmd1);
}
$cmd = "cat map_".$proj."[0-9]*.txt > map_".$proj.".txt";
$cmd2 = "cat ".$proj."[0-9]*.fna > ".$proj.".fna";
system($cmd);
system($cmd2);
$cmd = "rm $proj\[0-9\]*.fna";
system($cmd);
$cmd = "rm map_$proj\[0-9\]*.txt";
system($cmd);
###############################################################
#Run QIIME
###############################################################
print"###############################################################\n";
print"Run QIIME\n";
print"###############################################################\n";
print LOG "###############################################################\n";
print LOG "#Run QIIME\n";
print LOG "###############################################################\n";
#Setup variables
$mapf = "map_".$proj.".txt";
$fast = $proj.".fna";
if ($numruns == 1) {
	$barco = $barco; #6;
} elsif ($numruns <= 4) {
	$barco = $barco+1; #7;
} elsif ($numruns <= 16) {
	$barco = $barco+2; #8;
} else {
	$barco = $barco+3; #9;
}
#Write qiime_params.txt file
open ALPH, ">", "qiime_params.txt" or die "$!\n";
print ALPH "alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species\n";
print ALPH "beta_diversity:metrics weighted_unifrac,unweighted_unifrac,bray_curtis,abund_jaccard,binary_jaccard\n";
print ALPH "plot_taxa_summary:chart_type bar\n";
close ALPH;

#check_id_map
print "check_id_map.py\n";
my $out = `check_id_map.py -m $mapf -o map_aux`;
chomp($out);
if ($out=~/This script has been renamed validate_mapping_file\.py for clarity/) {
	$out = `validate_mapping_file.py -m $mapf -o map_aux`;
}
chomp($out);
if ($out ne "No errors or warnings were found in mapping file.") {
	print "Error in a mapfile $mapf\n";
	exit;
}
print LOG "tar -zcvf map_aux.tar.gz map_aux/\n";
`tar -zcvf map_aux.tar.gz map_aux/`;
print LOG "rm -rf map_aux/\n";
`rm -rf map_aux/`;
#split_libraries
my $splits = 1;
my @post_splitlib = ();
if ($splits) {
	print "split_libraries.py\n";
	$cmd = "";
	#if ($region eq "v3") {
	$cmd="split_libraries.py -m $mapf -f $fast -l $splitlib_min -L $splitlib_max -b $barco -M 25 -k -e 0 -H 20 -n 1 -o splits_dir/ 2>&1 | tee -a $err";
	if ($time eq "y") { $cmd="(time split_libraries.py -m $mapf -f $fast -l $splitlib_min -L $splitlib_max -b $barco -M 25 -k -e 0 -H 20 -n 1 -o splits_dir/) > $time_pwd/time_split_libraries.log 2>&1"; }
	#3} else {	$cmd="split_libraries.py -m $mapf -f $fast -l 200 -L 550 -b $barco -M 25 -k -e 0 -H 20 -n 1 -o splits_dir/ 2>&1 | tee -a $err";
	#	$cmd="split_libraries.py -m $mapf -f $fast -l 200 -L 550 -b $barco -M 25 -k -e 0 -H 20 -n 1 -o splits_dir/ 2>&1 | tee -a $err";
	#}
	system($cmd);
	if ($? == -1) {print "command failed: $!\n"; exit; }
	print LOG $cmd."\n";
	#Record sequence information for log_seq_info.txt
	for ($a = 0; $a <= $#sample_names; $a++) {
		$cmd = `grep '^$sample_names[$a]\\s' splits_dir/split_library_log.txt | tail -n 1`;
		$cmd=~/$sample_names[$a]\t(\d+).*?/;
		push @post_splitlib, $1;
	}
}
print SEQ "post_splitlib\t".join("\t", @post_splitlib)."\n";
###############################################################
#Run Modularized Code
###############################################################
#Chimera checking, if user asked for it
if ($chimeraChecking eq "y") {
	if ($gg eq "all") {
		&chimera_check("gg2011"); &chimera_check("gg2013"); &chimera_check("silva111");
	}
	elsif ($clust ne "uparse") { #uparse will do its own chimera checking within OTU picking
		&chimera_check($gg);
	}
} else {
	if ($gg eq "all") {
		&no_chimera_check("gg2011"); &no_chimera_check("gg2013"); &no_chimera_check("silva111");
	} else {
		&no_chimera_check($gg);
	}
}
#OTU picking algorithm
if ($clust eq "abundantotu") {
	if ($gg eq "all") {
		&abundantOTU("gg2011"); &abundantOTU("gg2013"); &abundantOTU("silva111");
	} else {
		&abundantOTU($gg);
	}
} elsif ($clust eq "uclust") {
	if ($gg eq "all") {
		&uclust("gg2011"); &uclust("gg2013"); &uclust("silva111");
	} else {
		&uclust($gg);
	}
} elsif ($clust eq "cdhit") {
	if ($gg eq "all") {
		&cdhit("gg2011"); &cdhit("gg2013"); &cdhit("silva111");
	} else {
		&cdhit($gg);
	}
} elsif ($clust eq "dnaclust") {
	if ($gg eq "all") {
		&dnaclust("gg2011"); &dnaclust("gg2013"); &dnaclust("silva111");
	} else {
		&dnaclust($gg);
	}
} elsif ($clust eq "mothur") {
	if ($linkage eq "all" && $gg eq "all") {
		&qiime_mothur("average", "gg2011"); &qiime_mothur("average", "gg2013"); &qiime_mothur("average", "silva111");
		&qiime_mothur("nearest", "gg2011"); &qiime_mothur("nearest", "gg2013"); &qiime_mothur("nearest", "silva111");
		&qiime_mothur("furthest", "gg2011"); &qiime_mothur("furthest", "gg2013"); &qiime_mothur("furthest", "silva111");
	} elsif ($linkage eq "all") {
		&qiime_mothur("average", $gg); &qiime_mothur("nearest", $gg); &qiime_mothur("furthest", $gg);
	} elsif ($gg eq "all") {
		&qiime_mothur($linkage, "gg2011"); &qiime_mothur($linkage, "gg2013"); &qiime_mothur($linkage, "silva111");
	} else {
		&qiime_mothur($linkage, $gg);
	}
} elsif ($clust eq "uclust-ref") {
	if ($gg eq "all") {
		&uclust_ref("gg2011"); &uclust_ref("gg2013"); &uclust_ref("silva111");
	} else {
		&uclust_ref($gg);
	}
} elsif ($clust eq "uclust-ref-strict") {
	if ($gg eq "all") {
		&uclust_ref_strict("gg2011"); &uclust_ref_strict("gg2013"); &uclust_ref_strict("silva111");
	} else {
		&uclust_ref_strict($gg);
	}
} elsif ($clust eq "blast") {
	if ($gg eq "all") {
		&blast("gg2011"); &blast("gg2013"); &blast("silva111");
	} else {
		&blast($gg);
	}
} elsif ($clust eq "uparse") {
	&uparse($gg);
} elsif ($clust eq "swarm") {
	&swarm($gg);
} elsif ($clust eq "usearch61") {
	&usearch61($gg);
}elsif ($clust eq "all") {
	if ($gg eq "all") {
		&abundantOTU("gg2011"); &abundantOTU("gg2013"); &abundantOTU("silva111");
		&uclust("gg2011"); &uclust("gg2013"); &uclust("silva111");
		&cdhit("gg2011"); &cdhit("gg2013"); &cdhit("silva111");
		&dnaclust("gg2011"); &dnaclust("gg2013"); &dnaclust("silva111");
		if ($linkage eq "all") {
			#&qiime_mothur("average", "gg2011"); &qiime_mothur("average", "gg2013"); &qiime_mothur("average", "silva111");
                        #&qiime_mothur("nearest", "gg2011"); &qiime_mothur("nearest", "gg2013"); &qiime_mothur("nearest", "silva111");
                        #&qiime_mothur("furthest", "gg2011"); &qiime_mothur("furthest", "gg2013"); &qiime_mothur("furthest", "silva111");
		} else {
			#&qiime_mothur($linkage, "gg2011"); &qiime_mothur($linkage, "gg2013"); &qiime_mothur($linkage, "silva111");
		}
		&uclust_ref("gg2011"); &uclust_ref("gg2013"); &uclust_ref("silva111");
		&blast("gg2011"); &blast("gg2013"); &blast("silva111");
		&uparse("gg2011"); &uparse("gg2013"); &uparse("silva111");
	} else {
		&abundantOTU($gg);
		&uclust($gg);
		&cdhit($gg);
		&dnaclust($gg);
		if ($linkage eq "all") {
			#&qiime_mothur("average", $gg); &qiime_mothur("nearest", $gg); &qiime_mothur("furthest", $gg);
		} else {
			#&qiime_mothur($linkage, $gg);
		}
		&uclust_ref($gg);
		&blast($gg);
		&uparse($gg);
	}
}
#Adjust $clust variable if mothur is used
if ($clust eq "mothur") {
	$clust = $clust."_".$linkage;
}
#Align Sequences
#if ($clust eq "all") {
#	&align_seqs("abundantotu", $gg);
#	&align_seqs("uclust", $gg);
#	&align_seqs("cdhit", $gg);
#	&align_seqs("dnaclust", $gg);
#	#&align_seqs("mothur_average", $gg);
#	#&align_seqs("mothur_nearest", $gg);
#	#&align_seqs("mothur_furthest", $gg);
#	&align_seqs("uclust-ref", $gg);
#	&align_seqs("blast", $gg);
#	&align_seqs("uparse", $gg);
#} else {
#	&align_seqs($clust, $gg);
#}
#Taxonomic assignment
if ($taxon eq "blast") {
	if (($clust eq "all") && ($gg eq "all")) {
		&taxa_blast("abundantotu", "gg2011"); &taxa_blast("abundantotu", "gg2013"); &taxa_blast("abundantotu", "silva111");
 		&taxa_blast("uclust", "gg2011"); &taxa_blast("uclust", "gg2013"); &taxa_blast("uclust", "silva111");
		&taxa_blast("cdhit", "gg2011"); &taxa_blast("cdhit", "gg2013"); &taxa_blast("cdhit", "silva111");
		&taxa_blast("dnaclust", "gg2011"); &taxa_blast("dnaclust", "gg2013"); &taxa_blast("dnaclust", "silva111");
		#&taxa_blast("mothur_average", "gg2011"); &taxa_blast("mothur_average", "gg2013"); &taxa_blast("mothur_average", "silva111");
		#&taxa_blast("mothur_nearest", "gg2011"); &taxa_blast("mothur_nearest", "gg2013"); &taxa_blast("mothur_nearest", "silva111");
		#&taxa_blast("mothur_furthest", "gg2011"); &taxa_blast("mothur_furthest", "gg2013"); &taxa_blast("mothur_furthest", "silva111");
		&taxa_blast("uclust-ref", "gg2011"); &taxa_blast("uclust-ref", "gg2013"); &taxa_blast("uclust-ref", "silva111");
		&taxa_blast("uclust-ref-strict", "gg2011"); &taxa_blast("uclust-ref-strict", "gg2013"); &taxa_blast("uclust-ref-strict", "silva111");
		&taxa_blast("blast", "gg2011"); &taxa_blast("blast", "gg2013"); &taxa_blast("blast", "silva111");
		&taxa_blast("uparse", "gg2011"); &taxa_blast("uparse", "gg2013"); &taxa_blast("uparse", "silva111");
	} elsif ($gg eq "all") {
		&taxa_blast($clust, "gg2011");
		&taxa_blast($clust, "gg2013");
		&taxa_blast($clust, "silva111");
	} elsif ($clust eq "all") {
		&taxa_blast("abundantotu", $gg);
		&taxa_blast("uclust", $gg);
		&taxa_blast("cdhit", $gg);
		&taxa_blast("dnaclust", $gg);
		#&taxa_blast("mothur_average", $gg);
		#&taxa_blast("mothur_nearest", $gg);
		#&taxa_blast("mothur_furthest", $gg);
		&taxa_blast("uclust-ref", $gg);
		&taxa_blast("uclust-ref-strict", $gg);
		&taxa_blast("blast", $gg);
		&taxa_blast("uparse", $gg);
	} else {
		&taxa_blast($clust, $gg);
	}
} elsif ($taxon eq "rdp-training") {
	if (($clust eq "all") && ($gg eq "all")) {
		&taxa_rdp_training("abundantotu", "gg2011"); &taxa_rdp_training("abundantotu", "gg2013"); &taxa_rdp_training("abundantotu", "silva111");
		&taxa_rdp_training("uclust", "gg2011"); &taxa_rdp_training("uclust", "gg2013"); &taxa_rdp_training("uclust", "silva111");
		&taxa_rdp_training("cdhit", "gg2011"); &taxa_rdp_training("cdhit", "gg2013"); &taxa_rdp_training("cdhit", "silva111");
		&taxa_rdp_training("dnaclust", "gg2011"); &taxa_rdp_training("dnaclust", "gg2013"); &taxa_rdp_training("dnaclust", "silva111");
		#&taxa_rdp_training("mothur_average", "gg2011"); &taxa_rdp_training("mothur_average", "gg2013"); &taxa_rdp_training("mothur_average", "silva111");
		#&taxa_rdp_training("mothur_nearest", "gg2011"); &taxa_rdp_training("mothur_nearest", "gg2013"); &taxa_rdp_training("mothur_nearest", "silva111");
		#&taxa_rdp_training("mothur_furthest", "gg2011"); &taxa_rdp_training("mothur_furthest", "gg2013"); &taxa_rdp_training("mothur_furthest", "silva111");
		&taxa_rdp_training("uclust-ref", "gg2011"); &taxa_rdp_training("uclust-ref", "gg2013"); &taxa_rdp_training("uclust-ref", "silva111");
		&taxa_rdp_training("uclust-ref-strict", "gg2011"); &taxa_rdp_training("uclust-ref-strict", "gg2013"); &taxa_rdp_training("uclust-ref-strict", "silva111");
		&taxa_rdp_training("blast", "gg2011"); &taxa_rdp_training("blast", "gg2013"); &taxa_rdp_training("blast", "silva111");
		&taxa_rdp_training("uparse", "gg2011"); &taxa_rdp_training("uparse", "gg2013"); &taxa_rdp_training("uparse", "silva111");
	} elsif ($gg eq "all") {
		&taxa_rdp_training($clust, "gg2011");
		&taxa_rdp_training($clust, "gg2013");
		&taxa_rdp_training($clust, "silva111");
	} elsif ($clust eq "all") {
                &taxa_rdp_training("abundantotu", $gg);
                &taxa_rdp_training("uclust", $gg);
                &taxa_rdp_training("cdhit", $gg);
                &taxa_rdp_training("dnaclust", $gg);
                #&taxa_rdp_training("mothur_average", $gg);
                #&taxa_rdp_training("mothur_nearest", $gg);
                #&taxa_rdp_training("mothur_furthest", $gg);
                &taxa_rdp_training("uclust-ref", $gg);
		&taxa_rdp_training("uclust-ref-strict", $gg);
                &taxa_rdp_training("blast", $gg);
                &taxa_rdp_training("uparse", $gg);
	} else {
		&taxa_rdp_training($clust, $gg);
	}
} elsif ($taxon eq "all") {
	if (($clust eq "all") && ($gg eq "all")) {
		&taxa_blast("abundantotu", "gg2011"); &taxa_blast("abundantotu", "gg2013"); &taxa_blast("abundantotu", "silva111");
                &taxa_blast("uclust", "gg2011"); &taxa_blast("uclust", "gg2013"); &taxa_blast("uclust", "silva111");
                &taxa_blast("cdhit", "gg2011"); &taxa_blast("cdhit", "gg2013"); &taxa_blast("cdhit", "silva111");
                &taxa_blast("dnaclust", "gg2011"); &taxa_blast("dnaclust", "gg2013"); &taxa_blast("dnaclust", "silva111");
                #&taxa_blast("mothur_average", "gg2011"); &taxa_blast("mothur_average", "gg2013"); &taxa_blast("mothur_average", "silva111");
                #&taxa_blast("mothur_nearest", "gg2011"); &taxa_blast("mothur_nearest", "gg2013"); &taxa_blast("mothur_nearest", "silva111");
                #&taxa_blast("mothur_furthest", "gg2011"); &taxa_blast("mothur_furthest", "gg2013"); &taxa_blast("mothur_furthest", "silva111");
                &taxa_blast("uclust-ref", "gg2011"); &taxa_blast("uclust-ref", "gg2013"); &taxa_blast("uclust-ref", "silva111");
		&taxa_blast("uclust-ref-strict", "gg2011"); &taxa_blast("uclust-ref-strict", "gg2013"); &taxa_blast("uclust-ref-strict", "silva111");
                &taxa_blast("blast", "gg2011"); &taxa_blast("blast", "gg2013"); &taxa_blast("blast", "silva111");
                &taxa_blast("uparse", "gg2011"); &taxa_blast("uparse", "gg2013"); &taxa_blast("uparse", "silva111");
		&taxa_rdp_training("abundantotu", "gg2011"); &taxa_rdp_training("abundantotu", "gg2013"); &taxa_rdp_training("abundantotu", "silva111");
                &taxa_rdp_training("uclust", "gg2011"); &taxa_rdp_training("uclust", "gg2013"); &taxa_rdp_training("uclust", "silva111");
                &taxa_rdp_training("cdhit", "gg2011"); &taxa_rdp_training("cdhit", "gg2013"); &taxa_rdp_training("cdhit", "silva111");
                &taxa_rdp_training("dnaclust", "gg2011"); &taxa_rdp_training("dnaclust", "gg2013"); &taxa_rdp_training("dnaclust", "silva111");
                #&taxa_rdp_training("mothur_average", "gg2011"); &taxa_rdp_training("mothur_average", "gg2013"); &taxa_rdp_training("mothur_average", "silva111");
                #&taxa_rdp_training("mothur_nearest", "gg2011"); &taxa_rdp_training("mothur_nearest", "gg2013"); &taxa_rdp_training("mothur_nearest", "silva111");
                #&taxa_rdp_training("mothur_furthest", "gg2011"); &taxa_rdp_training("mothur_furthest", "gg2013"); &taxa_rdp_training("mothur_furthest", "silva111");
                &taxa_rdp_training("uclust-ref", "gg2011"); &taxa_rdp_training("uclust-ref", "gg2013"); &taxa_rdp_training("uclust-ref", "silva111");
		&taxa_rdp_training("uclust-ref-strict", "gg2011"); &taxa_rdp_training("uclust-ref-strict", "gg2013"); &taxa_rdp_training("uclust-ref-strict", "silva111");
                &taxa_rdp_training("blast", "gg2011"); &taxa_rdp_training("blast", "gg2013"); &taxa_rdp_training("blast", "silva111");
                &taxa_rdp_training("uparse", "gg2011"); &taxa_rdp_training("uparse", "gg2013"); &taxa_rdp_training("uparse", "silva111");
	} elsif ($gg eq "all") {
		&taxa_blast($clust, "gg2011");
               	&taxa_blast($clust, "gg2013");
               	&taxa_blast($clust, "silva111");
		&taxa_rdp_training($clust, "gg2011");
               	&taxa_rdp_training($clust, "gg2013");
               	&taxa_rdp_training($clust, "silva111");
	} elsif ($clust eq "all") {
                &taxa_blast("abundantotu", $gg);
		&taxa_rdp_training("abundantotu", $gg);
		&taxa_blast("uclust", $gg);
		&taxa_rdp_training("uclust", $gg);
		&taxa_blast("cdhit", $gg);
		&taxa_rdp_training("cdhit", $gg);
		&taxa_blast("dnaclust", $gg);
		&taxa_rdp_training("dnaclust", $gg);
		#&taxa_blast("mothur_linkage", $gg);
		#&taxa_rdp_training("mothur_linkage", $gg);
		#&taxa_blast("mothur_nearest", $gg);
		#&taxa_rdp_training("mothur_nearest", $gg);
		#&taxa_blast("mothur_furthest", $gg);
		#&taxa_rdp_training("mothur_furthest", $gg);
		&taxa_blast("uclust-ref", $gg);
		&taxa_rdp_training("uclust-ref", $gg);
		&taxa_blast("uclust-ref-strict", $gg);
		&taxa_rdp_training("uclust-ref-strict", $gg);
		&taxa_blast("blast", $gg);
		&taxa_rdp_training("blast", $gg);
                #$clust = "uparse";
		&taxa_blast("uparse", $gg);
		&taxa_rdp_training("uparse", $gg);
	} else {
		&taxa_blast($clust, $gg);
		&taxa_rdp_training($clust, $gg);
	}
}
#Taxon assignment independent analysis
if (($clust eq "all") && ($gg eq "all")) {
	&align_seqs("abundantotu", "gg2011"); &align_seqs("abundantotu", "gg2013"); &align_seqs("abundantotu", "silva111");
	&filter_alignment("abundantotu", "gg2011"); &filter_alignment("abundantotu", "gg2013"); &filter_alignment("abundantotu", "silva111");
	&make_phylogeny("abundantotu", "gg2011"); &make_phylogeny("abundantotu", "gg2013"); &make_phylogeny("abundantotu", "silva111");
	&make_ggpruned("abundantotu", "gg2011"); &make_ggpruned("abundantotu", "gg2013"); &make_ggpruned("abundantotu", "silva111");
	
	&align_seqs("uclust", "gg2011"); &align_seqs("uclust", "gg2013"); &align_seqs("uclust", "silva111");
	&filter_alignment("uclust", "gg2011"); &filter_alignment("uclust", "gg2013"); &filter_alignment("uclust", "silva111");
	&make_phylogeny("uclust", "gg2011"); &make_phylogeny("uclust", "gg2013"); &make_phylogeny("uclust", "silva111");
	&make_ggpruned("uclust", "gg2011"); &make_ggpruned("uclust", "gg2013"); &make_ggpruned("uclust", "silva111");

	&align_seqs("cdhit", "gg2011"); &align_seqs("cdhit", "gg2013"); &align_seqs("cdhit", "silva111");
	&filter_alignment("cdhit", "gg2011"); &filter_alignment("cdhit", "gg2013"); &filter_alignment("cdhit", "silva111");
	&make_phylogeny("cdhit", "gg2011"); &make_phylogeny("cdhit", "gg2013"); &make_phylogeny("cdhit", "silva111");
	&make_ggpruned("cdhit", "gg2011"); &make_ggpruned("cdhit", "gg2013"); &make_ggpruned("cdhit", "silva111");
	
	&align_seqs("dnaclust", "gg2011"); &align_seqs("dnaclust", "gg2013"); &align_seqs("dnaclust", "silva111");
	&filter_alignment("dnaclust", "gg2011"); &filter_alignment("dnaclust", "gg2013"); &filter_alignment("dnaclust", "silva111");
	&make_phylogeny("dnaclust", "gg2011"); &make_phylogeny("dnaclust", "gg2013"); &make_phylogeny("dnaclust", "silva111");
	&make_ggpruned("dnaclust", "gg2011"); &make_ggpruned("dnaclust", "gg2013"); &make_ggpruned("dnaclust", "silva111");
	
	#&align_seqs("mothur_average", "gg2011"); &align_seqs("mothur_average", "gg2013"); &align_seqs("mothur_average", "silva111");
	#&filter_alignment("mothur_average", "gg2011"); &filter_alignment("mothur_average", "gg2013"); &filter_alignment("mothur_average", "silva111");
	#&make_phylogeny("mothur_average", "gg2011"); &make_phylogeny("mothur_average", "gg2013"); &make_phylogeny("mothur_average", "silva111");
	#&make_ggpruned("mothur_average", "gg2011"); &make_ggpruned("mothur_average", "gg2013"); &make_ggpruned("mothur_average", "silva111");
	
	#&align_seqs("mothur_nearest", "gg2011"); &align_seqs("mothur_nearest", "gg2013"); &align_seqs("mothur_nearest", "silva111");
	#&filter_alignment("mothur_nearest", "gg2011"); &filter_alignment("mothur_nearest", "gg2013"); &filter_alignment("mothur_nearest", "silva111");
	#&make_phylogeny("mothur_nearest", "gg2011"); &make_phylogeny("mothur_nearest", "gg2013"); &make_phylogeny("mothur_nearest", "silva111");
	#&make_ggpruned("mothur_nearest", "gg2011"); &make_ggpruned("mothur_nearest", "gg2013"); &make_ggpruned("mothur_nearest", "silva111");
	
	#&align_seqs("mothur_furthest", "gg2011"); &align_seqs("mothur_furthest", "gg2013"); &align_seqs("mothur_furthest", "silva111");
	#&filter_alignment("mothur_furthest", "gg2011"); &filter_alignment("mothur_furthest", "gg2013"); &filter_alignment("mothur_furthest", "silva111");
	#&make_phylogeny("mothur_furthest", "gg2011"); &make_phylogeny("mothur_furthest", "gg2013"); &make_phylogeny("mothur_furthest", "silva111");
	#&make_ggpruned("mothur_nearest", "gg2011"); &make_ggpruned("mothur_furthest", "gg2013"); &make_ggpruned("mothur_furthest", "silva111");

	
	&align_seqs("uclust-ref", "gg2011"); &align_seqs("uclust-ref", "gg2013"); &align_seqs("uclust-ref", "silva111");
	&filter_alignment("uclust-ref", "gg2011"); &filter_alignment("uclust-ref", "gg2013"); &filter_alignment("uclust-ref", "silva111");
	&make_phylogeny("uclust-ref", "gg2011"); &make_phylogeny("uclust-ref", "gg2013"); &make_phylogeny("uclust-ref", "silva111");
	&make_ggpruned("uclust-ref", "gg2011"); &make_ggpruned("uclust-ref", "gg2013"); &make_ggpruned("uclust-ref", "silva111");
	
	&align_seqs("uclust-ref-strict", "gg2011"); &align_seqs("uclust-ref-strict", "gg2013"); &align_seqs("uclust-ref-strict", "silva111");
	&filter_alignment("uclust-ref-strict", "gg2011"); &filter_alignment("uclust-ref-strict", "gg2013"); &filter_alignment("uclust-ref-strict", "silva111");
	&make_phylogeny("uclust-ref-strict", "gg2011"); &make_phylogeny("uclust-ref-strict", "gg2013"); &make_phylogeny("uclust-ref-strict", "silva111");
	&make_ggpruned("uclust-ref-strict", "gg2011"); &make_ggpruned("uclust-ref-strict", "gg2013"); &make_ggpruned("uclust-ref-strict", "silva111");

	&align_seqs("blast", "gg2011"); &align_seqs("blast", "gg2013"); &align_seqs("blast", "silva111");
	&filter_alignment("blast", "gg2011"); &filter_alignment("blast", "gg2013"); &filter_alignment("blast", "silva111");
	&make_phylogeny("blast", "gg2011"); &make_phylogeny("blast", "gg2013"); &make_phylogeny("blast", "silva111");
	&make_ggpruned("blast", "gg2011"); &make_ggpruned("blast", "gg2013"); &make_ggpruned("blast", "silva111");
	
	&align_seqs("uparse", "gg2011"); &align_seqs("uparse", "gg2013"); &align_seqs("uparse", "silva111");
	&filter_alignment("uparse", "gg2011"); &filter_alignment("uparse", "gg2013"); &filter_alignment("uparse", "silva111");
	&make_phylogeny("uparse", "gg2011"); &make_phylogeny("uparse", "gg2013"); &make_phylogeny("uparse", "silva111");
	&make_ggpruned("uparse", "gg2011"); &make_ggpruned("uparse", "gg2013"); &make_ggpruned("uparse", "silva111");
} elsif ($gg eq "all") {
	&align_seqs($clust, "gg2011");
	&filter_alignment($clust, "gg2011");
	&make_phylogeny($clust, "gg2011");
	&make_ggpruned($clust, "gg2011");

	&align_seqs($clust, "gg2013");
	&filter_alignment($clust, "gg2013");
	&make_phylogeny($clust, "gg2013");
	&make_ggpruned($clust, "gg2013");

	&align_seqs($clust, "silva111");
	&filter_alignment($clust, "silva111");
	&make_phylogeny($clust, "silva111");
	&make_ggpruned($clust, "silva111");
} elsif ($clust eq "all") {
	&align_seqs("abundantotu", $gg);
	&filter_alignment("abundantotu", $gg);
	&make_phylogeny("abundantotu", $gg);
	&make_ggpruned("abundantotu", $gg);

	&align_seqs("uclust", $gg);
	&filter_alignment("uclust", $gg);
	&make_phylogeny("uclust", $gg);
	&make_ggpruned("uclust", $gg);

	&align_seqs("cdhit", $gg);
	&filter_alignment("cdhit", $gg);
	&make_phylogeny("cdhit", $gg);
	&make_ggpruned("cdhit", $gg);	

	&align_seqs("dnaclust", $gg);
	&filter_alignment("dnaclust", $gg);
	&make_phylogeny("dnaclust", $gg);
	&make_ggpruned("dnaclust", $gg);
	
	#&align_seqs("mothur_average", $gg);
	#&filter_alignment("mothur_average", $gg);
	#&make_phylogeny("mothur_average", $gg);
	#&make_ggpruned("mothur_average", $gg);

	#&align_seqs("mothur_nearest", $gg);
	#&filter_alignment("mothur_nearest", $gg);
	#&make_phylogeny("mothur_nearest", $gg);
	#&make_ggpruned("mothur_nearest", $gg);

	#&align_seqs("mothur_furthest", $gg);
	#&filter_alignment("mothur_furthest", $gg);
	#&make_phylogeny("mothur_furthest", $gg);
	#&make_ggpruned("mothur_furthest", $gg);

 	&align_seqs("uclust-ref", $gg);
	&filter_alignment("uclust-ref", $gg);
	&make_phylogeny("uclust-ref", $gg);
	&make_ggpruned("uclust-ref", $gg);

	&align_seqs("uclust-ref-strict", $gg);
	&filter_alignment("uclust-ref-strict", $gg);
	&make_phylogeny("uclust-ref-strict", $gg);
	&make_ggpruned("uclust-ref-strict", $gg);

	&align_seqs("blast", $gg);
	&filter_alignment("blast", $gg);
	&make_phylogeny("blast", $gg);
	&make_ggpruned("blast", $gg);

	&align_seqs("uparse", $gg);
	&filter_alignment("uparse", $gg);
	&make_phylogeny("uparse", $gg);
	&make_ggpruned("uparse", $gg);
} else {
        &align_seqs($clust, $gg);
	&filter_alignment($clust, $gg);
	&make_phylogeny($clust, $gg);
	&make_ggpruned($clust, $gg);
}
#Analysis
if (($clust eq "all") && ($taxon eq "all") && ($gg eq "all")) {
	&make_otu_table("abundantotu", "blast", "gg2011"); &make_otu_table("abundantotu", "blast", "gg2013"); &make_otu_table("abundantotu", "blast", "silva111");
        &filter_otus("abundantotu", "blast", "gg2011"); &filter_otus("abundantotu", "blast", "gg2013"); &filter_otus("abundantotu", "blast", "silva111");
        &removeRoot("abundantotu", "blast", "gg2011"); &removeRoot("abundantotu", "blast", "gg2013"); &removeRoot("abundantotu", "blast", "silva111");
        &summarize_taxa("abundantotu", "blast", "gg2011"); &summarize_taxa("abundantotu", "blast", "gg2013"); &summarize_taxa("abundantotu", "blast", "silva111");
        &alpha_div("abundantotu", "blast", "gg2011"); &alpha_div("abundantotu", "blast", "gg2013"); &alpha_div("abundantotu", "blast", "silva111");
        &beta_div("abundantotu", "blast", "gg2011"); &beta_div("abundantotu", "blast", "gg2013"); &beta_div("abundantotu", "blast", "silva111");
	&RAnalysis("abundantotu", "blast", "gg2011"); &RAnalysis("abundantotu", "blast", "gg2013"); &RAnalysis("abundantotu", "blast", "silva111");
	&cleanUp("abundantotu", "blast", "gg2011"); &cleanUp("abundantotu", "blast", "gg2013"); &cleanUp("abundantotu", "blast", "silva111");
        &make_otu_table("abundantotu", "rdp-training", "gg2011"); &make_otu_table("abundantotu", "rdp-training", "gg2013"); &make_otu_table("abundantotu", "rdp-training", "silva111");
        &filter_otus("abundantotu", "rdp-training", "gg2011"); &filter_otus("abundantotu", "rdp-training", "gg2013"); &filter_otus("abundantotu", "rdp-training", "silva111");
        &removeRoot("abundantotu", "rdp-training", "gg2011"); &removeRoot("abundantotu", "rdp-training", "gg2013"); &removeRoot("abundantotu", "rdp-training", "silva111");
        &summarize_taxa("abundantotu", "rdp-training", "gg2011"); &summarize_taxa("abundantotu", "rdp-training", "gg2013"); &summarize_taxa("abundantotu", "rdp-training", "silva111");
        &alpha_div("abundantotu", "rdp-training", "gg2011"); &alpha_div("abundantotu", "rdp-training", "gg2013"); &alpha_div("abundantotu", "rdp-training", "silva111");
        &beta_div("abundantotu", "rdp-training", "gg2011"); &beta_div("abundantotu", "rdp-training", "gg2013"); &beta_div("abundantotu", "rdp-training", "silva111");
	&RAnalysis("abundantotu", "rdp-training", "gg2011"); &RAnalysis("abundantotu", "rdp-training", "gg2013"); &RAnalysis("abundantotu", "rdp-training", "silva111");
	&cleanUp("abundantotu", "rdp-training", "gg2011"); &cleanUp("abundantotu", "rdp-training", "gg2013"); &cleanUp("abundantotu", "rdp-training", "silva111");

        &make_otu_table("uclust", "blast", "gg2011"); &make_otu_table("uclust", "blast", "gg2013"); &make_otu_table("uclust", "blast", "silva111");
        &filter_otus("uclust", "blast", "gg2011"); &filter_otus("uclust", "blast", "gg2013"); &filter_otus("uclust", "blast", "silva111");
        &removeRoot("uclust", "blast", "gg2011"); &removeRoot("uclust", "blast", "gg2013"); &removeRoot("uclust", "blast", "silva111");
        &summarize_taxa("uclust", "blast", "gg2011"); &summarize_taxa("uclust", "blast", "gg2013"); &summarize_taxa("uclust", "blast", "silva111");
        &alpha_div("uclust", "blast", "gg2011"); &alpha_div("uclust", "blast", "gg2013"); &alpha_div("uclust", "blast", "silva111");
        &beta_div("uclust", "blast", "gg2011"); &beta_div("uclust", "blast", "gg2013"); &beta_div("uclust", "blast", "silva111");
	&RAnalysis("uclust", "blast", "gg2011"); &RAnalysis("uclust", "blast", "gg2013"); &RAnalysis("uclust", "blast", "silva111");
	&cleanUp("uclust", "blast", "gg2011"); &cleanUp("uclust", "blast", "gg2013"); &cleanUp("uclust", "blast", "silva111");
        &make_otu_table("uclust", "rdp-training", "gg2011"); &make_otu_table("uclust", "rdp-training", "gg2013"); &make_otu_table("uclust", "rdp-training", "silva111");
        &filter_otus("uclust", "rdp-training", "gg2011"); &filter_otus("uclust", "rdp-training", "gg2013"); &filter_otus("uclust", "rdp-training", "silva111");
        &removeRoot("uclust", "rdp-training", "gg2011"); &removeRoot("uclust", "rdp-training", "gg2013"); &removeRoot("uclust", "rdp-training", "silva111");
        &summarize_taxa("uclust", "rdp-training", "gg2011"); &summarize_taxa("uclust", "rdp-training", "gg2013"); &summarize_taxa("uclust", "rdp-training", "silva111");
        &alpha_div("uclust", "rdp-training", "gg2011"); &alpha_div("uclust", "rdp-training", "gg2013"); &alpha_div("uclust", "rdp-training", "silva111");
        &beta_div("uclust", "rdp-training", "gg2011"); &beta_div("uclust", "rdp-training", "gg2013"); &beta_div("uclust", "rdp-training", "silva111");
	&RAnalysis("uclust", "rdp-training", "gg2011"); &RAnalysis("uclust", "rdp-training", "gg2013"); &RAnalysis("uclust", "rdp-training", "silva111");
	&cleanUp("uclust", "rdp-training", "gg2011"); &cleanUp("uclust", "rdp-training", "gg2013"); &cleanUp("uclust", "rdp-training", "silva111");

        &make_otu_table("cdhit", "blast", "gg2011"); &make_otu_table("cdhit", "blast", "gg2013"); &make_otu_table("cdhit", "blast", "silva111");
        &filter_otus("cdhit", "blast", "gg2011"); &filter_otus("cdhit", "blast", "gg2013"); &filter_otus("cdhit", "blast", "silva111");
        &removeRoot("cdhit", "blast", "gg2011"); &removeRoot("cdhit", "blast", "gg2013"); &removeRoot("cdhit", "blast", "silva111");
        &summarize_taxa("cdhit", "blast", "gg2011"); &summarize_taxa("cdhit", "blast", "gg2013"); &summarize_taxa("cdhit", "blast", "silva111");
        &alpha_div("cdhit", "blast", "gg2011"); &alpha_div("cdhit", "blast", "gg2013"); &alpha_div("cdhit", "blast", "silva111");
        &beta_div("cdhit", "blast", "gg2011"); &beta_div("cdhit", "blast", "gg2013"); &beta_div("cdhit", "blast", "silva111");
	&RAnalysis("cdhit", "blast", "gg2011"); &RAnalysis("cdhit", "blast", "gg2013"); &RAnalysis("cdhit", "blast", "silva111");
	&cleanUp("cdhit", "blast", "gg2011"); &cleanUp("cdhit", "blast", "gg2013"); &cleanUp("cdhit", "blast", "silva111");
        &make_otu_table("cdhit", "rdp-training", "gg2011"); &make_otu_table("cdhit", "rdp-training", "gg2013"); &make_otu_table("cdhit", "rdp-training", "silva111");
        &filter_otus("cdhit", "rdp-training", "gg2011"); &filter_otus("cdhit", "rdp-training", "gg2013"); &filter_otus("cdhit", "rdp-training", "silva111");
        &removeRoot("cdhit", "rdp-training", "gg2011"); &removeRoot("cdhit", "rdp-training", "gg2013"); &removeRoot("cdhit", "rdp-training", "silva111");
        &summarize_taxa("cdhit", "rdp-training", "gg2011"); &summarize_taxa("cdhit", "rdp-training", "gg2013"); &summarize_taxa("cdhit", "rdp-training", "silva111");
        &alpha_div("cdhit", "rdp-training", "gg2011"); &alpha_div("cdhit", "rdp-training", "gg2013"); &alpha_div("cdhit", "rdp-training", "silva111");
        &beta_div("cdhit", "rdp-training", "gg2011"); &beta_div("cdhit", "rdp-training", "gg2013"); &beta_div("cdhit", "rdp-training", "silva111");
	&RAnalysis("cdhit", "rdp-training", "gg2011"); &RAnalysis("cdhit", "rdp-training", "gg2013"); &RAnalysis("cdhit", "rdp-training", "silva111");
	&cleanUp("cdhit", "rdp-training", "gg2011"); &cleanUp("cdhit", "rdp-training", "gg2013"); &cleanUp("cdhit", "rdp-training", "silva111");

        &make_otu_table("dnaclust", "blast", "gg2011"); &make_otu_table("dnaclust", "blast", "gg2013"); &make_otu_table("dnaclust", "blast", "silva111");
        &filter_otus("dnaclust", "blast", "gg2011"); &filter_otus("dnaclust", "blast", "gg2013"); &filter_otus("dnaclust", "blast", "silva111");
        &removeRoot("dnaclust", "blast", "gg2011"); &removeRoot("dnaclust", "blast", "gg2013"); &removeRoot("dnaclust", "blast", "silva111");
        &summarize_taxa("dnaclust", "blast", "gg2011"); &summarize_taxa("dnaclust", "blast", "gg2013"); &summarize_taxa("dnaclust", "blast", "silva111");
        &alpha_div("dnaclust", "blast", "gg2011"); &alpha_div("dnaclust", "blast", "gg2013"); &alpha_div("dnaclust", "blast", "silva111");
        &beta_div("dnaclust", "blast", "gg2011"); &beta_div("dnaclust", "blast", "gg2013"); &beta_div("dnaclust", "blast", "silva111");
	&RAnalysis("dnaclust", "blast", "gg2011"); &RAnalysis("dnaclust", "blast", "gg2013"); &RAnalysis("dnaclust", "blast", "silva111");
	&cleanUp("dnaclust", "blast", "gg2011"); &cleanUp("dnaclust", "blast", "gg2013"); &cleanUp("dnaclust", "blast", "silva111");
        &make_otu_table("dnaclust", "rdp-training", "gg2011"); &make_otu_table("dnaclust", "rdp-training", "gg2013"); &make_otu_table("dnaclust", "rdp-training", "silva111");
        &filter_otus("dnaclust", "rdp-training", "gg2011"); &filter_otus("dnaclust", "rdp-training", "gg2013"); &filter_otus("dnaclust", "rdp-training", "silva111");
        &removeRoot("dnaclust", "rdp-training", "gg2011"); &removeRoot("dnaclust", "rdp-training", "gg2013"); &removeRoot("dnaclust", "rdp-training", "silva111");
        &summarize_taxa("dnaclust", "rdp-training", "gg2011"); &summarize_taxa("dnaclust", "rdp-training", "gg2013"); &summarize_taxa("dnaclust", "rdp-training", "silva111");
        &alpha_div("dnaclust", "rdp-training", "gg2011"); &alpha_div("dnaclust", "rdp-training", "gg2013"); &alpha_div("dnaclust", "rdp-training", "silva111");
        &beta_div("dnaclust", "rdp-training", "gg2011"); &beta_div("dnaclust", "rdp-training", "gg2013"); &beta_div("dnaclust", "rdp-training", "silva111");
	&RAnalysis("dnaclust", "rdp-training", "gg2011"); &RAnalysis("dnaclust", "rdp-training", "gg2013"); &RAnalysis("dnaclust", "rdp-training", "silva111");
	&cleanUp("dnaclust", "rdp-training", "gg2011"); &cleanUp("dnaclust", "rdp-training", "gg2013"); &cleanUp("dnaclust", "rdp-training", "silva111");

        #&make_otu_table("mothur_average", "blast", "gg2011"); &make_otu_table("mothur_average", "blast", "gg2013"); &make_otu_table("mothur_average", "blast", "silva111");
        #&filter_otus("mothur_average", "blast", "gg2011"); &filter_otus("mothur_average", "blast", "gg2013"); &filter_otus("mothur_average", "blast", "silva111");
        #&removeRoot("mothur_average", "blast", "gg2011"); &removeRoot("mothur_average", "blast", "gg2013"); &removeRoot("mothur_average", "blast", "silva111");
        #&summarize_taxa("mothur_average", "blast", "gg2011"); &summarize_taxa("mothur_average", "blast", "gg2013"); &summarize_taxa("mothur_average", "blast", "silva111");
        #&alpha_div("mothur_average", "blast", "gg2011"); &alpha_div("mothur_average", "blast", "gg2013"); &alpha_div("mothur_average", "blast", "silva111");
        #&beta_div("mothur_average", "blast", "gg2011"); &beta_div("mothur_average", "blast", "gg2013"); &beta_div("mothur_average", "blast", "silva111");
	#&RAnalysis("mothur_average", "blast", "gg2011"); &RAnalysis("mothur_average", "blast", "gg2013"); &RAnalysis("mothur_average", "blast", "silva111");
	#&cleanUp("mothur_average", "blast", "gg2011"); &cleanUp("mothur_average", "blast", "gg2013"); &cleanUp("mothur_average", "blast", "silva111");
        #&make_otu_table("mothur_average", "rdp-training", "gg2011"); &make_otu_table("mothur_average", "rdp-training", "gg2013"); &make_otu_table("mothur_average", "rdp-training", "silva111");
        #&filter_otus("mothur_average", "rdp-training", "gg2011"); &filter_otus("mothur_average", "rdp-training", "gg2013"); &filter_otus("mothur_average", "rdp-training", "silva111");
        #&removeRoot("mothur_average", "rdp-training", "gg2011"); &removeRoot("mothur_average", "rdp-training", "gg2013"); &removeRoot("mothur_average", "rdp-training", "silva111");
        #&summarize_taxa("mothur_average", "rdp-training", "gg2011"); &summarize_taxa("mothur_average", "rdp-training", "gg2013"); &summarize_taxa("mothur_average", "rdp-training", "silva111");
        #&alpha_div("mothur_average", "rdp-training", "gg2011"); &alpha_div("mothur_average", "rdp-training", "gg2013"); &alpha_div("mothur_average", "rdp-training", "silva111");
        #&beta_div("mothur_average", "rdp-training", "gg2011"); &beta_div("mothur_average", "rdp-training", "gg2013"); &beta_div("mothur_average", "rdp-training", "silva111");
	#&RAnalysis("mothur_average", "rdp-training", "gg2011"); &RAnalysis("mothur_average", "rdp-training", "gg2013"); &RAnalysis("mothur_average", "rdp-training", "silva111");
	#&cleanUp("mothur_average", "rdp-training", "gg2011"); &cleanUp("mothur_average", "rdp-training", "gg2013"); &cleanUp("mother_average" , "rdp-training", "silva111");

        #&make_otu_table("mothur_nearest", "blast", "gg2011"); &make_otu_table("mothur_nearest", "blast", "gg2013"); &make_otu_table("mothur_nearest", "blast", "silva111");
        #&filter_otus("mothur_nearest", "blast", "gg2011"); &filter_otus("mothur_nearest", "blast", "gg2013"); &filter_otus("mothur_nearest", "blast", "silva111");
        #&removeRoot("mothur_nearest", "blast", "gg2011"); &removeRoot("mothur_nearest", "blast", "gg2013"); &removeRoot("mothur_nearest", "blast", "silva111");
        #&summarize_taxa("mothur_nearest", "blast", "gg2011"); &summarize_taxa("mothur_nearest", "blast", "gg2013"); &summarize_taxa("mothur_nearest", "blast", "silva111");
        #&alpha_div("mothur_nearest", "blast", "gg2011"); &alpha_div("mothur_nearest", "blast", "gg2013"); &alpha_div("mothur_nearest", "blast", "silva111");
        #&beta_div("mothur_nearest", "blast", "gg2011"); &beta_div("mothur_nearest", "blast", "gg2013"); &beta_div("mothur_nearest", "blast", "silva111");
	#&RAnalysis("mothur_nearest", "blast", "gg2011"); &RAnalysis("mothur_nearest", "blast", "gg2013"); &RAnalysis("mothur_nearest", "blast", "silva111");
	#&cleanUp("mothur_nearest", "blast", "gg2011"); &cleanUp("mothur_nearest", "blast", "gg2013"); &cleanUp("mothur_nearest", "blast", "silva111");
        #&make_otu_table("mothur_nearest", "rdp-training", "gg2011"); &make_otu_table("mothur_nearest", "rdp-training", "gg2013"); &make_otu_table("mothur_nearest", "rdp-training", "silva111");
        #&filter_otus("mothur_nearest", "rdp-training", "gg2011"); &filter_otus("mothur_nearest", "rdp-training", "gg2013"); &filter_otus("mothur_nearest", "rdp-training", "silva111");
        #&removeRoot("mothur_nearest", "rdp-training", "gg2011"); &removeRoot("mothur_nearest", "rdp-training", "gg2013"); &removeRoot("mothur_nearest", "rdp-training", "silva111");
        #&summarize_taxa("mothur_nearest", "rdp-training", "gg2011"); &summarize_taxa("mothur_nearest", "rdp-training", "gg2013"); &summarize_taxa("mothur_nearest", "rdp-training", "silva111");
        #&alpha_div("mothur_nearest", "rdp-training", "gg2011"); &alpha_div("mothur_nearest", "rdp-training", "gg2013"); &alpha_div("mothur_nearest", "rdp-training", "silva111");
        #&beta_div("mothur_nearest", "rdp-training", "gg2011"); &beta_div("mothur_nearest", "rdp-training", "gg2013"); &beta_div("mothur_nearest", "rdp-training", "silva111");
	#&RAnalysis("mothur_nearest", "rdp-training", "gg2011"); &RAnalysis("mothur_nearest", "rdp-training", "gg2013"); &RAnalysis("mothur_nearest", "rdp-training", "silva111");
	#&cleanUp("mothur_nearest", "rdp-training", "gg2011"); &cleanUp("mothur_nearest", "rdp-training", "gg2013"); &cleanUp("mothur_nearest", "rdp-training", "silva111");

        #&make_otu_table("mothur_furthest", "blast", "gg2011"); &make_otu_table("mothur_furthest", "blast", "gg2013"); &make_otu_table("mothur_furthest", "blast", "silva111");
        #&filter_otus("mothur_furthest", "blast", "gg2011"); &filter_otus("mothur_furthest", "blast", "gg2013"); &filter_otus("mothur_furthest", "blast", "silva111");
        #&removeRoot("mothur_furthest", "blast", "gg2011"); &removeRoot("mothur_furthest", "blast", "gg2013"); &removeRoot("mothur_furthest", "blast", "silva111");
        #&summarize_taxa("mothur_furthest", "blast", "gg2011"); &summarize_taxa("mothur_furthest", "blast", "gg2013"); &summarize_taxa("mothur_furthest", "blast", "silva111");
        #&alpha_div("mothur_furthest", "blast", "gg2011"); &alpha_div("mothur_furthest", "blast", "gg2013"); &alpha_div("mothur_furthest", "blast", "silva111");
        #&beta_div("mothur_furthest", "blast", "gg2011"); &beta_div("mothur_furthest", "blast", "gg2013"); &beta_div("mothur_furthest", "blast", "silva111");
	#&RAnalysis("mothur_furthest", "blast", "gg2011"); &RAnalysis("mothur_furthest", "blast", "gg2013"); &RAnalysis("mothur_furthest", "blast", "silva111");
	#&cleanUp("mothur_furthest", "blast", "gg2011"); &cleanUp("mothur_furthest", "blast", "gg2013"); &cleanUp("mothur_furthest", "blast", "silva111");
        #&make_otu_table("mothur_furthest", "rdp-training", "gg2011"); &make_otu_table("mothur_furthest", "rdp-training", "gg2013"); &make_otu_table("mothur_nearest", "rdp-training", "silva111");
        #&filter_otus("mothur_furthest", "rdp-training", "gg2011"); &filter_otus("mothur_furthest", "rdp-training", "gg2013"); &filter_otus("mothur_furthest", "rdp-training", "silva111");
        #&removeRoot("mothur_furthest", "rdp-training", "gg2011"); &removeRoot("mothur_furthest", "rdp-training", "gg2013"); &removeRoot("mothur_furthest", "rdp-training", "silva111");
        #&summarize_taxa("mothur_furthest", "rdp-training", "gg2011"); &summarize_taxa("mothur_furthest", "rdp-training", "gg2013"); &summarize_taxa("mothur_furthest", "rdp-training", "silva111");
        #&alpha_div("mothur_furthest", "rdp-training", "gg2011"); &alpha_div("mothur_furthest", "rdp-training", "gg2013"); &alpha_div("mothur_furthest", "rdp-training", "silva111");
        #&beta_div("mothur_furthest", "rdp-training", "gg2011"); &beta_div("mothur_furthest", "rdp-training", "gg2013"); &beta_div("mothur_furthest", "rdp-training", "silva111");
	#&RAnalysis("mothur_furthest", "rdp-training", "gg2011"); &RAnalysis("mothur_furthest", "rdp-training", "gg2013"); &RAnalysis("mothur_furthest", "rdp-training", "silva111");
	#&cleanUp("mothur_furthest", "rdp-training, "gg2011"); &cleanUp("mothur_furthest", "rdp-training", "gg2013"); &cleanUp("mothur_furthest", "rdp-training", "silva111");


	&make_otu_table("uclust-ref", "blast", "gg2011"); &make_otu_table("uclust-ref", "blast", "gg2013"); &make_otu_table("uclust-ref", "blast", "silva111");
        &filter_otus("uclust-ref", "blast", "gg2011"); &filter_otus("uclust-ref", "blast", "gg2013"); &filter_otus("uclust-ref", "blast", "silva111");
        &removeRoot("uclust-ref", "blast", "gg2011"); &removeRoot("uclust-ref", "blast", "gg2013"); &removeRoot("uclust-ref", "blast", "silva111");
        &summarize_taxa("uclust-ref", "blast", "gg2011"); &summarize_taxa("uclust-ref", "blast", "gg2013"); &summarize_taxa("uclust-ref", "blast", "silva111");
        &alpha_div("uclust-ref", "blast", "gg2011"); &alpha_div("uclust-ref", "blast", "gg2013"); &alpha_div("uclust-ref", "blast", "silva111");
        &beta_div("uclust-ref", "blast", "gg2011"); &beta_div("uclust-ref", "blast", "gg2013"); &beta_div("uclust-ref", "blast", "silva111");
	&RAnalysis("uclust-ref", "blast", "gg2011"); &RAnalysis("uclust-ref", "blast", "gg2013"); &RAnalysis("uclust-ref", "blast", "silva111");
	&cleanUp("uclust-ref", "blast", "gg2011"); &cleanUp("uclust-ref", "blast", "gg2013"); &cleanUp("uclust-ref", "blast", "silva111");
        &make_otu_table("uclust-ref", "rdp-training", "gg2011"); &make_otu_table("uclust-ref", "rdp-training", "gg2013"); &make_otu_table("uclust-ref", "rdp-training", "silva111");
        &filter_otus("uclust-ref", "rdp-training", "gg2011"); &filter_otus("uclust-ref", "rdp-training", "gg2013"); &filter_otus("uclust-ref", "rdp-training", "silva111");
        &removeRoot("uclust-ref", "rdp-training", "gg2011"); &removeRoot("uclust-ref", "rdp-training", "gg2013"); &removeRoot("uclust-ref", "rdp-training", "silva111");
        &summarize_taxa("uclust-ref", "rdp-training", "gg2011"); &summarize_taxa("uclust-ref", "rdp-training", "gg2013"); &summarize_taxa("uclust-ref", "rdp-training", "silva111");
        &alpha_div("uclust-ref", "rdp-training", "gg2011"); &alpha_div("uclust-ref", "rdp-training", "gg2013"); &alpha_div("uclust-ref", "rdp-training", "silva111");
        &beta_div("uclust-ref", "rdp-training", "gg2011"); &beta_div("uclust-ref", "rdp-training", "gg2013"); &beta_div("uclust-ref", "rdp-training", "silva111"); 
	&RAnalysis("uclust-ref", "rdp-training", "gg2011"); &RAnalysis("uclust-ref", "rdp-training", "gg2013"); &RAnalysis("uclust-ref", "rdp-training", "silva111");
	&cleanUp("uclust-ref", "rdp-training", "gg2011"); &cleanUp("uclust-ref", "rdp-training", "gg2013"); &cleanUp("uclust-ref", "rdp-training", "silva111");

	&make_otu_table("uclust-ref-strict", "blast", "gg2011"); &make_otu_table("uclust-ref-strict", "blast", "gg2013"); &make_otu_table("uclust-ref-strict", "blast", "silva111");
        &filter_otus("uclust-ref-strict", "blast", "gg2011"); &filter_otus("uclust-ref-strict", "blast", "gg2013"); &filter_otus("uclust-ref-strict", "blast", "silva111");
        &removeRoot("uclust-ref-strict", "blast", "gg2011"); &removeRoot("uclust-ref-strict", "blast", "gg2013"); &removeRoot("uclust-ref-strict", "blast", "silva111");
        &summarize_taxa("uclust-ref-strict", "blast", "gg2011"); &summarize_taxa("uclust-ref-strict", "blast", "gg2013"); &summarize_taxa("uclust-ref-strict", "blast", "silva111");
        &alpha_div("uclust-ref-strict", "blast", "gg2011"); &alpha_div("uclust-ref-strict", "blast", "gg2013"); &alpha_div("uclust-ref-strict", "blast", "silva111");
        &beta_div("uclust-ref-strict", "blast", "gg2011"); &beta_div("uclust-ref-strict", "blast", "gg2013"); &beta_div("uclust-ref-strict", "blast", "silva111");
	&RAnalysis("uclust-ref-strict", "blast", "gg2011"); &RAnalysis("uclust-ref-strict", "blast", "gg2013"); &RAnalysis("uclust-ref-strict", "blast", "silva111");
	&cleanUp("uclust-ref-stricit", "blast", "gg2011"); &cleanUp("uclust-ref-strict", "blast", "gg2013"); &cleanUp("uclust-ref-strict", "blast", "silva111");
        &make_otu_table("uclust-ref-strict", "rdp-training", "gg2011"); &make_otu_table("uclust-ref-strict", "rdp-training", "gg2013"); &make_otu_table("uclust-ref-strict", "rdp-training", "silva111");
        &filter_otus("uclust-ref-strict", "rdp-training", "gg2011"); &filter_otus("uclust-ref-strict", "rdp-training", "gg2013"); &filter_otus("uclust-ref-strict", "rdp-training", "silva111");
        &removeRoot("uclust-ref-strict", "rdp-training", "gg2011"); &removeRoot("uclust-ref-strict", "rdp-training", "gg2013"); &removeRoot("uclust-ref-strict", "rdp-training", "silva111");
        &summarize_taxa("uclust-ref-strict", "rdp-training", "gg2011"); &summarize_taxa("uclust-ref-strict", "rdp-training", "gg2013"); &summarize_taxa("uclust-ref-strict", "rdp-training", "silva111");
        &alpha_div("uclust-ref-strict", "rdp-training", "gg2011"); &alpha_div("uclust-ref-strict", "rdp-training", "gg2013"); &alpha_div("uclust-ref-strict", "rdp-training", "silva111");
        &beta_div("uclust-ref-strict", "rdp-training", "gg2011"); &beta_div("uclust-ref-strict", "rdp-training", "gg2013"); &beta_div("uclust-ref-strict", "rdp-training", "silva111");
	&RAnalysis("uclust-ref-strict", "rdp-training", "gg2011"); &RAnalysis("uclust-ref-strict", "rdp-training", "gg2013"); &RAnalysis("uclust-ref-strict", "rdp-training", "silva111");
	&cleanUp("uclust-ref-strict", "rdp-training", "gg2011"); &cleanUp("uclust-ref-strict", "rdp-training", "gg2013"); &cleanUp("uclust-ref-strict", "rdp-training", "silva111");

        &make_otu_table("blast", "blast", "gg2011"); &make_otu_table("blast", "blast", "gg2013"); &make_otu_table("blast", "blast", "silva111");
        &filter_otus("blast", "blast", "gg2011"); &filter_otus("blast", "blast", "gg2013"); &filter_otus("blast", "blast", "silva111");
        &removeRoot("blast", "blast", "gg2011"); &removeRoot("blast", "blast", "gg2013"); &removeRoot("blast", "blast", "silva111");
        &summarize_taxa("blast", "blast", "gg2011"); &summarize_taxa("blast", "blast", "gg2013"); &summarize_taxa("blast", "blast", "silva111");
        &alpha_div("blast", "blast", "gg2011"); &alpha_div("blast", "blast", "gg2013"); &alpha_div("blast", "blast", "silva111");
        &beta_div("blast", "blast", "gg2011"); &beta_div("blast", "blast", "gg2013"); &beta_div("blast", "blast", "silva111");
	&RAnalysis("blast", "blast", "gg2011"); &RAnalysis("blast", "blast", "gg2013"); &RAnalysis("blast", "blast", "silva111");
	&cleanUp("blast", "blast", "gg2011"); &cleanUp("blast", "blast", "gg2013"); &cleanUp("blast", "blast", "silva111");

        &make_otu_table("blast", "rdp-training", "gg2011"); &make_otu_table("blast", "rdp-training", "gg2013"); &make_otu_table("blast", "rdp-training", "silva111");
        &filter_otus("blast", "rdp-training", "gg2011"); &filter_otus("blast", "rdp-training", "gg2013"); &filter_otus("blast", "rdp-training", "silva111");
        &removeRoot("blast", "rdp-training", "gg2011"); &removeRoot("blast", "rdp-training", "gg2013"); &removeRoot("blast", "rdp-training", "silva111");
        &summarize_taxa("blast", "rdp-training", "gg2011"); &summarize_taxa("blast", "rdp-training", "gg2013"); &summarize_taxa("blast", "rdp-training", "silva111");
        &alpha_div("blast", "rdp-training", "gg2011"); &alpha_div("blast", "rdp-training", "gg2013"); &alpha_div("blast", "rdp-training", "silva111");
        &beta_div("blast", "rdp-training", "gg2011"); &beta_div("blast", "rdp-training", "gg2013"); &beta_div("blast", "rdp-training", "silva111");
	&RAnalysis("blast", "rdp-training", "gg2011"); &RAnalysis("blast", "rdp-training", "gg2013"); &RAnalysis("blast", "rdp-training", "silva111");
	&cleanUp("blast", "rdp-training", "gg2011"); &cleanUp("blast", "rdp-training", "gg2013"); &cleanUp("blast", "rdp-training", "silva111");

        &make_otu_table("uparse", "blast", "gg2011"); &make_otu_table("uparse", "blast", "gg2013"); &make_otu_table("uparse", "blast", "silva111");
        &filter_otus("uparse", "blast", "gg2011"); &filter_otus("uparse", "blast", "gg2013"); &filter_otus("uparse", "blast", "silva111");
        &removeRoot("uparse", "blast", "gg2011"); &removeRoot("uparse", "blast", "gg2013"); &removeRoot("uparse", "blast", "silva111");
        &summarize_taxa("uparse", "blast", "gg2011"); &summarize_taxa("uparse", "blast", "gg2013"); &summarize_taxa("uparse", "blast", "silva111");
        &alpha_div("uparse", "blast", "gg2011"); &alpha_div("uparse", "blast", "gg2013"); &alpha_div("uparse", "blast", "silva111");
        &beta_div("uparse", "blast", "gg2011"); &beta_div("uparse", "blast", "gg2013"); &beta_div("uparse", "blast", "silva111");
	&RAnalysis("uparse", "blast", "gg2011"); &RAnalysis("uparse", "blast", "gg2013"); &RAnalysis("uparse", "blast", "silva111");
	&cleanUp("uparse", "blast", "gg2011"); &cleanUp("uparse", "blast", "gg2013"); &cleanUp("uparse", "blast", "silva111");

        &make_otu_table("uparse", "rdp-training", "gg2011"); &make_otu_table("uparse", "rdp-training", "gg2013"); &make_otu_table("uparse", "rdp-training", "silva111");
        &filter_otus("uparse", "rdp-training", "gg2011"); &filter_otus("uparse", "rdp-training", "gg2013"); &filter_otus("uparse", "rdp-training", "silva111");
        &removeRoot("uparse", "rdp-training", "gg2011"); &removeRoot("uparse", "rdp-training", "gg2013"); &removeRoot("uparse", "rdp-training", "silva111");
        &summarize_taxa("uparse", "rdp-training", "gg2011"); &summarize_taxa("uparse", "rdp-training", "gg2013"); &summarize_taxa("uparse", "rdp-training", "silva111");
        &alpha_div("uparse", "rdp-training", "gg2011"); &alpha_div("uparse", "rdp-training", "gg2013"); &alpha_div("uparse", "rdp-training", "silva111");
        &beta_div("uparse", "rdp-training", "gg2011"); &beta_div("uparse", "rdp-training", "gg2013"); &beta_div("uparse", "rdp-training", "silva111");
	&RAnalysis("uparse", "rdp-training", "gg2011"); &RAnalysis("uparse", "rdp-training", "gg2013"); &RAnalysis("uparse", "rdp-training", "silva111");
	&cleanUp("uparse", "rdp-training", "gg2011"); &cleanUp("uparse", "rdp-training", "gg2013"); &cleanUp("uparse", "rdp-training", "silva111");
} elsif (($clust eq "all") && ($taxon eq "all")) {
	&make_otu_table("abundantotu", "blast", $gg);
        &filter_otus("abundantotu", "blast", $gg);
        &removeRoot("abundantotu", "blast", $gg);
        &summarize_taxa("abundantotu", "blast", $gg);
        &alpha_div("abundantotu", "blast", $gg);
        &beta_div("abundantotu", "blast", $gg);
	&RAnalysis("abundantotu", "blast", $gg);
	&cleanUp("abundantotu", "blast", $gg);
	&make_otu_table("abundantotu", "rdp-training", $gg);
        &filter_otus("abundantotu", "rdp-training", $gg);
        &removeRoot("abundantotu", "rdp-training", $gg);
        &summarize_taxa("abundantotu", "rdp-training", $gg);
        &alpha_div("abundantotu", "rdp-training", $gg);
        &beta_div("abundantotu", "rdp-training", $gg);
	&RAnalysis("abundantotu", "rdp-training", $gg);
	&cleanUp("abundantotu", "rdp-training", $gg);
	
	&make_otu_table("uclust", "blast", $gg);
        &filter_otus("uclust", "blast", $gg);
        &removeRoot("uclust", "blast", $gg);
        &summarize_taxa("uclust", "blast", $gg);
        &alpha_div("uclust", "blast", $gg);
        &beta_div("uclust", "blast", $gg);
	&RAnalysis("uclust", "blast", $gg);
	&cleanUp("uclust", "blast", $gg);
        &make_otu_table("uclust", "rdp-training", $gg);
        &filter_otus("uclust", "rdp-training", $gg);
        &removeRoot("uclust", "rdp-training", $gg);
        &summarize_taxa("uclust", "rdp-training", $gg);
        &alpha_div("uclust", "rdp-training", $gg);
        &beta_div("uclust", "rdp-training", $gg);
	&RAnalysis("uclust", "rdp-training", $gg);
	&cleanUp("uclust", "rdp-training", $gg);
	
	&make_otu_table("cdhit", "blast", $gg);
        &filter_otus("cdhit", "blast", $gg);
        &removeRoot("cdhit", "blast", $gg);
        &summarize_taxa("cdhit", "blast", $gg);
        &alpha_div("cdhit", "blast", $gg);
        &beta_div("cdhit", "blast", $gg);
	&RAnalysis("cdhit", "blast", $gg);
	&cleanUp("cdhit", "blast", $gg);
        &make_otu_table("cdhit", "rdp-training", $gg);
        &filter_otus("cdhit", "rdp-training", $gg);
        &removeRoot("cdhit", "rdp-training", $gg);
        &summarize_taxa("cdhit", "rdp-training", $gg);
        &alpha_div("cdhit", "rdp-training", $gg);
        &beta_div("cdhit", "rdp-training", $gg);
	&RAnalysis("cdhit", "rdp-training", $gg);
	&cleanUp("cdhit", "rdp-training", $gg);
	
	&make_otu_table("dnaclust", "blast", $gg);
        &filter_otus("dnaclust", "blast", $gg);
        &removeRoot("dnaclust", "blast", $gg);
        &summarize_taxa("dnaclust", "blast", $gg);
        &alpha_div("dnaclust", "blast", $gg);
        &beta_div("dnaclust", "blast", $gg);
	&RAnalysis("dnaclust", "blast", $gg);
	&cleanUp("dnaclust", "blast", $gg);
        &make_otu_table("dnaclust", "rdp-training", $gg);
        &filter_otus("dnaclust", "rdp-training", $gg);
        &removeRoot("dnaclust", "rdp-training", $gg);
        &summarize_taxa("dnaclust", "rdp-training", $gg);
        &alpha_div("dnaclust", "rdp-training", $gg);
        &beta_div("dnaclust", "rdp-training", $gg);
	&RAnalysis("dnaclust", "rdp-training", $gg);
	&cleanUp("dnaclust", "rdp-training", $gg);
	
	#&make_otu_table("mothur_average", "blast", $gg);
        #&filter_otus("mothur_average", "blast", $gg);
        #&removeRoot("mothur_average", "blast", $gg);
        #&summarize_taxa("mothur_average", "blast", $gg);
        #&alpha_div("mothur_average", "blast", $gg);
        #&beta_div("mothur_average", "blast", $gg);
	#&RAnalysis("mothur_average", "blast", $gg);
	#&cleanUp("mothur_average", "blast", $gg);
        #&make_otu_table("mothur_average", "rdp-training", $gg);
        #&filter_otus("mothur_average", "rdp-training", $gg);
        #&removeRoot("mothur_average", "rdp-training", $gg);
        #&summarize_taxa("mothur_average", "rdp-training", $gg);
        #&alpha_div("mothur_average", "rdp-training", $gg);
        #&beta_div("mothur_average", "rdp-training", $gg);
	#&RAnalysis("mothur_average", "rdp-training", $gg);
	#&cleanUp("mothur_average", "rdp-training", $gg);

	#&make_otu_table("mothur_nearest", "blast", $gg);
        #&filter_otus("mothur_nearest", "blast", $gg);
        #&removeRoot("mothur_nearest", "blast", $gg);
        #&summarize_taxa("mothur_nearest", "blast", $gg);
        #&alpha_div("mothur_nearest", "blast", $gg);
        #&beta_div("mothur_nearest", "blast", $gg);
	#&RAnalysis("mothur_nearest", "blast", $gg);
	#&cleanUp("mothur_nearest", "blast", $gg);
        #&make_otu_table("mothur_nearest", "rdp-training", $gg);
        #&filter_otus("mothur_nearest", "rdp-training", $gg);
        #&removeRoot("mothur_nearest", "rdp-training", $gg);
        #&summarize_taxa("mothur_nearest", "rdp-training", $gg);
        #&alpha_div("mothur_nearest", "rdp-training", $gg);
        #&beta_div("mothur_nearest", "rdp-training", $gg);
	#&RAnalysis("mothur_nearest", "rdp-training", $gg);
	#&cleanUp("mothur_nearest", "rdp-training", $gg);

	#&make_otu_table("mothur_furthest", "blast", $gg);
        #&filter_otus("mothur_furthest", "blast", $gg);
        #&removeRoot("mothur_furthest", "blast", $gg);
        #&summarize_taxa("mothur_furthest", "blast", $gg);
        #&alpha_div("mothur_furthest", "blast", $gg);
        #&beta_div("mothur_furthest", "blast", $gg);
	#&RAnalysis("mothur_furthest", "blast", $gg);
	#&cleanUp("mothur_furthest", "blast", $gg);
        #&make_otu_table("mothur_furthest", "rdp-training", $gg);
        #&filter_otus("mothur_furthest", "rdp-training", $gg);
        #&removeRoot("mothur_furthest", "rdp-training", $gg);
        #&summarize_taxa("mothur_furthest", "rdp-training", $gg);
        #&alpha_div("mothur_furthest", "rdp-training", $gg);
        #&beta_div("mothur_furthest", "rdp-training", $gg);
	#&RAnalysis("mothur_furthest", "rdp-training", $gg);
	#&cleanUp("mothur_furthest", "rdp-training", $gg);

	&make_otu_table("uclust-ref", "blast", $gg);
        &filter_otus("uclust-ref", "blast", $gg);
        &removeRoot("uclust-ref", "blast", $gg);
        &summarize_taxa("uclust-ref", "blast", $gg);
        &alpha_div("uclust-ref", "blast", $gg);
        &beta_div("uclust-ref", "blast", $gg);
	&RAnalysis("uclust-ref", "blast", $gg);
	&cleanUp("uclust-ref", "blast", $gg);
        &make_otu_table("uclust-ref", "rdp-training", $gg);
        &filter_otus("uclust-ref", "rdp-training", $gg);
        &removeRoot("uclust-ref", "rdp-training", $gg);
        &summarize_taxa("uclust-ref", "rdp-training", $gg);
        &alpha_div("uclust-ref", "rdp-training", $gg);
        &beta_div("uclust-ref", "rdp-training", $gg);
	&RAnalysis("uclust-ref", "rdp-training", $gg);
	&cleanUp("uclust-ref", "rdp-training", $gg);

	&make_otu_table("uclust-ref-strict", "blast", $gg);
        &filter_otus("uclust-ref-strict", "blast", $gg);
        &removeRoot("uclust-ref-strict", "blast", $gg);
        &summarize_taxa("uclust-ref-strict", "blast", $gg);
        &alpha_div("uclust-ref-strict", "blast", $gg);
        &beta_div("uclust-ref-strict", "blast", $gg);
	&RAnalysis("uclust-ref-strict", "blast", $gg);
	&cleanUp("uclust-ref-strict", "blast", $gg);
        &make_otu_table("uclust-ref-strict", "rdp-training", $gg);
        &filter_otus("uclust-ref-strict", "rdp-training", $gg);
        &removeRoot("uclust-ref-strict", "rdp-training", $gg);
        &summarize_taxa("uclust-ref-strict", "rdp-training", $gg);
        &alpha_div("uclust-ref-strict", "rdp-training", $gg);
        &beta_div("uclust-ref-strict", "rdp-training", $gg);
	&RAnalysis("uclust-ref-strict", "rdp-training", $gg);
	&cleanUp("uclust-ref-strict", "rdp-training", $gg);
		
	&make_otu_table("blast", "blast", $gg);
        &filter_otus("blast", "blast", $gg);
        &removeRoot("blast", "blast", $gg);
        &summarize_taxa("blast", "blast", $gg);
        &alpha_div("blast", "blast", $gg);
        &beta_div("blast", "blast", $gg);
	&RAnalysis("blast", "blast", $gg);
	&cleanUp("blast", "blast", $gg);
        &make_otu_table("blast", "rdp-training", $gg);
        &filter_otus("blast", "rdp-training", $gg);
        &removeRoot("blast", "rdp-training", $gg);
        &summarize_taxa("blast", "rdp-training", $gg);
        &alpha_div("blast", "rdp-training", $gg);
        &beta_div("blast", "rdp-training", $gg);
	&RAnalysis("blast", "rdp-training", $gg);
	&cleanUp("blast", "rdp-training", $gg);
	
	&make_otu_table("uparse", "blast", $gg);
        &filter_otus("uparse", "blast", $gg);
        &removeRoot("uparse", "blast", $gg);
        &summarize_taxa("uparse", "blast", $gg);
        &alpha_div("uparse", "blast", $gg);
        &beta_div("uparse", "blast", $gg);
	&RAnalysis("uparse", "blast", $gg);
	&cleanUp("uparse", "blast", $gg);
	
	&make_otu_table("uparse", "rdp-training", $gg);
        &filter_otus("uparse", "rdp-training", $gg);
        &removeRoot("uparse", "rdp-training", $gg);
        &summarize_taxa("uparse", "rdp-training", $gg);
        &alpha_div("uparse", "rdp-training", $gg);
        &beta_div("uparse", "rdp-training", $gg);
	&RAnalysis("uparse", "rdp-training", $gg);
	&cleanUp("uparse", "rdp-training", $gg);
} elsif (($clust eq "all") && ($gg eq "all")) {
    	&make_otu_table("abundantotu", $taxon, "gg2011");
    	&filter_otus("abundantotu", $taxon, "gg2011");
    	&removeRoot("abundantotu", $taxon, "gg2011");
    	&summarize_taxa("abundantotu", $taxon, "gg2011");
    	&alpha_div("abundantotu", $taxon, "gg2011");
    	&beta_div("abundantotu", $taxon, "gg2011");
	&RAnalysis("abundantotu", $taxon, "gg2011");
	&cleanUp("abundantotu", $taxon, "gg2011");
    	&make_otu_table("abundantotu", $taxon, "gg2013");
    	&filter_otus("abundantotu", $taxon, "gg2013");
    	&removeRoot("abundantotu", $taxon, "gg2013");
    	&summarize_taxa("abundantotu", $taxon, "gg2013");
    	&alpha_div("abundantotu", $taxon, "gg2013");
    	&beta_div("abundantotu", $taxon, "gg2013");
	&RAnalysis("abundantotu", $taxon, "gg2013");
	&cleanUp("abundantotu", $taxon, "gg2013");
    	&make_otu_table("abundantotu", $taxon, "silva111");
    	&filter_otus("abundantotu", $taxon, "silva111");
    	&removeRoot("abundantotu", $taxon, "silva111");
    	&summarize_taxa("abundantotu", $taxon, "silva111");
    	&alpha_div("abundantotu", $taxon, "silva111");
    	&beta_div("abundantotu", $taxon, "silva111");
	&RAnalysis("abundantotu", $taxon, "silva111");
	&cleanUp("abundantotu", $taxon, "silva111");
    	
    	&make_otu_table("uclust", $taxon, "gg2011");
    	&filter_otus("uclust", $taxon, "gg2011");
    	&removeRoot("uclust", $taxon, "gg2011");
    	&summarize_taxa("uclust", $taxon, "gg2011");
    	&alpha_div("uclust", $taxon, "gg2011");
    	&beta_div("uclust", $taxon, "gg2011");
	&RAnalysis("uclust", $taxon, "gg2011");
	&cleanUp("uclust", $taxon, "gg2011");
    	&make_otu_table("uclust", $taxon, "gg2013");
    	&filter_otus("uclust", $taxon, "gg2013");
    	&removeRoot("uclust", $taxon, "gg2013");
    	&summarize_taxa("uclust", $taxon, "gg2013");
    	&alpha_div("uclust", $taxon, "gg2013");
    	&beta_div("uclust", $taxon, "gg2013");
	&RAnalysis("uclust", $taxon, "gg2013");
	&cleanUp("uclust", $taxon, "gg2013");
    	&make_otu_table("uclust", $taxon, "silva111");
    	&filter_otus("uclust", $taxon, "silva111");
    	&removeRoot("uclust", $taxon, "silva111");
    	&summarize_taxa("uclust", $taxon, "silva111");
    	&alpha_div("uclust", $taxon, "silva111");
    	&beta_div("uclust", $taxon, "silva111");
	&RAnalysis("uclust", $taxon, "silva111");
	&cleanUp("uclust", $taxon, "silva111");
    	
    	&make_otu_table("cdhit", $taxon, "gg2011");
    	&filter_otus("cdhit", $taxon, "gg2011");
    	&removeRoot("cdhit", $taxon, "gg2011");
    	&summarize_taxa("cdhit", $taxon, "gg2011");
    	&alpha_div("cdhit", $taxon, "gg2011");
    	&beta_div("cdhit", $taxon, "gg2011");
	&RAnalysis("cdhit", $taxon, "gg2011");
	&cleanUp("cdhit", $taxon, "gg2011");
    	&make_otu_table("cdhit", $taxon, "gg2013");
    	&filter_otus("cdhit", $taxon, "gg2013");
    	&removeRoot("cdhit", $taxon, "gg2013");
    	&summarize_taxa("cdhit", $taxon, "gg2013");
    	&alpha_div("cdhit", $taxon, "gg2013");
    	&beta_div("cdhit", $taxon, "gg2013");
	&RAnalysis("cdhit", $taxon, "gg2013");
	&cleanUp("cdhit", $taxon, "gg2013");
    	&make_otu_table("cdhit", $taxon, "silva111");
    	&filter_otus("cdhit", $taxon, "silva111");
   	&removeRoot("cdhit", $taxon, "silva111");
   	&summarize_taxa("cdhit", $taxon, "silva111");
   	&alpha_div("cdhit", $taxon, "silva111");
    	&beta_div("cdhit", $taxon, "silva111");
	&RAnalysis("cdhit", $taxon, "silva111");
	&cleanUp("cdhit", $taxon, "silva111");
    	
    	&make_otu_table("dnaclust", $taxon, "gg2011");
    	&filter_otus("dnaclust", $taxon, "gg2011");
    	&removeRoot("dnaclust", $taxon, "gg2011");
    	&summarize_taxa("dnaclust", $taxon, "gg2011");
    	&alpha_div("dnaclust", $taxon, "gg2011");
    	&beta_div("dnaclust", $taxon, "gg2011");
	&RAnalysis("dnaclust", $taxon, "gg2011");
	&cleanUp("dnaclust", $taxon, "gg2011");
    	&make_otu_table("dnaclust", $taxon, "gg2013");
    	&filter_otus("dnaclust", $taxon, "gg2013");
    	&removeRoot("dnaclust", $taxon, "gg2013");
    	&summarize_taxa("dnaclust", $taxon, "gg2013");
    	&alpha_div("dnaclust", $taxon, "gg2013");
    	&beta_div("dnaclust", $taxon, "gg2013");
	&RAnalysis("dnaclust", $taxon, "gg2013");
	&cleanUp("dnaclust", $taxon, "gg2013");
    	&make_otu_table("dnaclust", $taxon, "silva111");
    	&filter_otus("dnaclust", $taxon, "silva111");
    	&removeRoot("dnaclust", $taxon, "silva111");
    	&summarize_taxa("dnaclust", $taxon, "silva111");
    	&alpha_div("dnaclust", $taxon, "silva111");
    	&beta_div("dnaclust", $taxon, "silva111");
	&RAnalysis("dnaclust", $taxon, "silva111");
	&cleanUp("dnaclust", $taxon, "silva111");
    	
    	#&make_otu_table("mothur_average", $taxon, "gg2011");
    	#&filter_otus("mothur_average", $taxon, "gg2011");
    	#&removeRoot("mothur_average", $taxon, "gg2011");
    	#&summarize_taxa("mothur_average", $taxon, "gg2011");
    	#&alpha_div("mothur_average", $taxon, "gg2011");
    	#&beta_div("mothur_average", $taxon, "gg2011");
	#&RAnalysis("mothur_average", $taxon, "gg2011");
	#&cleanUp("mothur_average", $taxon, "gg2011");
    	#&make_otu_table("mothur_average", $taxon, "gg2013");
    	#&filter_otus("mothur_average", $taxon, "gg2013");
    	#&removeRoot("mothur_average", $taxon, "gg2013");
    	#&summarize_taxa("mothur_average", $taxon, "gg2013");
    	#&alpha_div("mothur_average", $taxon, "gg2013");
    	#&beta_div("mothur_average", $taxon, "gg2013");
	#&RAnalysis("mothur_average", $taxon, "gg2013");
	#&cleanUp("mothur_average", $taxon, "gg2013");
    	#&make_otu_table("mothur_average", $taxon, "silva111");
    	#&filter_otus("mothur_average", $taxon, "silva111");
    	#&removeRoot("mothur_average", $taxon, "silva111");
    	#&summarize_taxa("mothur_average", $taxon, "silva111");
    	#&alpha_div("mothur_average", $taxon, "silva111");
    	#&beta_div("mothur_average", $taxon, "silva111");
	#&RAnalysis("mothur_average", $taxon, "silva111");
	#&cleanUp("mothur_average", $taxon, "silva111");
    	
    	#&make_otu_table("mothur_nearest", $taxon, "gg2011");
    	#&filter_otus("mothur_nearest", $taxon, "gg2011");
    	#&removeRoot("mothur_nearest", $taxon, "gg2011");
    	#&summarize_taxa("mothur_nearest", $taxon, "gg2011");
    	#&alpha_div("mothur_nearest", $taxon, "gg2011");
    	#&beta_div("mothur_nearest", $taxon, "gg2011");
	#&RAnalysis("mothur_nearest", $taxon, "gg2011");
	#&cleanUp("mothur_nearest", $taxon, "gg2011");
    	#&make_otu_table("mothur_nearest", $taxon, "gg2013");
    	#&filter_otus("mothur_nearest", $taxon, "gg2013");
    	#&removeRoot("mothur_nearest", $taxon, "gg2013");
    	#&summarize_taxa("mothur_nearest", $taxon, "gg2013");
    	#&alpha_div("mothur_nearest", $taxon, "gg2013");
    	#&beta_div("mothur_nearest", $taxon, "gg2013");
	#&RAnalysis("mothur_nearest", $taxon, "gg2013");
	#&cleanUp("mothur_nearest", $taxon, "gg2013");
    	#&make_otu_table("mothur_nearest", $taxon, "silva111");
    	#&filter_otus("mothur_nearest", $taxon, "silva111");
    	#&removeRoot("mothur_nearest", $taxon, "silva111");
    	#&summarize_taxa("mothur_nearest", $taxon, "silva111");
    	#&alpha_div("mothur_nearest", $taxon, "silva111");
    	#&beta_div("mothur_nearest", $taxon, "silva111");
	#&RAnalysis("mothur_nearest", $taxon, "silva111");
	#&cleanUp("mothur_nearest", $taxon, "silva111");
    	
    	#&make_otu_table("mothur_furthest", $taxon, "gg2011");
    	#&filter_otus("mothur_furthest", $taxon, "gg2011");
    	#&removeRoot("mothur_furthest", $taxon, "gg2011");
    	#&summarize_taxa("mothur_furthest", $taxon, "gg2011");
    	#&alpha_div("mothur_furthest", $taxon, "gg2011");
    	#&beta_div("mothur_furthest", $taxon, "gg2011");
	#&RAnalysis("mothur_furthest", $taxon, "gg2011");
	#&cleanUp("mothur_furthest", $taxon, "gg2011");
    	#&make_otu_table("mothur_furthest", $taxon, "gg2013");
    	#&filter_otus("mothur_furthest", $taxon, "gg2013");
    	#&removeRoot("mothur_furthest", $taxon, "gg2013");
    	#&summarize_taxa("mothur_furthest", $taxon, "gg2013");
    	#&alpha_div("mothur_furthest", $taxon, "gg2013");
    	#&beta_div("mothur_furthest", $taxon, "gg2013");
	#&RAnalysis("mothur_furthest", $taxon, "gg2013");
	#&cleanUp("mothur_furthest", $taxon, "gg2013");
    	#&make_otu_table("mothur_furthest", $taxon, "silva111");
    	#&filter_otus("mothur_furthest", $taxon, "silva111");
    	#&removeRoot("mothur_furthest", $taxon, "silva111");
    	#&summarize_taxa("mothur_furthest", $taxon, "silva111");
    	#&alpha_div("mothur_furthest", $taxon, "silva111");
    	#&beta_div("mothur_furthest", $taxon, "silva111");
	#&RAnalysis("mothur_furthest", $taxon, "silva111");
	#&cleanUp("mothur_furthest", $taxon, "silva111");
    	
    	&make_otu_table("uclust-ref", $taxon, "gg2011");
    	&filter_otus("uclust-ref", $taxon, "gg2011");
   	&removeRoot("uclust-ref", $taxon, "gg2011");
   	&summarize_taxa("uclust-ref", $taxon, "gg2011");
   	&alpha_div("uclust-ref", $taxon, "gg2011");
   	&beta_div("uclust-ref", $taxon, "gg2011");
	&RAnalysis("uclust-ref", $taxon, "gg2011");
	&cleanUp("uclust-ref", $taxon, "gg2011");
   	&make_otu_table("uclust-ref", $taxon, "gg2013");
    	&filter_otus("uclust-ref", $taxon, "gg2013");
    	&removeRoot("uclust-ref", $taxon, "gg2013");
    	&summarize_taxa("uclust-ref", $taxon, "gg2013");
    	&alpha_div("uclust-ref", $taxon, "gg2013");
    	&beta_div("uclust-ref", $taxon, "gg2013");
	&RAnalysis("uclust-ref", $taxon, "gg2013");
	&cleanUp("uclust-ref", $taxon, "gg2013");
    	&make_otu_table("uclust-ref", $taxon, "silva111");
    	&filter_otus("uclust-ref", $taxon, "silva111");
    	&removeRoot("uclust-ref", $taxon, "silva111");
    	&summarize_taxa("uclust-ref", $taxon, "silva111");
    	&alpha_div("uclust-ref", $taxon, "silva111");
    	&beta_div("uclust-ref", $taxon, "silva111");
	&RAnalysis("uclust-ref", $taxon, "silva111");
	&cleanUp("uclust-ref", $taxon, "silva111");
    	
    	&make_otu_table("uclust-ref-strict", $taxon, "gg2011");
    	&filter_otus("uclust-ref-strict", $taxon, "gg2011");
    	&removeRoot("uclust-ref-strict", $taxon, "gg2011");
    	&summarize_taxa("uclust-ref-strict", $taxon, "gg2011");
    	&alpha_div("uclust-ref-strict", $taxon, "gg2011");
    	&beta_div("uclust-ref-strict", $taxon, "gg2011");
	&RAnalysis("uclust-ref-strict", $taxon, "gg2011");
	&cleanUp("uclust-ref-strict", $taxon, "gg2011");
    	&make_otu_table("uclust-ref-strict", $taxon, "gg2013");
    	&filter_otus("uclust-ref-strict", $taxon, "gg2013");
    	&removeRoot("uclust-ref-strict", $taxon, "gg2013");
    	&summarize_taxa("uclust-ref-strict", $taxon, "gg2013");
    	&alpha_div("uclust-ref-strict", $taxon, "gg2013");
    	&beta_div("uclust-ref-strict", $taxon, "gg2013");
	&RAnalysis("uclust-ref-strict", $taxon, "gg2013");
	&cleanUp("uclust-ref-strict", $taxon, "gg2013");
    	&make_otu_table("uclust-ref-strict", $taxon, "silva111");
    	&filter_otus("uclust-ref-strict", $taxon, "silva111");
    	&removeRoot("uclust-ref-strict", $taxon, "silva111");
    	&summarize_taxa("uclust-ref-strict", $taxon, "silva111");
    	&alpha_div("uclust-ref-strict", $taxon, "silva111");
    	&beta_div("uclust-ref-strict", $taxon, "silva111");
	&RAnalysis("uclust-ref-strict", $taxon, "silva111");
	&cleanUp("uclust-ref-strict", $taxon, "silva111");
    
    	&make_otu_table("blast", $taxon, "gg2011");
    	&filter_otus("blast", $taxon, "gg2011");
    	&removeRoot("blast", $taxon, "gg2011");
    	&summarize_taxa("blast", $taxon, "gg2011");
    	&alpha_div("blast", $taxon, "gg2011");
    	&beta_div("blast", $taxon, "gg2011");
    	&RAnalysis("blast", $taxon, "gg2011");
	&cleanUp("blast", $taxon, "gg2011");
	&make_otu_table("blast", $taxon, "gg2013");
    	&filter_otus("blast", $taxon, "gg2013");
    	&removeRoot("blast", $taxon, "gg2013");
    	&summarize_taxa("blast", $taxon, "gg2013");
    	&alpha_div("blast", $taxon, "gg2013");
    	&beta_div("blast", $taxon, "gg2013");
	&RAnalysis("blast", $taxon, "gg2013");
	&cleanUp("blast", $taxon, "gg2013");
    	&make_otu_table("blast", $taxon, "silva111");
    	&filter_otus("blast", $taxon, "silva111");
    	&removeRoot("blast", $taxon, "silva111");
    	&summarize_taxa("blast", $taxon, "silva111");
    	&alpha_div("blast", $taxon, "silva111");
    	&beta_div("blast", $taxon, "silva111");
	&RAnalysis("blast", $taxon, "silva111");
	&cleanUp("blast", $taxon, "silva111");
    
    	&make_otu_table("uparse", $taxon, "gg2011");
    	&filter_otus("uparse", $taxon, "gg2011");
    	&removeRoot("uparse", $taxon, "gg2011");
    	&summarize_taxa("uparse", $taxon, "gg2011");
    	&alpha_div("uparse", $taxon, "gg2011");
    	&beta_div("uparse", $taxon, "gg2011");
	&RAnalysis("uparse", $taxon, "gg2011");
	&cleanUp("uparse", $taxon, "gg2011");
    	&make_otu_table("uparse", $taxon, "gg2013");
    	&filter_otus("uparse", $taxon, "gg2013");
    	&removeRoot("uparse", $taxon, "gg2013");
    	&summarize_taxa("uparse", $taxon, "gg2013");
    	&alpha_div("uparse", $taxon, "gg2013");
    	&beta_div("uparse", $taxon, "gg2013");
	&RAnalysis("uparse", $taxon, "gg2013");
	&cleanUp("uparse", $taxon, "gg2013");
    	&make_otu_table("uparse", $taxon, "silva111");
    	&filter_otus("uparse", $taxon, "silva111");
    	&removeRoot("uparse", $taxon, "silva111");
    	&summarize_taxa("uparse", $taxon, "silva111");
    	&alpha_div("uparse", $taxon, "silva111");
	&beta_div("uparse", $taxon, "silva111");
	&RAnalysis("uparse", $taxon, "silva111");
	&cleanUp("uparse", $taxon, "silva111");
} elsif (($taxon eq "all") && ($gg eq "all")) {
	&make_otu_table($clust, "blast", "gg2011");
    	&filter_otus($clust, "blast", "gg2011");
    	&removeRoot($clust, "blast", "gg2011");
    	&summarize_taxa($clust, "blast", "gg2011");
    	&alpha_div($clust, "blast", "gg2011");
    	&beta_div($clust, "blast", "gg2011");
	&RAnalysis($clust, "blast", "gg2011");
	&cleanUp($clust, "blast", "gg2011");
    	&make_otu_table($clust, "blast", "gg2013");
    	&filter_otus($clust, "blast", "gg2013");
    	&removeRoot($clust, "blast", "gg2013");
    	&summarize_taxa($clust, "blast", "gg2013");
    	&alpha_div($clust, "blast", "gg2013");
    	&beta_div($clust, "blast", "gg2013");
	&RAnalysis($clust, "blast", "gg2013");
	&cleanUp($clust, "blast", "gg2013");
    	&make_otu_table($clust, "blast", "silva111");
    	&filter_otus($clust, "blast", "silva111");
    	&removeRoot($clust, "blast", "silva111");
    	&summarize_taxa($clust, "blast", "silva111");
    	&alpha_div($clust, "blast", "silva111");
    	&beta_div($clust, "blast", "silva111");
	&RAnalysis($clust, "blast", "silva111");
	&cleanUp($clust, "blast", "silva111");

	&make_otu_table($clust, "rdp-training", "gg2011");
        &filter_otus($clust, "rdp-training", "gg2011");
        &removeRoot($clust, "rdp-training", "gg2011");
        &summarize_taxa($clust, "rdp-training", "gg2011");
        &alpha_div($clust, "rdp-training", "gg2011");
        &beta_div($clust, "rdp-training", "gg2011");
	&RAnalysis($clust, "rdp-training", "gg2011");
	&cleanUp($clust, "rdp-training", "gg2011");
        &make_otu_table($clust, "rdp-training", "gg2013");
        &filter_otus($clust, "rdp-training", "gg2013");
        &removeRoot($clust, "rdp-training", "gg2013");
        &summarize_taxa($clust, "rdp-training", "gg2013");
        &alpha_div($clust, "rdp-training", "gg2013");
        &beta_div($clust, "rdp-training", "gg2013");
	&RAnalysis($clust, "rdp-training", "gg2013");
	&cleanUp($clust, "rdp-training", "gg2013");
        &make_otu_table($clust, "rdp-training", "silva111");
        &filter_otus($clust, "rdp-training", "silva111");
        &removeRoot($clust, "rdp-training", "silva111");
        &summarize_taxa($clust, "rdp-training", "silva111");
        &alpha_div($clust, "rdp-training", "silva111");
        &beta_div($clust, "rdp-training", "silva111");
	&RAnalysis($clust, "rdp-training", "silva111");
	&cleanUp($clust, "rdp-training", "silva111");
} elsif ($clust eq "all") {
        &make_otu_table("abundantotu", $taxon, $gg);
        &filter_otus("abundantotu", $taxon, $gg);
        &removeRoot("abundantotu", $taxon, $gg);
        &summarize_taxa("abundantotu", $taxon, $gg);
        &alpha_div("abundantotu", $taxon, $gg);
        &beta_div("abundantotu", $taxon, $gg);
	&RAnalysis("abundantotu", $taxon, $gg);
	&cleanUp("abundantotu", $taxon, $gg);
        &make_otu_table("uclust", $taxon, $gg);
        &filter_otus("uclust", $taxon, $gg);
        &removeRoot("uclust", $taxon, $gg);
        &summarize_taxa("uclust", $taxon, $gg);
        &alpha_div("uclust", $taxon, $gg);
        &beta_div("uclust", $taxon, $gg);
	&RAnalysis("uclust", $taxon, $gg);
	&cleanUp("uclust", $taxon, $gg);
        &make_otu_table("cdhit", $taxon, $gg);
        &filter_otus("cdhit", $taxon, $gg);
        &removeRoot("cdhit", $taxon, $gg);
        &summarize_taxa("cdhit", $taxon, $gg);
        &alpha_div("cdhit", $taxon, $gg);
        &beta_div("cdhit", $taxon, $gg);
	&RAnalysis("cdhit", $taxon, $gg);
	&cleanUp("cdhit", $taxon, $gg);
        &make_otu_table("dnaclust", $taxon, $gg);
        &filter_otus("dnaclust", $taxon, $gg);
        &removeRoot("dnaclust", $taxon, $gg);
        &summarize_taxa("dnaclust", $taxon, $gg);
        &alpha_div("dnaclust", $taxon, $gg);
        &beta_div("dnaclust", $taxon, $gg);
	&RAnalysis("dnaclust", $taxon, $gg);
	&cleanUp("dnaclust", $taxon, $gg);
        #&make_otu_table("mothur_average", $taxon, $gg);
        #&filter_otus("mothur_average", $taxon, $gg);
        #&removeRoot("mothur_average", $taxon, $gg);
        #&summarize_taxa("mothur_average", $taxon, $gg);
        #&alpha_div("mothur_average", $taxon, $gg);
        #&beta_div("mothur_average", $taxon, $gg);
	#&RAnalysis("mothur_average", $taxon, $gg);
	#&cleanUp("mothur_average", $taxon, $gg);
        #&make_otu_table("mothur_nearest", $taxon, $gg);
        #&filter_otus("mothur_nearest", $taxon, $gg);
        #&removeRoot("mothur_nearest", $taxon, $gg);
        #&summarize_taxa("mothur_nearest", $taxon, $gg);
        #&alpha_div("mothur_nearest", $taxon, $gg);
        #&beta_div("mothur_nearest", $taxon, $gg);
	#&RAnalysis("mothur_nearest", $taxon, $gg);
	#&cleanUp("mothur_nearest", $taxon, $gg);
        #&make_otu_table("mothur_furthest", $taxon, $gg);
        #&filter_otus("mothur_furthest", $taxon, $gg);
        #&removeRoot("mothur_furthest", $taxon, $gg);
        #&summarize_taxa("mothur_furthest", $taxon, $gg);
        #&alpha_div("mothur_furthest", $taxon, $gg);
        #&beta_div("mothur_furthest", $taxon, $gg);
	#&RAnalysis("mothur_furthest", $taxon, $gg);
	#&cleanUp("mothur_furthest", $taxon, $gg);
        &make_otu_table("uclust-ref", $taxon, $gg);
        &filter_otus("uclust-ref", $taxon, $gg);
        &removeRoot("uclust-ref", $taxon, $gg);
        &summarize_taxa("uclust-ref", $taxon, $gg);
        &alpha_div("uclust-ref", $taxon, $gg);
        &beta_div("uclust-ref", $taxon, $gg);
	&RAnalysis("uclust-ref", $taxon, $gg);
	&cleanUp("uclust-ref", $taxon, $gg);
	&make_otu_table("uclust-ref-strict", $taxon, $gg);
	&filter_otus("uclust-ref-strict", $taxon, $gg);
	&removeRoot("uclust-ref-strict", $taxon, $gg);
	&summarize_taxa("uclust-ref-strict", $taxon, $gg);
	&alpha_div("uclust-ref-strict", $taxon, $gg);
	&beta_div("uclust-ref-strict", $taxon, $gg);
	&RAnalysis("uclust-ref-strict", $taxon, $gg);
	&cleanUp("uclust-ref-strict", $taxon, $gg);
        &make_otu_table("blast", $taxon, $gg);
        &filter_otus("blast", $taxon, $gg);
        &removeRoot("blast", $taxon, $gg);
        &summarize_taxa("blast", $taxon, $gg);
        &alpha_div("blast", $taxon, $gg);
        &beta_div("blast", $taxon, $gg);
	&RAnalysis("blast", $taxon, $gg);
	&cleanUp("blast", $taxon, $gg);
	&make_otu_table("uparse", $taxon, $gg);
        &filter_otus("uparse", $taxon, $gg);
        &removeRoot("uparse", $taxon, $gg);
        &summarize_taxa("uparse", $taxon, $gg);
        &alpha_div("uparse", $taxon, $gg);
        &beta_div("uparse", $taxon, $gg);
	&RAnalysis("uparse", $taxon, $gg);
	&cleanUp("uparse", $taxon, $gg);
} elsif ($taxon eq "all") {
	&make_otu_table($clust, "blast", $gg);
        &filter_otus($clust, "blast", $gg);
        &removeRoot($clust, "blast", $gg);
        &summarize_taxa($clust, "blast", $gg);
        &alpha_div($clust, "blast", $gg);
        &beta_div($clust, "blast", $gg);
	&RAnalysis($clust, "blast", $gg);
	&cleanUp($clust, "blast", $gg);
	&make_otu_table($clust, "rdp-training", $gg);
        &filter_otus($clust, "rdp-training", $gg);
        &removeRoot($clust, "rdp-training", $gg);
        &summarize_taxa($clust, "rdp-training", $gg);
        &alpha_div($clust, "rdp-training", $gg);
        &beta_div($clust, "rdp-training", $gg);
	&RAnalysis($clust, "rdp-training", $gg);
	&cleanUp($clust, "rdp-training", $gg);
} elsif ($gg eq "all") {
	&make_otu_table($clust, $taxon, "gg2011"); &make_otu_table($clust, $taxon, "gg2013"); &make_otu_table($clust, $taxon, "silva111");
	&filter_otus($clust, $taxon, "gg2011"); &filter_otus($clust, $taxon, "gg2013"); &filter_otus($clust, $taxon, "silva111");
	&removeRoot($clust, $taxon, "gg2011"); &removeRoot($clust, $taxon, "gg2013"); &removeRoot($clust, $taxon, "silva111");
	&summarize_taxa($clust, $taxon, "gg2011"); &summarize_taxa($clust, $taxon, "gg2013"); &summarize_taxa($clust, $taxon, "silva111");
	&alpha_div($clust, $taxon, "gg2011"); &alpha_div($clust, $taxon, "gg2013"); &alpha_div($clust, $taxon, "silva111");
	&beta_div($clust, $taxon, "gg2011"); &beta_div($clust, $taxon, "gg2013"); &beta_div($clust, $taxon, "silva111");
	&RAnalysis($clust, $taxon, "gg2011"); &RAnalysis($clust, $taxon, "gg2013"); &RAnalysis($clust, $taxon, "silva111");
	&cleanUp($clust, $taxon, "gg2011"); &cleanUp($clust, $taxon, "gg2013"); &cleanUp($clust, $taxon, "silva111");
} else {
	&make_otu_table($clust, $taxon, $gg);
	&filter_otus($clust, $taxon, $gg);
	&removeRoot($clust, $taxon, $gg);
	&summarize_taxa($clust, $taxon, $gg);
	&alpha_div($clust, $taxon, $gg);
	&beta_div($clust, $taxon, $gg);
	&RAnalysis($clust, $taxon, $gg);
	&cleanUp($clust, $taxon, $gg);
}
#Clean up the files in the main folder.
$cmd = "tar -zcvf splits_dir.tar.gz splits_dir/ && rm -R splits_dir/";
system($cmd);
$cmd = "tar -zcvf $proj.fna.tar.gz $proj.fna && rm $proj.fna";
system($cmd);
close LOG;
close ERR;
exit;

###############################################################
#Subroutines
###############################################################

####SUBROUTINE: getprimers####
sub getPrimers($$) {
        my ($barloc, $region) = @_;
	$fwd_primer = "";
        $rev_primer = "";
        $fwd_revcomp_primer = "";
        $rev_revcomp_primer = "";
        $overlap = 1;
        if ($barloc eq "rev") {
                        if ($region eq "v3") {
                                $fwd_primer = "CCTACGGGAGGCAGCAG";
                                $rev_primer = "ATTACCGCGGCTGCTGG";
                                $fwd_revcomp_primer = "CTGCTGCCTCCCGTAGG";
                                $rev_revcomp_primer = "CCAGCAGCCGCGGTAAT";
                                $overlap = 1;
				#Length of amplicon is expected to be 161bp (341F to 518R with a 16bp forward primer)
				#Give 50bp leeway on either side.
				$splitlib_min = 111;
				$splitlib_max = 211;
                        } elsif ($region eq "v34") {
                                $fwd_primer = "CCTACGGGAGGCAGCAG";
                                $rev_primer = "GGACTACHVGGGTWTCTAAT";
                                $fwd_revcomp_primer = "CTGCTGCCTCCCGTAGG";
                                $rev_revcomp_primer = "ATTAGAWACCCBDGTAGTCC";
                                $overlap = 1;
				#Length of amplicon is expected to be 481bp (341F to 806R with a 16bp forward primer)
				$splitlib_min = 380;
				$splitlib_max = 480;
                        } elsif ($region eq "v4") {
				$fwd_primer = "GTGYCAGCMGCCGCGGTAA";
				$rev_primer = "GGACTACNVGGGTWTCTAAT";
				$fwd_revcomp_primer = "TTACCGCGGCKGCTCRCAC";
				$rev_revcomp_primer = "ATTAGAWACCCBNGTAGTCC";
				$overlap = 1;
				#Length of amplicon is expected to be 275bp (515F to 806R with a 16bp forward primer)
				$splitlib_min = 225;
				$splitlib_max = 325;
			}
	} elsif ($barloc eq "fwd") {
                        if ($region eq "v3") {
                                $fwd_primer = "ATTACCGCGGCTGCTGG";
                                $rev_primer = "CCTACGGGAGGCAGCAG";
                                $fwd_revcomp_primer = "CCAGCAGCCGCGGTAAT";
                                $rev_revcomp_primer = "CTGCTGCCTCCCGTAGG";
                                $overlap = 1;
				$splitlib_min = 111;
				$splitlib_max = 211;
                        } elsif ($region eq "v34") {
                                $fwd_primer = "GGACTACHVGGGTWTCTAAT";
                                $rev_primer = "CCTACGGGAGGCAGCAG";
                                $fwd_revcomp_primer = "ATTAGAWACCCBDGTAGTCC";
                                $rev_revcomp_primer = "CTGCTGCCTCCCGTAGG";
                                $overlap = 1;
				$splitlib_min = 380;
				$splitlib_max = 480;
                        } elsif ($region eq "v4") {
				$fwd_primer = "GGACTACNVGGGTWTCTAAT";
				$rev_primer = "GTGYCAGCMGCCGCGGTAA";
				$fwd_revcomp_primer = "ATTAGAWACCCBNGTAGTCC";
				$rev_revcomp_primer = "TTACCGCGGCKGCTCRCAC";
				$overlap = 1;
				$splitlib_min = 225;
				$splitlib_max = 325;
			}
        }
        return ($fwd_primer, $rev_primer, $fwd_revcomp_primer, $rev_revcomp_primer, $overlap, $splitlib_min, $splitlib_max);
}

#### SUBROUTINE: getPrimersFile ####
sub getPrimersFile($) {
	my ($seqinfofile) = @_;
	open(SEQIN, "<", $seqinfofile) or die "Cannot open -s sequence information file";
	my @in = <SEQIN>;
	chomp @in;
	close SEQIN;
	$fwd_primer = $in[0];
	chomp($fwd_primer);
	$fwd_revcomp_primer = reverse($fwd_primer);
	$fwd_revcomp_primer =~ tr/ACGTacgt/TGCAtgca/;
	chomp($fwd_revcomp_primer);
	$rev_primer = $in[1];
	chomp($rev_primer);
	$rev_revcomp_primer = reverse($rev_primer);
	$rev_revcomp_primer =~ tr/ACGTacgt/TGCAtgca/;
	chomp($rev_revcomp_primer);
	$overlap = 1;
	$splitlib_min = $in[2];
	chomp($splitlib_min);
	$splitlib_max = $in[3];
	chomp($splitlib_max);
	$barloc = $in[4];
	chomp($barloc);
	return($fwd_primer, $rev_primer, $fwd_revcomp_primer, $rev_revcomp_primer, $overlap, $splitlib_min, $splitlib_max, $barloc);
}

#### SUBROUTINE: AbundantOTU+ ####
sub abundantOTU {
	my ($gg) = @_;
	print "pick otus: AbundantOTU+\n";
	chdir($time_pwd);
	unless (-e "picked_otus_abundantotu_$gg") {
        	system("mkdir picked_otus_abundantotu_$gg");
       	}    
        if ($? == -1) {print "command failed: $!\n"; exit; }
        chdir("picked_otus_abundantotu_$gg");
        if ($? == -1) {print "command failed: $!\n"; exit; }
        my $dis = (100 - $clustThres)/100;
        $cmd = "$bin/AbundantOTU+0.93b/bin/AbundantOTU+ -d $dis -i ../splits_dir/seqs_$gg.fna -o rep_set 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time $bin/AbundantOTU+0.93b/bin/AbundantOTU+ -d $dis -i ../splits_dir/seqs_$gg.fna -o rep_set) > $time_pwd/time_pick_otus_abundantotu_$taxon.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) {print "command failed: $!\n"; exit; }
        if ($? == -1) {print "command failed: $!\n"; exit; }
        $cmd = "$bin/abundantOTU_to_qiime.pl rep_set.clust seqs_".$gg."_otus.txt 2>&1 | tee -a $err";
	system($cmd);
        print LOG $cmd."\n";
	if ($? == -1) {print "command failed: $!\n"; exit; }
	$cmd = "mv rep_set.cons rep_set.fna";
        system($cmd);
	print LOG $cmd."\n";
        if ($? == -1) {print "command failed: $!\n"; exit; }
	$cmd = "perl -i -pe 's/Consensus//g' rep_set.fna";
        system($cmd);
	print LOG $cmd."\n";
        if ($? == -1) {print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: UCLUST ####
sub uclust {
	my ($gg) = @_;
	print "pick otus: UCLUST\n";
        chdir($time_pwd);
	my $dis = $clustThres/100;
	$cmd="pick_otus.py -m uclust -s $dis -o picked_otus_uclust_$gg -i splits_dir/seqs_$gg.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd="(time pick_otus.py -m uclust -s $dis -o picked_otus_uclust_$gg -i splits_dir/seqs_$gg.fna) > $time_pwd/time_pick_otus_uclust_".$taxon.".log 2>&1"; }
        system($cmd);
	print LOG $cmd."\n";
	chdir("picked_otus_uclust_$gg");
	print "pick_rep_set.py\n";
        $cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
	if ($? == -1) {print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: CDHIT ####
sub cdhit {
	my ($gg) = @_;
	print "pick otus: cdhit\n";
	chdir($time_pwd);
	my $dis = $clustThres/100;
        $cmd="pick_otus.py -m cdhit -s $dis -o picked_otus_cdhit_$gg -i splits_dir/seqs_$gg.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd="(time pick_otus.py -m cdhit -s $dis -o picked_otus_cdhit_$gg -i splits_dir/seqs_$gg.fna) > $time_pwd/time_pick_otus_cdhit_".$taxon.".log 2>&1"; }
        system($cmd);
	print LOG $cmd."\n";
	chdir("picked_otus_cdhit_$gg");
	print "pick_rep_set.py\n";
        $cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
	if ($? == -1) {print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: DNACLUST ####
sub dnaclust {
	my ($gg) = @_;
	print "pick otus: dnaclust\n";
	chdir($time_pwd);
	unless (-e "picked_otus_dnaclust_$gg") {
        	system("mkdir picked_otus_dnaclust_$gg");
        }
        if ($? == -1) { print "command failed: $!\n"; exit; }
        chdir("picked_otus_dnaclust_$gg");
        if ($? == -1) { print "command failed: $!\n"; exit; }
	my $dis = $clustThres/100;
        $cmd = "$bin/dnaclust_linux_release3/dnaclust -s $dis -i ../splits_dir/seqs_$gg.fna -t $threads > dnaclustoutfile";
        if ($time eq "y") { $cmd = "(time $bin/dnaclust_linux_release3/dnaclust -s $dis -i ../splits_dir/seqs_$gg.fna -t $threads > dnaclustoutfile) > $time_pwd/time_picked_otus_dnaclust_$taxon.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
        system("$bin/editDnaClustForQiime.pl ../splits_dir/seqs_$gg.fna dnaclustoutfile $gg");
        print LOG "$bin/editDnaClustForQiime.pl ../splits_dir/seqs_$gg.fna dnaclustoutfile $gg\n";
	if ($? == -1) { print "command failed: $!\n"; exit; }
        $cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
	if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
	system($cmd);
	print LOG $cmd."\n";
	if ($? == -1) {print "command failed $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: mothur ####
#sub mothur {
#
#}

#### SUBROUTINE: qiime's mothur ####
sub qiime_mothur {
	my ($linkage, $gg) = @_;
	print "pick otus: mothur\n";
	chdir($time_pwd);
	$cmd = "";
	if ($gg eq "gg2011") {
		$cmd = "align_seqs.py -i splits_dir/seqs_$gg.fna -t $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011_aligned.fasta -o seqs_aligned_$gg/ | tee -a $err";
		if ($time eq "y") { $cmd = "(time align_seqs.py -i splits_dir/seqs_$gg.fna -t $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011_aligned.fasta -o seqs_aligned_$gg/) > $time_pwd/time_pick_otus_mothur_alignseqs_".$taxon."_".$linkage.".log 2>&1"; }
	} elsif ($gg eq "gg2013") {
		$cmd = "align_seqs.py -i splits_dir/seqs_$gg.fna -t $db/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o seqs_aligned_$gg/ | tee -a $err";
		if ($time eq "y") { $cmd = "(time align_seqs.py -i splits_dir/seqs_$gg.fna -t $db/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o seqs_aligned_$gg/) > $time_pwd/time_pick_otus_mothur_alignseqs_".$taxon."_".$linkage.".log 2>&1"; }
	} elsif ($gg eq "silva111") {
		$cmd = "align_seqs.py -i splits_dir/seqs_$gg.fna -t $db/Silva_111_post/rep_set_aligned/97_Silva_111_rep_set.fasta -o seqs_aligned_$gg/ | tee -a $err";
        	if ($time eq "y") { $cmd = "(time align_seqs.py -i splits_dir/seqs_$gg.fna -t $db/Silva_111_post/rep_set_aligned/97_Silva_111_rep_set.fasta -o seqs_aligned_$gg/) > $time_pwd/time_pick_otus_mothur_alignseqs_".$taxon."_".$linkage.".log 2>&1"; }
	}
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
	my $dis = $clustThres/100;
        $cmd = "pick_otus.py -m mothur -s $dis -c $linkage -o picked_otus_mothur_".$linkage."_$gg -i seqs_aligned_$gg/seqs_".$gg."_aligned.fasta | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_otus.py -m mothur -s $dis -c $linkage -o picked_otus_mothur_".$linkage."_$gg -i seqs_aligned_$gg/seqs_".$gg."_aligned.fasta) > $time_pwd/time_pick_otus_mothur_".$taxon."_".$linkage.".log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
        #rename seqs_aligned_otus.txt to seqs_otus.txt
        $cmd = "picked_otus_mothur_".$linkage."_$gg";
        chdir($cmd);
        $cmd = "mv seqs_".$gg."_aligned_otus.txt seqs_".$gg."_otus.txt";
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
        print "pick_rep_set.py\n";
        $cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
	if ($? == -1) {print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: UCLUST-ref ####
sub uclust_ref {
	my ($gg) = @_;
	print "pick otus: uclust-ref\n";
	chdir($time_pwd);
	my $dis = $clustThres/100;
	$cmd = "";
	if ($gg eq "gg2011") {
		$cmd = "pick_otus.py -m uclust_ref -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref_$gg -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta --enable_rev_strand_match | tee -a $err 2>&1";
        	if ($time eq "y") { $cmd = "(time pick_otus.py -m uclust_ref -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref_$gg -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta --enable_rev_strand_match) > $time_pwd/time_pick_otus_uclust-ref_$taxon.log 2>&1"; }
        } elsif ($gg eq "gg2013") {
		$cmd = "pick_otus.py -m uclust_ref -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref_$gg -r $db/gg_13_8_otus/rep_set/97_otus.fasta --enable_rev_strand_match | tee -a $err 2>&1";
		if ($time eq "y") { $cmd = "(time pick_otus.py -m uclust_ref -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref_$gg -r $db/gg_13_8_otus/rep_set/97_otus.fasta --enable_rev_strand_match) > $time_pwd/time_pick_otus_uclust-ref_$taxon.log 2>&1"; }
	} elsif ($gg eq "silva111") {
		$cmd = "pick_otus.py -m uclust_ref -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref_$gg -r $db/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta --enable_rev_strand_match | tee -a $err 2>&1";
		if ($time eq "y") { $cmd = "(time pick_otus.py -m uclust_ref -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref_$gg -r $db/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta --enable_rev_strand_match) > $time_pwd/time_pick_otus_uclust-ref-strict_$taxon.log 2>&1"; }
	}
	system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
	$dir = "picked_otus_uclust-ref_$gg";
	chdir($dir);
	print "pick_rep_set.py\n";
        $cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
	if ($? == -1) {print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: UCLUST-ref-strict ####
sub uclust_ref_strict {
        my ($gg) = @_;
        print "pick otus: uclust-ref-strict\n";
        chdir($time_pwd);
        my $dis = $clustThres/100;
        $cmd = "";
        if ($gg eq "gg2011") {
                $cmd = "pick_otus.py -m uclust_ref -C -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref-strict_$gg -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta --enable_rev_strand_match | tee -a $err 2>&1"; 
                if ($time eq "y") { $cmd = "(time pick_otus.py -m uclust_ref -C -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref-strict_$gg -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta --enable_rev_strand_match) > $time_pwd/time_pick_otus_uclust-ref-strict_".$taxon."_".$gg.".log 2>&1"; }
        } elsif ($gg eq "gg2013") {
                $cmd = "pick_otus.py -m uclust_ref -C -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref-strict_$gg -r $db/gg_13_8_otus/rep_set/97_otus.fasta --enable_rev_strand_match | tee -a $err 2>&1";
                if ($time eq "y") { $cmd = "(time pick_otus.py -m uclust_ref -C -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref-strict_$gg -r $db/gg_13_8_otus/rep_set/97_otus.fasta --enable_rev_strand_match) > $time_pwd/time_pick_otus_uclust-ref-strict_i".$taxon."_".$gg.".log 2>&1"; }
        } elsif ($gg eq "silva111") {
                $cmd = "pick_otus.py -m uclust_ref -C -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref-strict_$gg -r $db/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta --enable_rev_strand_match | tee -a $err 2>&1";  
                if ($time eq "y") { $cmd = "(time pick_otus.py -m uclust_ref -C -s $dis -i splits_dir/seqs_$gg.fna -o picked_otus_uclust-ref-strict_$gg -r $db/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta --enable_rev_strand_match) > $time_pwd/time_pick_otus_uclust-ref-strict_".$taxon."_".$gg.".log 2>&1"; }
        }
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
        $dir = "picked_otus_uclust-ref-strict_$gg";
        chdir($dir);
        print "pick_rep_set.py\n";
        $cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) {print "command failed: $!\n"; exit; }
        chdir('..');
        return;
}

#### SUBROUTINE: BLAST ####
sub blast {
	my ($gg) = @_;
	print "pick otus: BLAST\n";
	chdir($time_pwd);
	$cmd = "";
	if ($gg eq "gg2011") {
		$cmd = "pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_blast_$gg -m blast -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta | tee -a $err 2>&1";
        	if ($time eq "y") { $cmd = "(time pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_blast_$gg -m blast -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta) > $time_pwd/time_pick_otus_blast_$taxon.log 2>&1"; }
	} elsif ($gg eq "gg2013") {
		$cmd = "pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_blast_$gg -m blast -r $db/gg_13_8_otus/rep_set/97_otus.fasta | tee -a $err 2>&1";
		if ($time eq "y") { $cmd = "(time pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_blast_$gg -m blast -r $db/gg_13_8_otus/rep_set/97_otus.fasta) > $time_pwd/time_pick_otus_blast_$taxon.log 2>&1"; }
	} elsif ($gg eq "silva111") {
		$cmd = "pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_blast_$gg -m blast -r $db/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta | tee -a $err 2>&1";
		if ($time eq "y") { $cmd = "(time pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_blast_$gg -m blast -r $db/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta) > $time_pwd/time_pick_otus_blast_$taxon.log 2>&1"; }
	}
        system($cmd);
        print LOG $cmd."\n";
        if ($? == -1) { print "command failed: $!\n"; exit; }
	print "pick_rep_set.py\n";
	$dir = "picked_otus_blast_$gg";
        chdir($dir);
	$cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
        if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
        system($cmd);
        print LOG $cmd."\n";
	if ($? == -1) {print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: swarm ####
sub swarm {
	my ($gg) = @_;
	print "pick otus: swarm\n";
	chdir($time_pwd);
	$cmd = "";
	if ($gg eq "gg2011") {
		$cmd = "pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_swarm_$gg -m swarm -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta | tee -a $err 2>&1     ";
		if ($time eq "y") { $cmd = "(time pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_swarm_$gg -m swarm -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb201     1.fasta) > $time_pwd/time_pick_otus_swarm_$taxon.log 2>&1"; }
	}
	system($cmd);
	print LOG $cmd."\n";
	if ($? == -1) { print "command failed: $!\n"; exit; }
	print "pick_rep_set.py\n";
	$dir = "picked_otus_swarm_$gg";
	chdir($dir);
	$cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
	if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
	system($cmd);
	print LOG $cmd."\n";
	if ($? == -1) { print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: usearch61 ####
sub usearch61 {
	my ($gg) = @_;
	print "pick otus: usearch61\n";
	chdir($time_pwd);
	$cmd = "";
	if ($gg eq "gg2011") {
		$cmd = "pick_otus.py -i splits_dir/seqs_$gg.fna -o picked_otus_usearch61_$gg -m usearch61 -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta | tee -a $err 2>&1";
		if ($time eq "y") { $cmd = "(time pick_otus.py -i splits/dir/seqs_$gg.fna -o picked_otus_usearch61_$gg -m usearch61 -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta) > $time_pwd/time_pick_otus_usearch61_$taxon.log 2>&1"; }
	}
	system($cmd);
	print LOG $cmd."\n";
	if ($? == -1) { print "command failed: $!\n"; exit; }
	print "pick_rep_set.py\n";
	$dir = "picked_otus_usearch61_$gg";
	chdir($dir);
	$cmd = "pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna 2>&1 | tee -a $err";
	if ($time eq "y") { $cmd = "(time pick_rep_set.py -m most_abundant -s otu -i seqs_".$gg."_otus.txt -f ../splits_dir/seqs_$gg.fna -o rep_set.fna) > $time_pwd/time_pick_rep_set.log 2>&1"; }
	system($cmd);
	print LOG $cmd."\n";
	if ($? == -1) { print "command failed: $!\n"; exit; }
	chdir('..');
	return;
}

#### SUBROUTINE: UPARSE ####
sub uparse {
	my ($gg) = @_;
	print "pick otus: UPARSE\n";
	chdir("$time_pwd");
	unless (-e "picked_otus_uparse_$gg") {
        	system("mkdir picked_otus_uparse_$gg");
        }
        if ($? == -1) {print "command failed: $!\n"; exit; }
        chdir("picked_otus_uparse_$gg");
        if ($? == -1) {print "command failed: $!\n"; exit; }
        #$cmd = "$bin/usearch -derep_fulllength ../splits_dir/seqs.fna -output ../splits_dir/unique_seqs.fna  -sizeout -minseqlength 64";
        $cmd = "$bin/usearch -derep_fulllength $time_pwd/splits_dir/seqs.fna -fastaout $time_pwd/splits_dir/unique_seqs.fna -sizeout -threads $threads | tee -a $err";
	if ($time eq "y") {
        	$cmd = "(time $bin/usearch -derep_fulllength $time_pwd/splits_dir/seqs.fna -fastaout $time_pwd/splits_dir/unique_seqs.fna -sizeout -threads $threads) > $time_pwd/time_pick_otus_uparse_".$taxon."0.log 2>&1";
        }
        system($cmd);
        print LOG $cmd."\n";
        #$cmd = "$bin/usearch -sortbysize ../splits_dir/unique_seqs.fna -output ../splits_dir/sorted_unique_seqs.fna -minsize 2";
        $cmd = "$bin/usearch -sortbysize $time_pwd/splits_dir/unique_seqs.fna -fastaout $time_pwd/splits_dir/sorted_unique_seqs.fna -minsize 2 | tee -a $err";
	if ($time eq "y") {
        	$cmd = "(time $bin/usearch -sortbysize $time_pwd/splits_dir/unique_seqs.fna -fastaout $time_pwd/splits_dir/sorted_unique_seqs.fna -minsize 2) > $time_pwd/time_pick_otus_uparse_".$taxon."1.log 2>&1";
        }
        system($cmd);
        print LOG $cmd."\n";
        #$cmd = "$bin/usearch -cluster_otus ../splits_dir/sorted_unique_seqs.fna -otus seqs_otus.otu";
	$cmd = "$bin/usearch -cluster_otus $time_pwd/splits_dir/sorted_unique_seqs.fna -otus $time_pwd/splits_dir/otus1.fa -relabel OTU_ -uparseout results.txt | tee -a $err";
        if ($time eq "y") {
        	$cmd = "(time $bin/usearch -cluster_otus $time_pwd/splits_dir/sorted_unique_seqs.fna -otus $time_pwd/splits_dir/otus1.fa -relabel OTU_ -uparseout results.txt) > $time_pwd/time_pick_otus_uparse_".$taxon."2.log 2>&1";
        }
        system($cmd);
        print LOG $cmd."\n";
        #$cmd = "python $bin/python_scripts/fasta_number.py ../splits_dir/sorted_unique_seqs.fna OTU_ > otus_otu.fna"; #flag!
        #if ($time eq "y") {
       # 	$cmd = "(time python $bin/python_scripts/fasta_number.py ../splits_dir/sorted_unique_seqs.fna OTU_ > otus_otu.fna) > time_pick_otus_uparse_".$taxon."3.log 2>&1";
        #}
        #system($cmd);
        #print LOG $cmd."\n";
        #$cmd = "$bin/usearch -usearch_global ../splits_dir/seqs.fna -db otus_otu.fna -strand plus -id 0.97 -uc map.uc";
	my $clusty = $clustThres/100;
	$cmd = "$bin/usearch -usearch_global $time_pwd/splits_dir/seqs.fna -db $time_pwd/splits_dir/otus1.fa -strand plus -id $clusty -uc map.uc -top_hit_only | tee -a $err";
        if ($time eq "y") {
        	$cmd = "(time $bin/usearch -usearch_global $time_pwd/splits_dir/seqs.fna -db $time_pwd/splits_dir/otus1.fa -strand plus -id $clusty -uc map.uc -top_hit_only | tee -a $err) > $time_pwd/time_pick_otus_uparse_".$taxon."4.log 2>&1";
        }
        system($cmd);
        print LOG $cmd."\n";
	#call parser to parse map.uc into seq_otus.txt
	
	#`perl -i -pe s'/(.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t)/\$1barcodelabel=/g' map.uc`;
        #$cmd = "python $bin/python_scripts/uc2otutab.py map.uc > otu_table.txt";
        #if ($time eq "y") {
        #	$cmd = "(time python $bin/python_scripts/uc2otutab.py map.uc > otu_table.txt) > time_pick_otus_uparse_".$taxon."5.log 2>&1";
        #}
        #system($cmd);
        #print LOG $cmd."\n";
        #$cmd = "biom convert -i otu_table.txt -o otu_table_uparse.biom --table-type=\"otu table\"";
        #system($cmd);
        #print LOG $cmd."\n";
	$cmd = "mv $time_pwd/splits_dir/otus1.fa $time_pwd/picked_otus_uparse_$gg/rep_set.fna";
	system($cmd);
	print LOG $cmd."\n";
	#call UPARSE parser to convert results.txt into seqs_otus.txt
	$cmd = "perl $bin/editUPARSEForQiime.pl results.txt $gg";
	system($cmd);
	print LOG $cmd."\n";
	chdir('..');
	return;
}

#### SUBROUTINE: align_seqs ####
sub align_seqs {
        my ($clust, $gg) = @_;
	print "align_seqs.py\n";
	chdir("$time_pwd/picked_otus_".$clust."_".$gg);
	if (-e "rep_set.fna") {
                if ($gg eq "gg2011") {
                        $cmd = "align_seqs.py -m pynast -a uclust -p 0.75 -i rep_set.fna -t $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011_aligned.fasta -o align_$gg/ -e 100 2>&1 | tee -a $err";
                        if ($time eq "y") { $cmd = "(time align_seqs.py -m pynast -a uclust -p 0.75 -i rep_set.fna -t $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011_aligned.fasta -o align_$gg/ -e 100) > $time_pwd/time_align_seqs.log 2>&1"; }
                } elsif ($gg eq "gg2013") {
                        $cmd = "align_seqs.py -m pynast -a uclust -p 0.75 -i rep_set.fna -t $db/gg_13_8_otus/rep_set_aligned/".$clustThres."_otus.fasta -o align_$gg/ -e 100 2>&1 | tee -a $err";
			if ($time eq "y") { $cmd = "(time align_seqs.py -m pynast -a uclust -p 0.75 -i rep_set.fna -t $db/gg_13_8_otus/rep_set_aligned/".$clustThres."_otus.fasta -o align_$gg/ -e 100) > $time_pwd/time_align_seqs.log 2>&1"; }
                } else { #silva111
			$cmd = "align_seqs.py -m pynast -a uclust -p 0.75 -i rep_set.fna -t $db/Silva_111_post/rep_set_aligned/90_Silva_111_rep_set.fasta -o align_$gg/ -e 100 2>&1 | tee -a $err";
			if ($time eq "y") { $cmd = "(time align_seqs.py -m pynast -a uclust -p 0.75 -i rep_set.fna -t $db/Silva_111_post/rep_set_aligned/90_Silva_111_rep_set.fasta -o align_$gg/ -e 100) > $time_pwd/time_align_seqs.log 2>&1"; }
		}
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "rep_set.fna DNE; aborting align_seqs.py\n";
                print ERR "rep_set.fna DNE; aborting align_seqs.py\n";
        }
	chdir('..');
	return;
}

### SUBROUTINE: chimera_check ###
sub chimera_check {
	my ($gg) = @_;
	print "identify_chimeric_seqs.py\n";
	chdir("$time_pwd/splits_dir");
	#if (-e "align/rep_set_aligned.fasta") {
	if (-e "seqs.fna") {
		if ($gg eq "gg2011") {
			#$cmd = "identify_chimeric_seqs.py -m usearch61 -i align/rep_set_aligned.fasta -r $db/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -o align/chimeric/ 2>&1 | tee -a $err";
			$cmd = "identify_chimeric_seqs.py -m usearch61 --non_chimeras_retention union -i seqs.fna -r $db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta -o chimeric_gg2011/ 2>&1 | tee -a $err";
			#--threads $threads
		} elsif ($gg eq "gg2013") {
			#$cmd = "identify_chimeric_seqs.py -m usearch61 -i align/rep_set_aligned.fasta -r $db/gg_13_8_otus/rep_set/97_otus.fasta -o align/chimeric 2>&1 | tee -a $err";
			$cmd = "identify_chimeric_seqs.py -m usearch61 --non_chimeras_retention union -i seqs.fna -r $db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta -o chimeric_gg2013/ 2>&1 | tee -a $err";
			#--threads $threads
		} else { #silva111
			$cmd = "identify_chimeric_seqs.py -m usearch61 --non_chimeras_retention union -i seqs.fna -r $db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta -o chimeric_silva111/ 2>&1 | tee -a $err";
			#--threads $threads
		}
		system($cmd);
		print LOG $cmd."\n";
		if ($? == -1) { print "command failed: $!\n"; exit; }
		#$cmd = "filter_fasta.py -f align/rep_set_aligned.fasta -o align_rep_set_aligned_nochim.fasta -s align/chimeric/chimeras.txt -n 2>&1 | tee -a $err";
		$cmd = "filter_fasta.py -f seqs.fna -o seqs_nochim.fna -s chimeric_$gg/chimeras.txt -n 2>&1 | tee -a $err";
		system($cmd);
		print LOG $cmd."\n";
		if ($? == -1) { print "command failed: $!\n"; exit; }
		`cp seqs.fna seqs_beforechim.fna`;
		`mv seqs_nochim.fna seqs_$gg.fna`;
	} else {
		#print "align/rep_set_aligned.fasta DNE; aborting identify_chimeric_seqs.py\n";
		#print ERR "align/rep_set_aligned.fasta DNE; aborting identify_chimeric_seqs.py\n";
		print "seqs.fna DNE; aborting identify_chimeric_seqs.py\n";
		print ERR "seqs.fna DNE; aborting identify_chimeric_seqs.py\n";
	}
	chdir('..');
	return;
}

#### SUBROUTINE: no_chimera_check ####
sub no_chimera_check {
	my ($gg) = @_;
	chdir("$time_pwd/splits_dir");
	if (-e "seqs.fna") {
		`mv seqs.fna seqs_$gg.fna`;
	} else {
		print "seqs.fna DNE; aborting no_chimera_check\n";
		print ERR "seqs.fna DNE; aborting no_chimera_check\n";
	}
	chdir('..');
	return;
}

#### SUBROUTINE: taxa-blast ####
sub taxa_blast() {
	my ($clust, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	print "assign_taxonomy.py\n";
        if (-e "rep_set.fna") {
		if ($gg eq "gg2011") {
        		$cmd="assign_taxonomy.py -m blast -i rep_set.fna -r $db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta -t $db/gg_otus_4feb2011/taxonomies/greengenes_tax.txt -o blast_".$gg."_assigned_taxonomy_tr 2>&1 | tee -a $err";
                	if ($time eq "y") { $cmd="(time assign_taxonomy.py -m blast  -i rep_set.fna -r $db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta -t $db/gg_otus_4feb2011/taxonomies/greengenes_tax.txt -o blast_".$gg."_assigned_taxonomy_tr) > $time_pwd/time_assign_taxonomy_blast.log 2>&1"; }
		} elsif ($gg eq "gg2013") {
			$cmd="assign_taxonomy.py -m blast -i rep_set.fna -r $db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta -t $db/gg_13_8_otus/taxonomy/".$clustThres."_otu_taxonomy.txt -o blast_".$gg."_assigned_taxonomy_tr 2>&1 | tee -a $err";
			if ($time eq "y") { $cmd="(time assign_taxonomy.py -m blast  -i rep_set.fna -r $db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta -t $db/gg_13_8_otus/taxonomy/".$clustThres."_otu_taxonomy.txt -o blast_".$gg."_assigned_taxonomy_tr) > $time_pwd/time_assign_taxonomy_blast.log 2>&1"; }
		} elsif ($gg eq "silva111") {
			$cmd="assign_taxonomy.py -m blast -i rep_set.fna -r $db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta -t $db/Silva_111_post/taxonomy/".$clustThres."_Silva_111_taxa_map_RDP_6_levels.txt -o blast_".$gg."_assigned_taxonomy_tr 2>&1 | tee -a $err";
			if ($time eq "y") { $cmd="(time assign_taxonomy.py -m blast  -i rep_set.fna -r $db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta -t $db/Silva_111_post/taxonomy/".$clustThres."_Silva_111_taxa_map_RDP_6_levels.txt -o blast_".$gg."_assigned_taxonomy_tr) > $time_pwd/time_assign_taxonomy_blast.log 2>&1"; }
		}
		system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
	} else {
                print "rep_set.fna DNE; aborting assign_taxonomy.py\n";
                print ERR "rep_set.fna DNE; aborting assign_taxonomy.py\n";
        }
	chdir('..');
	return;
}

####SUBROUTINE: taxa-rdp-training ####
sub taxa_rdp_training() {
	my ($clust, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	print "assign_taxonomy.py\n";
	if (-e "rep_set.fna") {
		if ($gg eq "gg2011") {
                	$cmd="assign_taxonomy.py -m rdp -c $taxon_cutoff --rdp_max_memory 15000 -i rep_set.fna -r $db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta -t $db/gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train_genus.txt -o rdp-training_".$gg."_assigned_taxonomy_tr 2>&1 | tee -a $err";
                        if ($time eq "y") { $cmd="(time assign_taxonomy.py -m rdp -c $taxon_cutoff --rdp_max_memory 8000 -i rep_set.fna -r $db/gg_otus_4feb2011/rep_set/gg_".$clustThres."_otus_4feb2011.fasta -t $db/gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train_genus.txt -o rdp-training_".$gg."_assigned_taxonomy_tr) > $time_pwd/time_assign_taxonomy_rdp-training.log 2>&1"; }
                } elsif ($gg eq "gg2013") {
                        $cmd="assign_taxonomy.py -m rdp -c $taxon_cutoff --rdp_max_memory 8000 -i rep_set.fna -r $db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta -t $db/gg_13_8_otus/taxonomy/".$clustThres."_otu_taxonomy.txt -o rdp-training_".$gg."_assigned_taxonomy_tr 2>&1 | tee -a $err";
			if ($time eq "y") { $cmd="(time assign_taxonomy.py -m rdp -c $taxon_cutoff --rdp_max_memory 8000 -i rep_set.fna -r $db/gg_13_8_otus/rep_set/".$clustThres."_otus.fasta -t $db/gg_13_8_otus/taxonomy/".$clustThres."_otu_taxonomy.txt -o rdp-training_".$gg."_assigned_taxonomy_tr) > $time_pwd/time_assign_taxonomy_rdp-training.log 2>&1"; }
                } elsif ($gg eq "silva111") {
			$cmd="assign_taxonomy.py -m rdp -c $taxon_cutoff --rdp_max_memory 8000 -i rep_set.fna -r $db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta -t $db/Silva_111_post/taxonomy/".$clustThres."_Silva_111_taxa_map_RDP_6_levels.txt -o rdp-training_".$gg."_assigned_taxonomy_tr 2>&1 | tee -a $err";
			if ($time eq "y") { $cmd="(time assign_taxonomy.py -m rdp -c $taxon_cutoff --rdp_max_memory 8000 -i rep_set.fna -r $db/Silva_111_post/rep_set/".$clustThres."_Silva_111_rep_set.fasta -t $db/Silva_111_post/taxonomy/".     $clustThres."_Silva_111_taxa_map_RDP_6_levels.txt -o rdp-training_".$gg."_assigned_taxonomy_tr) > $time_pwd/time_assign_taxonomy_rdp-training.log 2>&1"; }
		}
		system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }	
	} else {
		print "rep_set.fna DNE; aborting assign_taxonomy.py\n";
                print ERR "rep_set.fna DNE; aborting assign_taxonomy.py\n";
	}
	chdir('..');
	return;
}

#### SUBROUTINE: make_otu_table ####
sub make_otu_table {
	my ($clust, $taxon, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
        print "make_otu_table.py\n";
        if (!-z "seqs_".$gg."_otus.txt") {
                $cmd="make_otu_table.py -i seqs_".$gg."_otus.txt -t ".$taxon."_".$gg."_assigned_taxonomy_tr/rep_set_tax_assignments.txt -o otu_table_".$taxon."_".$gg.".biom 2>&1 | tee -a $err";
                if ($time eq "y") { $cmd="(time make_otu_table.py -i seqs_".$gg."_otus.txt -t ".$taxon."_".$gg."_assigned_taxonomy_tr/rep_set_tax_assignments.txt -o otu_table_".$taxon."_".$gg.".biom) > $time_pwd/time_make_otu_table_".$taxon.".log 2>&1"; }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
		#if chimera checking
#		if ($chimeraChecking eq "y") {
#			$cmd = "make_otu_table.py -i seqs_otus.txt -t ".$taxon."_assigned_taxonomy_tr/rep_set_tax_assignments.txt -o otu_table_".$taxon."_nochim.biom -e align/chimeric/chimeras.txt 2>&1 | tee -a $err";
#			system($cmd);
#			print LOG $cmd."\n";
#			if ($? == -1) { print "command failed: $!\n"; exit; }
#		}
        } else {
                print "seqs_".$gg."_otus.txt DNE or is empty; aborting make_otu_table.py\n";
                print ERR "seqs_".$gg."_otus.txt DNE or is empty; aborting make_otu_table.py\n";
        }
	chdir('..');
}

#### SUBROUTINE: filter_otus ####
sub filter_otus {
	my ($clust, $taxon, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
        print "filter_otus_from_otu_table\n";
        if (-e "otu_table_".$taxon."_".$gg.".biom") {
                $cmd = "filter_otus_from_otu_table.py -i otu_table_".$taxon."_".$gg.".biom -o otu_table_".$taxon."_".$gg."_n1.biom -n 2 2>&1 | tee -a $err";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
                $cmd = "filter_otus_from_otu_table.py -i otu_table_".$taxon."_".$gg.".biom -o otu_table_".$taxon."_".$gg."_n2.biom -n 3 2>&1 | tee -a $err";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) {print "command failed: $!\n"; exit; }
#		if ($chimeraChecking eq "y") {
#			$cmd = "filter_otus_from_otu_table.py -i otu_table_".$taxon."_nochim.biom -o otu_table_".$taxon."_nochim_n1.biom -n 2 2>&1 | tee -a $err";
#			system($cmd);
#			print LOG $cmd."\n";
#			if ($? == -1) { print "command failed: $!\n"; exit; }
#			$cmd = "filter_otus_from_otu_table.py -i otu_table_".$taxon."_nochim.biom -o otu_table_".$taxon."_nochim_n2.biom -n 3 2>&1 | tee -a $err";
#			system($cmd);
#			print LOG $cmd."\n";
#			if ($? == -1) { print "command failed: $!\n"; exit; }
#		}
        } else {
                print "otu_table_".$taxon."_".$gg.".biom DNE; aborting filter_otus_from_otu_table.py\n";
                print ERR "otu_table_".$taxon."_".$gg.".biom DNE; aborting filter_otus_from_otu_table.py\n";
        }
	#convert_biom
        print "convert_biom.py | biom convert\n";
	my $convert_biom; my $biom_ver; my $biom_ver2;
        if (-e "otu_table_".$taxon."_".$gg."_n1.biom") {
                $convert_biom = `which biom`;
		$biom_ver = `biom --version 2>&1`;
		#$biom_ver2 = `biom convert --version`;
		chomp($biom_ver);
		#$biom_ver=~/biom, version (.*)/;
		#$biom_ver = $1; $biom_ver=~s/\.//g;
                if (!$convert_biom) {
                        $cmd = "convert_biom.py -i otu_table_".$taxon."_".$gg."_n1.biom -o otu_table_".$taxon."_".$gg."_n1.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } elsif ($biom_ver!~/biom,\s*version.*/) {
                        $cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n1.biom -o otu_table_".$taxon."_".$gg."_n1.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                #} elsif ($biom_ver2 =~/Version: biom convert 2\.1\.4/) {
		#	$cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n1.biom -o otu_table_".$taxon."_".$gg."_n1.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		} else {
			$cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n1.biom -o otu_table_".$taxon."_".$gg."_n1.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		}
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg."_n1.biom DNE; aborting biom convert\n";
                print ERR "otu_table_".$taxon."_".$gg."_n1.biom DNE; aborting biom convert\n";
        }
        if (-e "otu_table_".$taxon."_".$gg.".biom") {
                if (!$convert_biom) {
                        $cmd = "convert_biom.py -i otu_table_".$taxon."_".$gg.".biom -o otu_table_".$taxon."_".$gg.".txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } elsif ($biom_ver!~/biom,\s*version.*/) {
                $cmd ="biom convert -i otu_table_".$taxon."_".$gg.".biom -o otu_table_".$taxon."_".$gg.".txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		#} elsif ($biom_ver2 =~/Version: biom convert 2\.1\.4/) {
                #        $cmd ="biom convert -i otu_table_".$taxon."_".$gg.".biom -o otu_table_".$taxon."_".$gg.".txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } else {
			$cmd ="biom convert -i otu_table_".$taxon."_".$gg.".biom -o otu_table_".$taxon."_".$gg.".txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		}
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg.".biom DNE; aborting biom conert\n";
                print ERR "otu_table_".$taxon."_".$gg.".biom DNE; aborting biom conert\n";
        }
        if (-e "otu_table_".$taxon."_".$gg."_n2.biom") {
                if (!$convert_biom) {
                        $cmd = "convert_biom.py -i otu_table_".$taxon."_".$gg."_n2.biom -o otu_table_".$taxon."_".$gg."_n2.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } elsif ($biom_ver!~/biom,\s*version.*/) {
                        $cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n2.biom -o otu_table_".$taxon."_".$gg."_n2.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		#} elsif ($biom_ver2 =~/Version: biom convert 2\.1\.4/) {
                #        $cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n2.biom -o otu_table_".$taxon."_".$gg."_n2.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } else {
			$cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n2.biom -o otu_table_".$taxon."_".$gg."_n2.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		}
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
	} else {
                print "otu_table_".$taxon."_".$gg."_n2.biom DNE; aborting biom convert\n";
                print ERR "otu_table_".$taxon."_".$gg."_n2.biom DNE; aborting biom convert\n";
        }
	chdir('..');
}

#### SUBROUTINE: filter_alignment ####
sub filter_alignment() {
        my ($clust, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	print "filter_alignment.py\n";
        if (-e "align_$gg/rep_set_aligned.fasta") {
                $cmd = "filter_alignment.py -i align_$gg/rep_set_aligned.fasta -o align_".$gg."_filtered/ 2>&1 | tee -a $err";
                if ($time eq "y") { $cmd = "(time filter_alignment.py -i align_$gg/rep_set_aligned.fasta -o align_".$gg."_filtered/) > $time_pwd/time_filter_alignment.log 2>&1"; }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "align_$gg/rep_set_aligned.fasta DNE; aborting filter_alignment.py\n";
                print ERR "align_$gg/rep_set_aligned.fasta DNE; aborting filter_alignment.py\n";
        }
	chdir('..');
}

#### SUBROUTINE: make_phylogeny ####
sub make_phylogeny() {
        my ($clust, $gg) = @_;
	chdir("picked_otus_".$clust."_$gg");
	print "make_phylogeny.py\n";
        if (-e "align_".$gg."_filtered/rep_set_aligned_pfiltered.fasta") {
                $cmd = "make_phylogeny.py -i align_".$gg."_filtered/rep_set_aligned_pfiltered.fasta -r midpoint -o rep_phylo_$gg.tre 2>&1 | tee -a $err";
                if ($time eq "y") { $cmd = "(time make_phylogeny.py -i align_".$gg."_filtered/rep_set_aligned_pfiltered.fasta -r midpoint -o rep_phylo_$gg.tre) > $time_pwd/time_make_phylogeny.log 2>&1"; }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "align_".$gg."_filtered/rep_set_aligned_pfiltered.fasta DNE; aborting make_phylogeny\n";
                print ERR "align_".$gg."_filtered/rep_set_aligned_pfiltered.fasta DNE; aborting make_phylogeny\n";
        }
	chdir('..');
}

#### SUBROUTINE: make_ggpruned ####
sub make_ggpruned() {
        my ($clust, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	print "make GG pruned phylogeny\n";
        if (-e "align_$gg/rep_set_log.txt") {
		$cmd = "cut -f 4 align_$gg/rep_set_log.txt > gg_taxa.tmp";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
		#remove header line
                $cmd = "perl -i -pe \"s/template ID\n//g\" gg_taxa.tmp";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
		# doublecheck; get rid of anything that's not a digit
                $cmd = "perl -i -pe 's/\\D+\n//g' gg_taxa.tmp";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
		# get rid of any blank lines (no search results for the OTU in the db)
                $cmd = "perl -i -pe \"s/^\n//g\" gg_taxa.tmp";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
                $cmd = "perl -i -pe \"s/^\\D+.*?\n//g\" gg_taxa.tmp";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
                $cmd = "perl -i -pe \"s/^\\D+.*?\$//g\" gg_taxa.tmp";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
                $cmd = "perl -i -pe \"s/.*:.*\n//g\" gg_taxa.tmp";
                print "filter_tree.py\n";
                $cmd = "";
		my $prunedtree = "";
		if (($gg eq "gg2011") || ($gg eq "all")) {
			$cmd = "filter_tree.py -i $db/gg_otus_4feb2011/trees/gg_".$clustThres."_otus_4feb2011.tre -t gg_taxa.tmp -o gg2011_".$clustThres."_pruned.tre 2>&1 | tee -a $err";
			system($cmd);
			print LOG $cmd."\n";
			if ($? == -1) { print "command failed: $!\n"; exit; }
			$prunedtree = "gg2011_".$clustThres."_pruned.tre";
		} elsif (($gg eq "gg2013") || ($gg eq "all")) {
			$cmd = "filter_tree.py -i $db/gg_13_8_otus/trees/".$clustThres."_otus.tree -t gg_taxa.tmp -o gg2013_".$clustThres."_pruned.tre 2>&1 | tee -a $err";
			system($cmd);
			print LOG $cmd."\n";
			if ($? == -1) { print "command failed: $!\n"; exit; }
			$prunedtree = "gg2013_".$clustThres."_pruned.tre";
		} elsif (($gg eq "silva111") || ($gg eq "all")) {
			$cmd = ""; #There is no accompanying silv111 phylogeny to filter; thus gg_pruned tree is not created.
			print LOG "filter_tree.py: `There is no accompany silva111 phylogeny to filter; thus gg_pruned tree is not created.\n";
			return;
		}
                #`cp gg_97_pruned.tre gg_97_pruned.bck.tre`;
                #rename branches on the gg_97_pruned.tre to OTU numbers
                open (FIN, "<", "gg_taxa.tmp") or die "$!\n";
                my @fin = <FIN>;
                $i = 0;
                my $woline = $fin[$i];
                while ($i <= $#fin) {
                        chop($woline);
                        #given greengenes id, retrieve otu
                        #won't necessarily be unique, retrieve every line in rep_set_log with $woline
                        #$linematches = `strings align_$gg/rep_set_log.txt | egrep '.*?\\s$woline\\s.*?'`; This is not specific to finding a taxa id but instead anything in the line that matches the taxa id for e.g. as in the length of the sequence input!
			#the column in rep_set_log.txt depends on whether OTUs were picked reference or open
			my $var = "";
			if ($clust eq "uparse") {
				$var = "3";
			} elsif (($clust eq "uclust") | ($clust eq "cdhit") | ($clust eq "dnaclust") | ($clust eq "uclust-ref") | ($clust eq "uclust-ref-strict") | ($clust eq "uparse")) {
				$var = "4";
			} else {
				$var = "5";
			}
			my $linematches = `strings align_$gg/rep_set_log.txt | awk '{if (\$$var == "$woline") print \$0;}'`;
                        my @matches = split('\n', $linematches);
			my $otus = "";
			$matches[0]=~/(.*?)\s.*/;
			if ($#matches+1 == 1) {
				if ($clust eq "uparse") {
					$otus = "'\\''$1'\\''";
				} elsif ($clust eq "uclust") {
					$otus = "$1";
				} else {
					$otus = "OTU".$1."";
				}
			} else {
                        	#$matches[0]=~/(.*?)\s.*/;
				if ($clust eq "uparse") {
					$otus = "('\\''$1'\\'':0";
				} elsif ($clust eq "uclust") {
					$otus = "($1:0";
				} else {
                        		$otus = "(OTU".$1.":0";
				}
                        	for($a=1; $a <=$#matches; $a++) {
                        	        $matches[$a]=~/(.*?)\s.*/;
					if ($clust eq "uparse") {
						$otus = $otus.",'\\''".$1."'\\'':0";
					} elsif ($clust eq "uclust") {
						$otus = $otus.",".$1.":0";
					} else {
                        	        	$otus = $otus.","."OTU".$1.":0";
					}
                        	}
                        	$otus = $otus.")";
			}
                        $cmd = "perl -i -pe 's/([;,(])$woline:/\\1$otus:/g' $prunedtree";
                        system($cmd);
                        #$cmd = "perl -i -pe 's/,$woline:/,$otus:/g' gg_97_pruned.tre";
                        #print $cmd."\n";
                        #system($cmd);
                        print LOG $cmd."\n";
                        if ($? == -1) { print "command failed: $!\n"; exit; }
                        $i += 1;
                        $woline = $fin[$i];
                }
		close FIN;
		#Remove OTUs from phylogeny
		if ($clust ne "uparse") {
			$cmd = "perl -i -pe 's/OTU//g' $prunedtree";
			system($cmd);
		}
		$cmd = "rm gg_taxa.tmp";
		system($cmd);
		print LOG $cmd."\n";
		if ($> == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "align_$gg/rep_set_log.txt DNE; aborting make GG pruned phylogeny\n";
                print ERR "align_$gg/rep_set_log.txt DNE; aborting make GG pruned phylogeny\n";
        }
	chdir('..');
}

#### SUBROUTINE: removeRoot ####
sub removeRoot() {
        my ($clust, $taxon, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	$pwd=`pwd`;
	chomp($pwd);
	print "removeRoot.pl\n";
        #call removeRoot on n1 biom file
	if ((-e "otu_table_".$taxon."_".$gg."_n1.biom") && ("otu_table_".$taxon."_".$gg."_n1.txt")) {
                $cmd = "$bin/removeRoot.pl ".$pwd."/otu_table_".$taxon."_".$gg."_n1 ".$pwd."/otu_table_".$taxon."_".$gg."_n1_noRoot 2>&1 | tee -a $err";
                system($cmd);
                if ($? == -1) { print "command failed: $!\n"; exit; }
                print LOG $cmd."\n";
        } else {
                print "otu_table_".$taxon."_".$gg."_n1.biom or otu_table_".$taxon."_".$gg."_n1.txt DNE; aborting removeRoot.pl\n";
                print ERR "otu_table_".$taxon."_".$gg."_n1.biom or otu_table_".$taxon."_".$gg."_n1.txt DNE; aborting removeRoot.pl\n";
        }
        #call removeRoot on original biom file
        if ((-e "otu_table_".$taxon."_".$gg.".biom") && ("otu_table_".$taxon."_".$gg.".txt")) {
                $cmd = "$bin/removeRoot.pl ".$pwd."/otu_table_".$taxon."_".$gg." ".$pwd."/otu_table_".$taxon."_".$gg."_noRoot 2>&1 | tee -a $err";
                system($cmd);
                print LOG $cmd."\n";
		if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg.".biom or otu_table_".$taxon."_".$gg.".txt DNE; aborting removeRoot.pl\n";
                print ERR "otu_table_".$taxon."_".$gg.".biom or otu_table_".$taxon."_".$gg.".txt DNE; aborting removeRoot.pl\n";
        }
        #call removeRoot on n2 biom file
        if ((-e "otu_table_".$taxon."_".$gg."_n2.biom") && ("otu_table_".$taxon."_".$gg."_n2.txt")) {
                $cmd = "$bin/removeRoot.pl ".$pwd."/otu_table_".$taxon."_".$gg."_n2 ".$pwd."/otu_table_".$taxon."_".$gg."_n2_noRoot 2>&1 | tee -a $err";
                system($cmd);
                if ($? == -1) {print "command failed: $!\n"; exit; }
                print LOG $cmd."\n";
        } else {
                print "otu_table_".$taxon."_".$gg."_n2.biom or otu_table_".$taxon."_".$gg."_n2.txt DNE; aborting removeRoot.pl\n";
                print ERR "otu_table_".$taxon."_".$gg."_n2.biom or otu_table_".$taxon."_".$gg."_n2.txt DNE; aborting removeRoot.pl\n";
        }
	#convert_biom
	my $convert_biom; my $biom_ver; my $biom_ver2;
	$convert_biom = `which biom`;
	$biom_ver = `biom --version 2>&1`;
	#$biom_ver2 = `biom convert --version`;
        print "convert_biom.py | biom convert\n";
	if (-e "otu_table_".$taxon."_".$gg."_n2_noRoot.biom") {
		if (!$convert_biom) {
			$cmd = "convert_biom.py -i otu_table_".$taxon."_".$gg."_n2_noRoot.biom -o otu_table_".$taxon."_".$gg."_n2_noRoot.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		} elsif ($biom_ver!~/biom,\s*version.*/) {
			$cmd = "biom convert -i otu_table_".$taxon."_".$gg."_n2_noRoot.biom -o otu_table_".$taxon."_".$gg."_n2_noRoot.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		#} elsif ($biom_ver2 =~/Version: biom convert 2\.1\.4/) {
                #        $cmd = "biom convert -i otu_table_".$taxon."_".$gg."_n2_noRoot.biom -o otu_table_".$taxon."_".$gg."_n2_noRoot.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		} else {
			$cmd = "biom convert -i otu_table_".$taxon."_".$gg."_n2_noRoot.biom -o otu_table_".$taxon."_".$gg."_n2_noRoot.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		}
		system($cmd);
		if ($? == -1) { print "command failed: $!\n"; exit; }
		print LOG $cmd."\n";
	} else {
		print "otu_table_".$taxon."_".$gg."_n2_noRoot.biom DNE; aborting biom convert\n";
		print ERR "otu_table_".$taxon."_".$gg."_n2_noRoot.biom DNE; aborting biom convert\n";
	}
        if (-e "otu_table_".$taxon."_".$gg."_n1_noRoot.biom") {
                if (!$convert_biom) {
                        $cmd = "convert_biom.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -o otu_table_".$taxon."_".$gg."_n1_noRoot.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } elsif ($biom_ver!~/biom,\s*version.*/) {
                        $cmd = "biom convert -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -o otu_table_".$taxon."_".$gg."_n1_noRoot.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		#} elsif ($biom_ver2 =~/Version: biom convert 2\.1\.4/) {
                #        $cmd ="biom convert -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -o otu_table_".$taxon."_".$gg."_n1_noRoot.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } else {
			$cmd = "biom convert -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -o otu_table_".$taxon."_".$gg."_n1_noRoot.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		}
                system($cmd);
                if ($? == -1) { print "command failed: $!\n"; exit; }
                print LOG $cmd."\n";
        } else {
                print "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting biom convert\n";
                print ERR "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting biom convert\n";
        }
        if (-e "otu_table_".$taxon."_".$gg."_noRoot.biom") {
                if (!$convert_biom) {
                        $cmd = "convert_biom.py -i otu_table_".$taxon."_".$gg."_noRoot.biom -o otu_table_".$taxon."_".$gg."_noRoot.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } elsif ($biom_ver!~/biom,\s*version.*/) {
                        $cmd = "biom convert -i otu_table_".$taxon."_".$gg."_noRoot.biom -o otu_table_".$taxon."_".$gg."_noRoot.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		#} elsif ($biom_ver2 =~/Version: biom convert 2\.1\.4/) {
                #        $cmd ="biom convert -i otu_table_".$taxon."_".$gg."_noRoot.biom -o otu_table_".$taxon."_".$gg."_noRoot.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
                } else {
			$cmd = "biom convert -i otu_table_".$taxon."_".$gg."_noRoot.biom -o otu_table_".$taxon."_".$gg."_noRoot.txt --to-tsv --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
		}
                system($cmd);
                if ($? == -1) { print "command failed: $!\n"; exit; }
                print LOG $cmd."\n";
        } else {
                print "otu_table_".$taxon."_".$gg."_noRoot.biom DNE; aborting biom convert\n";
                print ERR "otu_table_".$taxon."_".$gg."_noRoot.biom DNE; aborting biom convert\n";
        }
#	if (-e "otu_table_".$taxon."_chim_noRoot.biom") {
#		if (!$convert_biom) {
#			$cmd = "convert_biom.py -i otu_table_".$taxon."_chim_noRoot.biom -o otu_table_".$taxon."_chim_noRoot.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
#		} else {
#			$cmd = "biom convert -i otu_table_".$taxon."_chim_noRoot.biom -o otu_table_".$taxon."_chim_noRoot.xt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
#		}
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
#	if (-e "otu_table_".$taxon."_chim_n1_noRoot.biom") {
#		if (!$convert_biom) {
#			$cmd = "convert_biom.py -i otu_table_".$taxon."_chim_n1_noRoot.biom -o otu_table_".$taxon."_chim_n1_noRoot.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
#		} else {
#			$cmd = "biom convert -i otu_table_".$taxon."_chim_n1_noRoot.biom -o otu_table_".$taxon."_chim_n1_noRoot.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
#		}
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
#	if (-e "otu_table_".$taxon."_chim_n2_noRoot.biom") {
#		if (!$convert_biom) {
#			$cmd = "convert_biom.py -i otu_table_".$taxon."_chim_n2_noRoot.biom -o otu_table_".$taxon."_chim_n2_noRoot.txt -b --header_key=taxonomy --output_metadata_id=\"Consensus Lineage\" 2>&1 | tee -a $err";
#		} else {
#			$cmd = "biom convert -i otu_table_".$taxon."_chim_n2_noRoot.biom -o otu_table_".$taxon."_chim_n2_noRoot.txt -b --header-key=taxonomy --output-metadata-id=\"Consensus Lineage\" 2>&1 | tee -a $err";
#		}
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
	#Create a relative abundance table
	$cmd = "Rscript $bin/sl1p_relAbund.R -o otu_table_".$taxon."_".$gg."_n1_noRoot.txt -m ../$mapf -r otu_table_".$taxon."_".$gg."_n1_noRoot_relAbund.txt";
	print $cmd."\n";
	system($cmd);
	$cmd = "Rscript $bin/sl1p_relAbund.R -o otu_table_".$taxon."_".$gg.".txt -m ../$mapf -r otu_table_".$taxon."_".$gg."_relAbund.txt";
	#Organize files into otu_aux folder
	mkdir('otus_aux');
	$cmd = "mv otu_table_".$taxon."_".$gg."_n2* otus_aux";
	system($cmd);
	$cmd = "mv otu_table_".$taxon."_".$gg."_n1.* otus_aux";
	system($cmd);
	$cmd = "mv otu_table_".$taxon."_".$gg."_noRoot.* otus_aux";
	system($cmd);
	chdir('..');
}

#### SUBROUTINE: summarize_taxa ####
sub summarize_taxa() {
        print "summarize_taxa_through_plots.py\n";
        my ($clust, $taxon, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	if (-e "otu_table_".$taxon."_".$gg.".biom") {
                $cmd = "summarize_taxa_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -f -o wf_taxa_summary_".$taxon."_".$gg." -m ../$mapf -p ../qiime_params.txt 2>&1 | tee -a $err";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg.".biom DNE; aborting summarize_taxa_through_plots.py\n";
                print ERR "otu_table_".$taxon."_".$gg.".biom DNE; aborting summarize_taxa_through_plots.py\n";
        }
        if (-e "otu_table_".$taxon."_".$gg."_n1_noRoot.biom") {
                $cmd = "summarize_taxa_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -f -o wf_taxa_summary_".$taxon."_".$gg."_n1_noRoot -m ../$mapf -p ../qiime_params.txt 2>&1 | tee -a $err";
                if ($time eq "y") { $cmd = "(time summarize_taxa_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -f -o wf_taxa_summary_".$taxon."_".$gg."_n1_noRoot/ -m ../$mapf -p ../qiime_params.txt) > $time_pwd/time_sum_taxa_thru_plots_".$taxon.".log 2>&1"; }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting summarize_taxa_through_plots.py\n";
                print ERR "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting summarize_taxa_through_plots.py\n";
	}
#	if (-e "otu_table_".$taxon."_chim.biom") {
#		$cmd = "summarize_taxa_through_plots.py -i otu_table_".$taxon."_chim.biom -o wf_taxa_summary_".$taxon."_chim -m ../$mapf -p ../qiime_params.txt 2>&1 | tee -a $err";
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
#	if (-e "otu_table_".$taxon."_chim_n1_noRoot.biom") {
#		$cmd = "summarize_taxa_through_plots.py -i otu_table_".$taxon."_chim_n1_noRoot.biom -o wf_taxa_summary_".$taxon."_chim_n1_noRoot -m ../$mapf -p ../qiime_params.txt 2>&1 | tee -a $err";
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
	chdir('..');
}

#### SUBROUTINE: alpha_div ####
sub alpha_div() {
        print "alpha_rarefaction.py\n";
	my ($clust, $taxon, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
        if (-e "otu_table_".$taxon."_".$gg.".biom") {
                $cmd = "alpha_rarefaction.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_arare_".$taxon."_".$gg."/ -p ../qiime_params.txt -t rep_phylo_".$gg.".tre -a -O $threads 2>&1 | tee -a $err";
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) {print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg.".biom DNE; aborting alpha_rarefaction.py\n";
                print ERR "otu_table_".$taxon."_".$gg.".biom DNE; aborting alpha_rarefaction.py\n";
        }
        if (-e "otu_table_".$taxon."_".$gg."_n1_noRoot.biom") {
                $cmd = "alpha_rarefaction.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -f -o wf_arare_".$taxon."_".$gg."_n1_noRoot/ -p ../qiime_params.txt -t rep_phylo_".$gg.".tre -a -O $threads 2>&1 | tee -a $err";
                if ($time eq "y") { $cmd = "(time alpha_rarefaction.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -f -o wf_arare_".$taxon."_".$gg."_n1_noRoot/ -p ../qiime_params.txt -t rep_phylo_".$gg.".tre -a -O $threads) > $time_pwd/time_alpha_raref.log 2>&1"; }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) {print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting alpha_rarefaction.py\n";
                print ERR "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting alpha_rarefaction.py\n";
        }
#	if (-e "otu_table_".$taxon."_chim.biom") {
#		$cmd = "alpha_rarefaction.py -i otu_table_".$taxon."_chim.biom -m ../".$mapf." -o wf_arare_".$taxon."_chim/ -p ../qiime_params.txt -t rep_phylo_chim.tre -a -O $threads 2>&1 | tee -a $err";
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
#	if (-e "otu_table_".$taxon."_chim_n1_noRoot.biom") {
#		$cmd = "alpha_rarefaction.py -i otu_table_".$taxon."_chim_n1_noRoot.biom -m ..$mapf -o wf_arare_".$taxon."_chim_n1_noroot/ -p ../qiime_params.txt -t rep_phylo_chim.tre -a -O $threads 2>&1 | tee -a $err";
#		system($cmd);
 #               if ($? == -1) { print "command failed: $!\n"; exit; }
  #              print LOG $cmd."\n";
#	}
	chdir('..');
}

#### SUBROUTINE: beta_div ####
sub beta_div() {
	my ($clust, $taxon, $gg) = @_;
	chdir("$time_pwd/picked_otus_".$clust."_$gg");
	my $otu = "otu_table_".$taxon."_".$gg.".biom";
	my $lib_stats; my $minsamples; my @minsplit; my $evalue; my $evalue1;
        if (-e $otu) {
                $lib_stats = `which biom`;
		$out = "library_stats_".$taxon."_$gg.txt";
                if ($lib_stats) {
                        print "biom summarize-table\n";
                        $cmd = "biom summarize-table -i $otu -o $out 2>&1 | tee -a $err";
			print LOG $cmd."\n";
                        system($cmd);
                        $minsamples = `grep 'Min:' $out`;
                        @minsplit = split (' ', $minsamples);
                        $evalue = $minsplit[1];
                        $evalue =~ /(.*)\..*/;
                        $evalue = $1;
			$evalue =~ s/,//g;
                        print "evalue: $evalue \n";
                        #$cmd = "biom summarize-table -i $otu -o library_stats.txt 2>&1 | tee -a $err";
                } else {
                        print "per_library_stats.py\n";
                        $evalue = `per_library_stats.py -i $otu | awk '/core_qiime_analyses.py \\(just a suggestion\\): /{print \$5}'`;
                        chomp($evalue);
                        $evalue=~ /(.*)\..*/;
                        $evalue = $1;
			$evalue =~ s/,//g;
                        print "evalue: $evalue \n";
                        $cmd = "per_library_stats.py -i $otu > $out";
                        system($cmd);
                        if ($? == -1) { print "command failed: $!\n"; exit; }
                        print LOG $cmd."\n";
                }
        } else {
                print "$otu DNE; aborting biom summarize-table\n";
                print ERR "$otu DNE; aborting biom summarize-table\n";
        }
        $otu = "otu_table_".$taxon."_".$gg."_n1_noRoot.biom";
	if (-e $otu) {
                $lib_stats = `which biom`;
		$out = "library_stats_".$taxon."_".$gg."_n1_noRoot.txt";
                if ($lib_stats) {
                        print "biom summarize-table\n";
                        $cmd = "biom summarize-table -i $otu -o $out 2>&1 | tee -a $err";
                        system($cmd);
                        print LOG $cmd."\n";
			$minsamples = `grep 'Min:' $out`;
                        @minsplit = split (' ', $minsamples);
                        $evalue1 = $minsplit[1];
                        $evalue1 =~ /(.*)\..*/;
                        $evalue1 = $1;
                        $evalue1 =~ s/,//g;
			print "evalue: $evalue1 \n";
                        #$cmd = "biom summarize-table -i $otu -o library_stats.txt 2>&1 | tee -a $err";
                } else {
                        print "per_library_stats.py\n";
                        $evalue1 = `per_library_stats.py -i $otu | awk '/core_qiime_analyses.py \\(just a suggestion\\): /{print \$5}'`;
                        chomp($evalue1);
                        $evalue1=~ /(.*)\..*/;
                        $evalue1 = $1;
                        print "evalue: $evalue1 \n";
                        $cmd = "per_library_stats.py -i $otu > $out";
                        system($cmd);
                        if ($? == -1) { print "command failed: $!\n"; exit; }
                        print LOG $cmd."\n";
                }
        } else {
                print "$otu DNE; aborting biom summarize-table\n";
                print ERR "$otu DNE; aborting biom summarize-table\n";
        }
	print "beta_diversity_through_plots.py\n";
        if (-e "otu_table_".$taxon."_".$gg.".biom") {
		if (-e $gg."_97_pruned.tre") {
                	#$cmd = "beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue." -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue." -a -O $threads 2>&1 | tee -a $err";
			$cmd = "beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue." -t ".$gg."_97_pruned.tre -p ../qiime_params.txt -e ".$evalue." -a -O $threads 2>&1 | tee -a $err";
		} else {
			$cmd = "beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue. " -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue." -a -O $threads 2>&1 | tee -a $err";
		}
                if ($time eq "y") {
			if (-e $gg."_97_pruned.tre") {
                        	#$cmd = "(time beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue." -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue." -a -O $threads) > $time_pwd/time_beta_diversity.log 2>&1";
				$cmd = "(time beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue." -t ".$gg."_97_pruned.tre -p ../qiime_params.txt -e ".$evalue." -a -O $threads) > $time_pwd/time_beta_diversity.log 2>&1";
			} else {
				$cmd = "(time beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg.".biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue." -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue." -a -O $threads) > $time_pwd/time_beta_diversity.log 2>&1";
			}
                }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg.".biom DNE; aborting beta_diversity_through_plots.py\n";
                print ERR "otu_table_".$taxon."_".$gg.".biom DNE; aborting beta_diversity_through_plots.py\n";
        }
        if (-e "otu_table_".$taxon."_".$gg."_n1_noRoot.biom") {
		if (-e $gg."_97_pruned.tre") {
                	#$cmd = "beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue1."_n1_noRoot -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue1." -a -O $threads 2>&1 | tee -a $err";
			$cmd = "beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue1."_n1_noRoot -t ".$gg."_97_pruned.tre -p ../qiime_params.txt -e ".$evalue1." -a -O $threads 2>&1 | tee -a $err";
		} else {
			$cmd = "beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue1."_n1_noRoot -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue1." -a -O $threads 2>&1 | tree -a $err";
		}
                if ($time eq "y") {
			if (-e $gg."_97_pruned.tre") {
                        	#$cmd = "(time beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue1."_n1_noRoot -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue1." -a -O $threads) > $time_pwd/time_beta_diversity.log 2>&1";
				$cmd = "(time beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue1."_n1_noRoot -t ".$gg."_97_pruned.tre -p ../qiime_params.txt -e ".$evalue1." -a -O $threads) > $time_pwd/time_beta_diversity.log 2>&1";
			} else {
				$cmd = "(time beta_diversity_through_plots.py -i otu_table_".$taxon."_".$gg."_n1_noRoot.biom -m ../".$mapf." -o wf_beta_".$taxon."_".$gg."_e".$evalue1."n1_noRoot -t rep_phylo_".$gg.".tre -p ../qiime_params.txt -e ".$evalue1." -a -O $threads) > $time_pwd/time_beta_diversity.log 2>&1";
			}
                }
                system($cmd);
                print LOG $cmd."\n";
                if ($? == -1) { print "command failed: $!\n"; exit; }
        } else {
                print "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting beta_diversity_through_plots.py\n";
                print ERR "otu_table_".$taxon."_".$gg."_n1_noRoot.biom DNE; aborting beta_diversity_through_plots.py\n";
        }
        chdir('..');
        return;
}

#### SUBROUTINE: RAnalysis ####
sub RAnalysis() {
	my ($clust, $taxon, $gg) = @_;
	###############################################################
	#R Analysis
	###############################################################
	print "###############################################################\n";
	print "R Analysis\n";
	print "###############################################################\n";
	#Copy sl1p_analysis.Rmd to working directory so user has a copy
	chdir("$time_pwd/picked_otus_".$clust."_".$gg);
	$cmd = "cp $bin/sl1p_analysis.Rmd .";
	system($cmd);
	$cmd = "R -e \"params= list( otutable='otu_table_".$taxon."_".$gg."_n1_noRoot.txt', mapfile='../map_".$proj.".txt',trefile='".$gg."_97_pruned.tre',pwd='$time_pwd/picked_otus_".$clust."_".$gg."', L6file='wf_taxa_summary_".$taxon."_".$gg."_n1_noRoot/otu_table_".$taxon."_".$gg."_n1_noRoot_L6.txt', L2file='wf_taxa_summary_".$taxon."_".$gg."/otu_table_".$taxon."_".$gg."_L2.txt', libstats='library_stats_".$taxon."_".$gg.".txt'); rmarkdown::render('sl1p_analysis.Rmd')\"";
	system($cmd);
}

#### SUBROUTINE: cleanUp ####
sub cleanUp() {
	my ($clust, $taxon, $gg) = @_;
	###############################################################
	#Cleaning up
	###############################################################
	print LOG "###############################################################\n";
	print LOG "#Cleaning up\n";
	print LOG "###############################################################\n";
        print "###############################################################\n";
        print "#Cleaning up\n";
        print "###############################################################\n";
	chdir("$time_pwd/picked_otus_".$clust."_".$gg);	
	$cmd = "tar -zcvf align_$gg.tar.gz align_$gg/ && rm -R align_$gg/";
	system($cmd);
	$cmd = "tar -zcvf align_".$gg."_filtered.tar.gz align_".$gg."_filtered/ && rm -R align_".$gg."_filtered/";
	system($cmd);
	$cmd = "tar -zcvf ".$taxon."_".$gg."_assigned_taxonomy_tr.tar.gz ".$taxon."_".$gg."_assigned_taxonomy_tr/ && rm -R ".$taxon."_".$gg."_assigned_taxonomy_tr/";
	system($cmd);
	$cmd = "tar -zcvf seqs_".$gg."_otus.txt.tar.gz seqs_".$gg."_otus.txt && rm seqs_".$gg."_otus.txt";
	system($cmd);
	mkdir('rep_set_aux');
	$cmd = "mv rep_set.* rep_set_aux";
	system($cmd);
	$cmd = "tar -zcvf rep_set_aux.tar.gz rep_set_aux && rm -R rep_set_aux";
	system($cmd);
	mkdir('phylo_aux');
	$cmd = "mv rep_phylo_$gg.tre phylo_aux";
	system($cmd);
	$cmd = "tar -zcvf phylo_aux.tar.gz phylo_aux/ && rm -R phylo_aux/";
	system($cmd);
	$cmd = "tar -zcvf otus_aux.tar.gz otus_aux/ && rm -R otus_aux/";
	system($cmd);
	chdir('..');
	#Archival of mapping files and splits_dir are done separately right before exit call incase there are multiple picked_ folders.
}
