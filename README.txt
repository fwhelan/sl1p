# ======================================================== # 
# Fiona Whelan <whelanfj@mcmaster.ca>			   #
#							   #
# 11.12.13; 01.14.15; 05.12.16				   #
# The short-read library 16S rRNA gene seqing pipeline	   #
# sl1p.pl						   #
#							   #
# readme_sl1p.txt					   #
#							   #
# ======================================================== #


Contents
========
* What's new
* Installation
* Overview
    * Quick Start
    * User Input
    * Output
* Required Software
* Default processing
* Non-default processing

# ===================================================== #
# What's new						#
# ===================================================== #
-----------
Version 4.2
-----------
-fixes to the sl1p_analysis.Rmd headers to match Rmd formatting best practices
-adding -e functionality to allow the user to choose a taxonomic confidence cutoff
-----------
Version 4.1
-----------
-bug fixes in gg_pruned when >1 OTU maps to the same Greengenes Identifier
-rep_phylo.tre is now midpoint rooted
-implemented strict uclust closed-reference OTU picking (uclust-ref-strict)
-fastq input compressed with gzip now accepted as input
-existence of each fastq input in fofn is checked for before processing begins
-the minimum and maximum read lengths for input into split_lib has been adjusted to be a more stringent +-50bp from the expected amplicon (based on the primer locations)
-outputs an n1 no Root OTU table in relative abundance
-structure of output reorganized to help the user distinguish important, final files from intermediates
-fixed some bugs with the gg_pruned tree being opened in FigTree and R
-beta diversity metrics now run with pruned.tre as opposed to the rep_phylo.tre
-----------
Version 4.0
-----------
-fixed bug in -u thread option
-added a -l all option
-added functionality for additional primer sets that the user provides with a sequence information file (details below).
-sl1p calls multiple threads for each software it uses that has a multithreaded option.
-filter_tree now works with reference phylogenies other than gg2011 97.
-sl1p can now accept up to 64 fofn files
-bug fix (issue #84) that prevented uclust-ref with gg2013
-fofn files can now contain full paths to .fastqs in folders other than the one sl1p is run from
-Defaults updated: chimera checking to yes, and quality threshold to 30.
-----------
Version 3.9
-----------
-corrected bug in chimera checking
-fixed output to log_seq_info.txt
-added sickle quality trimming
-----------
Version 3.8
-----------
-Fixed #72: fixed a bug in the calls to cutadapt which should only trim primers and adapters, not trim based on quality.
-----------
Version 3.7
-----------
-there is now an option to assign taxonomy against all 3 databases (gg2011, gg2013, silva111) per run (see -d all)
-added the SILVA111 database as an option to assign taxonomy against
-assigning taxonomy used to have 3 options: blast, rdp, or trained rdp. Because we're using various databases, our
 rdp is always retrained; therefore the rdp option was removed.
-added UPARSE functionality.
-added chimera checking. sl1p now includes OTU table outputs that do not have chimeras in addition to the original
 tables. More details can be found below.
----------
Version 3.6
-----------
-changed the usage of sl1p to include flags as opposed to questions
-added the option of multiple OTU picking clustering thresholds
-added output file log_seq_info.txt which shows where sequences are culled along the pipeline
-----------
Version 3.5
-----------
-sl1p now creates 2 phylogenies; The first is created de novo based on alignments of the input sequences to the
 greengenes database with QIIME's make_phylogeny. This phylogeny is called rep_phylo.tre and is what sl1p has
 been producing since v1.0.  The second is created by pruning the greengenes reference
 tree down to only those nodes whom taxonomy has been assigned in the dataset; the nodes of this tree are subsequently
 renamed with the OTUs which were assigned to them. This phylogeny is called gg_97_pruned.tre.
-sl1p now chooses the most abundant read from a given OTU to use as its representative sequence
-change ConsensusLineage to Consensus Lineage to make pulling data into phyloseq easier.
-intermediate files are now removed to conserve space
-----------
Version 3.4
-----------
-taxa plots generated for .biom as well as _n1_noRoot.biom
-updated shebang
-updated QIIME's assign_taxonomy to be run with --rdp_max_memory 8000, 2x the original, to aid in running large datasets
-removed bug associated with having mixed barcodes
-corrected typo in removeRoot that prevented all Roots from being removed
-corrected the creation of multiple .err files
-edited the remove Root script in order to be able to handle the greengenes 2013 database release
-added alpha and beta outputs on the base OTU table
-----------
Version 3.3
-----------
-incorporated functionality for v34 region primers
-----------
Version 3.2
-----------
-updated log file to include usage information and details decided on at startup from the user
-log file now has a much easier to interpret name
-log file is much more verbose; every command is listed and comments explain what each block of code is doing
-a new err file is created with the same name as the log (but ending in err); this file captures all STDOUT and STDERR to better interpret errors
-----------
Version 3.1
-----------
-improved read me organization and version tracking
-updated shebang
-updated test of whether splits_dir should be overwritten (similar to behaviour for picked_otu_
-----------
Version 3.0
-----------
-incorporated handling of barcoded forward and reverse reads
-2 greengenes databases are now available: Feb 2011 and May 2013; user gets a choice upon startup
-sl1p now fails more gracefully

# ===================================================== #
# Installation						#
# ===================================================== #
To install sl1p on your machine, simply run the setup pipeline:
$sh setup_pipeline.sh
This will install the required software to your system. Note that administrator access may be necessary for setup_pipeline.sh to operate correctly. Please contact your system administrator as necessary.
In addition to the setup script, there are a few dependencies. Please see INSTALL.txt for more details.

# ===================================================== #
# Overview						#
# ===================================================== #
sl1p (pronounced "slip") is a tool designed to process 16S rRNA paired-end Illumina gene sequencing straight from the sequencer, through OTU picking, taxonomic assignment, and preliminary QIIME analyses. This tool is designed for the biologist with an emphasis on avoiding advanced knowledge of the unix commandline, without compromising detailed processing, including reproducibility and feasibility. Furthermore, sl1p acts as an academic exploratory tool that allows the user to compare differences in OTU picking, taxonomic assignment, and phylogeny creation while keeping all other elements in the processing consistent.

# ===================================================== #
# Quick start						#
# ===================================================== #
sl1p requires a File Of File Names (FOFN) file which lists the filenames of the paired-end fastq files, split by sample, of the samples that you wish to run through the pipeline. These filenames must be one file per line, and pairs must be in order (alphabetical order is recommended). Any overlapping barcodes (for e.g., from 2 separate Illumina sequencing runs) must be in different FOFN files.

Usage: sl1p.pl <# of runs> <FOFN for each run> <PROJNAME> [-r region] [-b barcode location] [-d taxonomy database] [-p OTU picking algorithm] [-l linkage] [-t taxonomic assignment method] [-f force overwrite]

    where the # of FOFNs reflects how many individual FOFN files are input to sl1p.pl
              FOFN for each run lists the name of each FOFN file, split by a blank space (no commas); the order of  these files does not matter
              PROJNAME is a project name of your choosing. This name will be used to identify the fasta and map files (see below). Note that special characters are not allowed in project names.
	      -r which region of the 16S rRNA gene sequenced [v3/v34] (default: v3)
	      -b location of the sequence barcode [fwd/rev/mixed] (default: fwd)
	      -d the greengenes database that you would like to use for taxonomic assignment [gg_2011/gg_2013] (default: gg2011)
	      -p OTU picking algorithm [abundantotu/uclust/cdhit/dnaclust/uclust-ref/uclust-ref-strict/mothur/blast/uparse/all] (default: abundantotu)
	      -c OTU picking clustering threshold [0-100] (default: 97)
	      -l the linkage method to use with mothur [nearest/average/furthest] (default: average)
	      -m whether the user would like to employ chimera checking [y/n] (default: y
	      -t the taxonomic assignment method [blast/rdp/rdp-training/all] (default: rdp-training)
	      -x do you want timed results for each step of the analysis [y/n] (default: n)
	      -u how many CPU threads you would like to envoke [0-9] (default: 1)
	      -f force overwrite .fna, map, picked_otu, and split_dir files from previous runs [y/n] (default: n)

# ===================================================== #
# User Input						#
# ===================================================== #
sl1p takes the raw FWD and REV .fastq files produced by Illumina sequencing.  These files must all be in the same directory - the directory sl1p is called from - and must end in the standard .fastq.  Additionally, a File Of File Names (FOFN) file must also be present in the same directory.  For example, if the user has 2 samples across 4 files, Sample1_R1_L001.fastq Sample1_R2_L001.fastq Sample2_R1_L001.fastq Sample2_R2_L001.fastq, these file names should be in the fofn file, each on their own line, such as:
Sample1_R1_L001.fastq
Sample1_R2_L001.fastq
Sample2_R1_L001.fastq
Sample2_R2_L001.fastq
Please ensure that there are no blank lines at the end of your FOFN file and that they appear in alphabetical order (it is essential that the R1 file directly preceeds the R2 file for a given sample).
The easiest way to create an FOFN is to pipe a list of all .fastq files in your directory into a text file. For example:
$ls *.fastq > FOFN.txt

Usage: sl1p.pl <# of FOFNs> <FOFN for each run> <PROJNAME>
where the # of FOFNs reflects how many individual FOFN files are input to sl1p.pl
FOFN for each run lists the name of each FOFN file, split by a blank space (no commas); the order of  these files does not matter
PROJNAME is a project name of your choosing. This name will be used to identify the fasta and map files (see below).
<<expand on this once I change it; omitting mention of questions now since they'll be gone soon>>

Sequence information file
-------------------------
If the user wishes to use primers other than those for the v3 and v34 region, they can be provided to sl1p in the form of a sequence information file. This file contains:
line 1: the forward primer sequence
line 2: the reverse primer sequence
line 3: the minimum sequence length included in the analysis
line 4: the maximum sequence length included in the analysis
line 5: the barcode location (fwd/rev)

Note: This file MUST have Unix line endings. This file cannot be generated in any sort of Document Editor such as Microsoft Word or Libre Office.

# ===================================================== #
# Output						#
# ===================================================== #
sl1p produces a series of output files, which are explained within:
log_sl1p_<date&time>.log -> logs how the user invoked sl1p as well as all commandline commands that sl1p performed based on user input; records all system commands as well as STDOUT.
log_sl1p_<date&time>.err -> logs any potential errors that may have occurred during the running of sl1p; captures STDERR.
<projname>.fna -> all trimmed, and aligned sequences from all samples as input for OTU picking
map_<projname>.txt -> map file generated by sl1p for input into QIIME
map_<projname>[_corrected.txt | .html | .log] -> output of QIIME's check_id_map.py
splits_dir/ -> output of QIIME's split_libraries.py
picked_otu_abundantou/
align/ -> output of QIIME's align_seqs.py
align_filtered/ -> output of QIIME's filter_alignment.py
library_stats.txt -> output of biom summarize-table -i otu_table_rdpgenus.biom
library_stats_n1_noRoot.txt -> output of biom summarize-table -i otu_table_rdpgenus_n1_noRoot.biom
otu_table_rdpgenus[.biom | .txt] -> output of QIIME's make_otu_table.py
otu_table_rdpgenus_nochim[.biom | .txt] -> otu_table_rdpgenus with any chimeras removed
otu_table_rdpgenus_noRoot[.biom | .txt] -> otu_table_rdpgenus with any OTUs assigned to Root or Unclassified removed using removeRoot.pl
otu_table_rdpgenus_n1[.biom | .txt] -> otu_table_rdpgenus with any singleton OTUs removed with QIIME's filter_otus_from_otu_table.py
otu_table_rdpgenus_n1_noRoot[.biom | .txt] -> otu_table_rdpgenus_n1 with any OTUs assigned to Root or Unclassified removed using removeRoot.pl
otu_table_rdpgenus_nochim_n1[.biom | .txt] -> otu_table_rdpgenus_nochim with any singletons removed with QIIME's filter_otus_from_otu_table.py
otu_table_rdpgenus_n2[.biom | .txt] -> otu_table_rdpgenus with any singleton and doubleton OTUs removed with QIIME's filter_otus_from_otu_table.py
otu_table_rdpgenus_nochim_n2[.biom | .txt] -> otu_table_rdpgenus_nochim with any singleton and doubleton OTUs removed with QIIME's filter_otus_from_otu_table.py
otu_table_rdpgenus_n2_noRoot[.biom | .txt] -> otu_table_rdpgenus_n2 with any OTUs assigned to Root or Unclassified removed using removeRoot.pl
rdpgenus_assigned_taxonomy_tr/ -> output of QIIME's assign_taxonomy.py
rep_phylo.tre -> output of QIIME's make_phylogeny.py
rep_set[.clust | .clustsize | .fna] -> output of OTU picking by AbundantOTU+
seqs_otus.txt -> rep_set.clust formatted for input into QIIME using abundantOTU_to_qiime.pl
wf_arare_rdpgenus/ -> output of QIIME's alpha_rarefaction.py with otu_table_rdpgenus.biom as input
wf_arare_rdpgenus_n1_noRoot/ -> output of QIIME's alpha_rarefaction.py with otu_table_rdpgenus_n1_noRoot.biom as input
wf_beta_rdpgenus_e[###]/ -> output of QIIME's beta_diversity_through_plots.py with otu_table_rdpgenus.biom as input; the e-value indicated in the folder name is that recommended by QIIME (see library_stats.txt for more information)
wf_beta_rdpgenus_e[###]_n1_noRoot/ -> output of QIIME's beta_diversity_through_plots.py with otu_table_rdpgenus_n1_noRoot.biom as input; the e-value indicated in the folder name is that recommended by QIIME (see library_stats_n1_noRoot.txt for more information)
wf_taxa_summary_rdpgenus/ -> output of QIIME's summarize_taxa_through_plots.py with otu_table_rdpgenus.biom as input
wf_taxa_summary_rdpgenus_n1_noRoot/ -> output of QIIME's summarize_taxa_through_plots.py with otu_table_rdpgenus_n1_noRoot.biom as input
qiime_params.txt -> a parameter file generated by sl1p for input into QIIME for preliminary analyses

# ===================================================== #
# Required software					#
# ===================================================== #
sl1p requires that QIIME (1.6.0 - 1.8.0) is installed on the system, with the optional RDP Classifier. See http://qiime.org/install/install.html for details.
All other necessary software to run sl1p with default settings can be installed via sh setup_pipeline.sh. These tools and the version's that sl1p has been tested with are listed below for completeness:
        Greengenes database (version 4Feb2011 and May2013)
        pandaseq 2.7
        AbundantOTU+
        cutadapt 1.4.2
<<the above isn't strickly true at the moment>>
Once this setup script has been run successfully, the user has the option to install additional software in order to give sl1p full functionality. sl1p was designed in order to allow the user the most choice of softwares possible; the full functionality thus relies on multiple different software and tools. These tools can be installed via <<add a separate install script for all of these>> and are listed below for completeness:
<<list all tools here>>

# ===================================================== #
# Default Processing					#
# ===================================================== #
Running sl1p with the defaults results in the below processing of the input samples. Appropriate changes to this process occur when non-default setting are chosen by the user.

#Cutadapt
The first step of the sl1p pipline is to trim the FWD (R1) and REV (R2) reads at the point where they run into the opposing primer.  This is an issue when the length of sequencing on the Illumina machine exceeds the length of the variable region used; for example, 200bp sequencing of the V3 region (~150bps) exceeds this region. Without trimming, the FWD and REV reads will not align properly with each other during the PANDAseq alignment (below).
input: R1 and R2 files within FOFNs.
output: clipped#/R1 and R2 files.

#PANDAseq
PANDAseq is used to align the FWD (R1) and REV (R2) reads for each sample individually. PANDAseq only aligns sequence pairs for which the FWD and REV primers are present (and match at a probability of 0.7); these primers are removed from the resulting aligned sequence. If the alignment contains unknown bases (Ns) it is culled from the output. For more information, please see the PANDAseq documentation: https://github.com/neufeld/pandaseq.
input: clipped#/R1 and R2 files.
output: clipped#/<sample name>.fasta; <sample name>.log.bz2

#Cutadapt
Cutadapt is used to remove any sequences which still contain either the FWD or REV primers. The presence of any primers may indicate mis-aligned sequences or PCR error. Additionally, any paired-end sequences which contain Illumina sequencing primers are culled. Illumina sequencing primers should not be included in any sequence; the presence of them here indicates a PCR error which are filtered out.
input: clipped#/<projname>#.fasta files
output: clipped#/<projname>#.fasta files

#Flip any Forward Barcoded Reads
In the rare causes of projects with a sub-set of samples barcoded on the FWD read and others barcoded on the REV read, the FWD-barcoded reads are reverse complemented here to mirror their REV-barcoded counterparts. This is necessary, not for taxonomic assignment, but in order to allow sequences from both of these types of barcodes to be clustered into the same OTU if sequence identity calls for it.
input: clipped#/<projname>#.fasta files
output: clipped#/<projname>#.fasta files

#Map file generated
A map file, necessary for QIIME analyses, is generated. This includes adding an additional nucleotide to the barcode if there are any overlapping barcodes used in the input set of samples.
output: map_<projname>.txt

#QIIME: check_id_map.py or validate_mapping_file.py (depending on the verion of QIIME)
Map file is verified by QIIME before use.
input: map_<projname>.txt
output: map_<projname>[.txt | .log | .html | _corrected.txt]

#QIIME: split_libraries.py
split_libraries.py does further quality filtering of the input sequences. Any sequences < 100 bps or > 250bps are removed; any sequences with any errors in the barcode are removed; any sequences with a homopolymer run of > 20bps is removed (this is the most homopolymer runs that is observed in the Greengenes 4Feb2011 database).
input: map_<projname>.txt, <projname>.fna
output: splits_dir/

#Cluster OTUs
sl1p includes multiple options for the clustering of OTUs, including greedy heuristic clustering and hierarchical methods; both taxonomy-independent (open-reference OTU picking) and -dependent methods (closed-reference OTU picking) are available.  It is up to the user and the data as to which they would like to envoke. We recommend AbundantOTU+ which, in our in silico replications, appears to be the most accurate. To use AbundantOTU+, choose the default settings upon startup, or choose abundantotu from the list of available clustering algorithms. After the OTUs have been clustered, a representative sequence is chosen from each cluster for further steps in the pipeline; how this is done is dependent on the OTU picking algorithm employed.
input: splits_dir/seqs.fna
output: picked_otus_<OTUpickingalgo>/rep_set.fna, picked_otus_<OTUpickingalgo>/seqs_otus.txt

#QIIME: align_seqs.py
QIIME's script aligns a representative of each OTU to a pre-calculated alignment of the 16S rRNA Greengenes database (options of 4Feb2011 or May2013 are chosen by the user at startup). This is necessary for the creation of a phylogenetic tree later in the pipeline (see make_phylogeny.py below).
input: picked_otus_<OTUpickingalgo>/rep_set.fna
output: picked_otus_<OTUpickingalgo>/align/

#Assign taxonomy
Similar to OTU clustering, sl1p includes multiple options for the assignment of taxonomy, the choice of which is up to the user. These methods have been tested within our laboratory, and we believe the RDP classifier trained on genus level data to be the most accurate.  To use RDP at the genus level, choose the default settings upon startup, or choose rdpgenus from the list of available taxonomic assignment methods.
input: picked_otus_<OTUpickingalgo>/rep_set.fna
output: picked_otus_<OTUpickingalgo>/rdpgenus_assigned_taxonomy_tr/

#QIIME: filter_alignment.py
Use QIIME's filter_alignment.py script to filter the alignment resulting from align_seqs.py for gap columns. Because the Greengenes database is vast compared to the input sequences, there may be uninformative gaps in the resulting alignment of the user's input sequences that are best removed.
input: picked_otus_<OTUpickingalgo>/align/
output: picked_otus_<OTUpickingalgo>/align_filtered/

#QIIME: make_phylogeny.py
Create a phylogeny from the filtered aslignment of all representative sequences of each OTU. This phylogeny is created de novo and is not built based on the Greengenes phylogeny.
input: picked_otus_<OTUpickingalgo>/align_filtered/
output: picked_otus_<OTUpickingalgo>/rep_phylo.tre

#QIIME: filter_tree.py
Create a filtered version of the gg2011 or gg2013 phylogeny based on the align_seq output for each OTU. Note: this cannot be done with the silva111 database.
input: rep_phylo.tre
output: gg[2011/2013]_clustThres_pruned.tre

#QIIME: make_otu_table.py
Tabulates the number of times an OTU is found in each sample and outputs into a .biom file along with the taxonomic assignment of each OTU.
input: picked_otus_<OTUpickingalgo>/seqs_otus.txt, picked_otus_<OTUpickingalgo>/rdpgenus_assigned_taxonomy_tr/rep_set_tax_assignments.txt
output: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus.biom

#QIIME: filter_otus_from_otu_table.py
sl1p uses QIIME's filter_otus_from_otu_table.py script to filter out any OTUs with <= 1 or <= 2 sequences in it (i.e. singleons and doubletons). It is believed that these OTUs represent unequal sampling of the rare taxa long tail and thus should not be used in comparisons between samples.
input: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus.biom
output: otu_table_taxon_n1.biom and _n2.biom

#removeRoot.pl
We have decided against the inclusion of any sequences which are not able to be assigned as Bacteria in our analyses.  This script removes any unassignable sequences and creates otu_table_taxon_noRoot.biom, otu_table_taxon_n1_noRoot.biom, .txt, and noRoot_filtered_otus.txt which lists the removed OTU #s.

#biom convert
sl1p uses biom convert multiple times in order to make .txt versions of the .biom OTU tables.
input: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus[.biom | _n1.biom | _n1_noRoot.biom | _n2.biom | _n2_noRoot.biom]
output: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus[.txt | _n1.txt | _n1_noRoot.txt | _n2.txt | _n2_noRoot.txt]

Even though sl1p is a processing pipeline, a few preliminary analyses are performed as described below. The authors of sl1p feel strongly that any analyses should be understood by the user and are unique to the dataset in question. That being said, its always satisfying to be able to take a quick look at one's data in the form of taxa summaries, rarefaction and beta diversity before performing more detailed and dataset-specific analyses.

#QIIME: summarize_taxa_through_plots.py
Summarizes and plots taxonomic information at varies taxonomic levels (from Phylum to Genus).
input: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus[.biom | _n1_noRoot.biom]
output: picked_otus_<OTUpickingalgo>/wf_taxa_summary_rdpgenus[_n1_noRoot]/

#QIIME: alpha_rarefaction.py
Rarefies the input OTU table; computes alpha diversity metrics; collates results; generates alpha rarefation plots.
input: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus[.biom | _n1_noRoot.biom]
output: picked_otus_<OTUpickingalgo>/wf_arare_rdpgenus[_n1_noRoot]/

#QIIME: per_library_stats.py or biom summarize-table (depending on the version of QIIME)
Summarizes the number of sequences per sample explained in the input OTU table. sl1p uses the suggested e-value (i.e. # of sequences per sample in the dataset) for calculations of beta diversity (see below).
input: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus[.biom | _n1_noRoot.biom]
output: picked_otus_<OTUpickingalgo>/library_stats[_n1_noRoot].txt

#QIIME: beta_diversity_through_plots.py
Performs beta diversity on the input OTU table and creates PCoA plots using the evalue calculated above. sl1p asks QIIME to use the Binary Jaccard, Bray Curtis, Unweighted Unifrac, and Weighted Unifrac metrics to produce 4 separate groups of plots.
input: picked_otus_<OTUpickingalgo>/otu_table_rdpgenus[.biom | _n1_noRoot.biom]
output: picked_otus_<OTUpickingalgo>/wf_beta_rdpgenus_e###[_n1_noRoot]/

QIIME has a great overview tutorial and great documentation for each individual script used in sl1p. Additionally, the log_sl1p_<date&time>.log file provides the user with the QIIME command line system calls that are used from within sl1p that may help the user gain a better understanding of the data processing steps taken by sl1p.

# ===================================================== #
# Non-default Processing                                #   
# ===================================================== #

#Taxonomic Reference Databases (-d)
sl1p provides options for which reference database the user would like to use to assign a taxonomic identifier to their OTUs.
	gg2011 -> 4 Feb 2011 release of greengenes
	gg2013 -> Aug 2013 release of greengenes

#OTU Picking (-p)
sl1p provides many options of OTU picking algorithms that the user can choose from.
	abundantotu -> AbundantOTU+, http://omics.informatics.indiana.edu/AbundantOTU/
	uclust -> as implemented in QIIME
	cdhit -> as implemented in QIIME
	dnaclust -> 
	clustom ->
	uclust-ref -> as implemented in QIIME
	mothur -> as implemented in QIIME; a secondary option, -l, allows the user to use different linkages
	blast -> as implemented in QIIME
	uparse

#Chimera checking (-m)
sl1p recommends checking for and removing chimeras from any 16S rRNA gene sequencing dataset. This option uses USEARCH to identify chimeras using a reference database (see -d) and using the de novo approach. Only those chimeras identified by both methods are culled.

#Taxonomic assignment

# ===================================================== #
# Addition of new databases/third party methods		#
# ===================================================== #

sl1p is publicly available, open-source software which is presented as-is. That being said, sl1p has an active author-base and will be maintained to suit user need and to adjust to new third party methods and databases as they become available. New versions, pull requests, issues, bugs, and suggested improvements can be viewed/requested at http://bitbucket.org/fwhelan/sl1p. Further, the user is welcome to edit or improve this software as per GNU copyright law. The sl1p code has been modularized and includes a separate subroutine per-third party option. A new algorithm or database option can be added by including a new subroutine and replacing the appropriate call to it in the Main method code.

