# The short-read library 16S rRNA gene sequencing pipeline #

![sl1p_graphic.jpg](https://bitbucket.org/repo/x8Ardxn/images/2793309766-sl1p_graphic.jpg)

## (sl1p, pronounced "slip"). ##
sl1p was developed in the <a href="http://surettelab.ca">Surette laboratory</a>.

sl1p is a software pipeline designed to simply the processing of 16S rRNA (and other marker gene) sequencing data by automating quality control, the picking of operational taxonomic units, and taxonomic assignment. sl1p also performs some basic analyses of the dataset, including an R markdown file which can aid the user in furthering these analyses. sl1p is publically available for any use and is maintained.

For information on installation and requirements, see `INSTALL.txt`.

For a comprehensive manual of sl1p, see `README.txt`. For any additional inquires or to report any bugs or errors, please use the Issues tab.

Citation information: Whelan FJ & Surette MG. 'A comprehensive evaluation of the sl1p pipeline for 16S rRNA gene sequencing analysis.' (2017) Microbiome. https://doi.org/10.1186/s40168-017-0314-2

Please note that sl1p wraps third party software into its pipeline. If you cite sl1p, you should also cite the tools that it uses. These can vary according to the parameters used; for the default parameters, sl1p uses:

-cutadapt

-sickle

-pandaseq

-AbundantOTU+

-uchime

-Qiime

-R's phyloseq, ggplot2, & getopt packages
