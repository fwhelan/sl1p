#!/usr/bin/Rscript

#Author: Fiona J. Whelan
#Date: December 7th 2016

#Load libraries
library(getopt)
library(phyloseq)

#Get user input
spec = matrix(c('otutable', 'o', 1, "character", 'mapfile', 'm', 1, "character", 'relfile', 'r', 1, "character"), byrow=TRUE, ncol=4)
#spec = matrix(c('otutable', 'o', 1, "character", 'relfile' , 'r', 1, "character"), byroow=TRUE, ncol=4)
opt = getopt(spec)
otuFile <- opt$otutable
mapFile <- opt$mapfile
relFile <- opt$relfile

#Load phyloseq object
qd <- import_qiime(otuFile, mapFile)
#qd <- import_biom(otuFile, parseFunction=parse_taxonomy_greengenes)

#Do proportional normalization
#Code a modified form of that supplied in the supplemental data of Waste Not, Want No
normf = function(x) {
	x/sum(x)
}
qd.norm = transform_sample_counts(qd, normf)

#Output qd.norm into a relative abundance table
qd.df <- otu_table(qd.norm)
write.table(qd.df, relFile, sep="\t")
