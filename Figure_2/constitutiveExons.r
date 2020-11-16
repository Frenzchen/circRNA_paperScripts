#!/usr/bin/env Rscript

# R 3.1
#
# Figure 2D
#
# Calculates how often (%) circRNA is dominant circRNA if expressed in hotspot and multiple tissues


# libraries
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3
library(stringr) # v1.4.0


# functions
##########################################################################################

# Calculates the mean number of tissues, exon types occur in
# exon types: c = contained, o = outside, u = utr
#
# Arguments
# x: list (split by parental gene ID) of exons and number of tissues they are expressed in
#
# Returns
# data frame with exon types and mean number of tissue it is expressed in
tissue.ratio <- function(x) {
	ratios <- c()
	for (i in c("c", "o", "u")) {
		ratios <- c(ratios, mean(x$tissue[x$type == i]))
	}
	return(data.frame(type=c("c", "o", "u"), ratio=ratios))
}


# run
##########################################################################################

scratch <- "/archive/cig/kaessmann/fgruhl" # needs to be adapted to git folder
species <- c("opossum", "mouse", "rat", "rhesus", "human")
df <- data.frame(type=NULL, ratio=NULL, species=NULL)

for (s in species) {
	d <- paste(scratch, s, "exonCounts/", sep="/")
	setwd(d) # hardcoded
	
	exons <- read.table("exonCounts.txt", sep="\t", as.is=TRUE, header=TRUE)
	index <- read.table("mapped_exons.txt", sep="\t", as.is=TRUE)
	index$type <- do.call(rbind, str_split(index$V2, "\\|"))[,2]
	index$gene_id <- do.call(rbind, str_split(index$V2, "\\|"))[,1]
	colnames(index)[1] <- "exon_id"

	m <- merge(exons, index[,c(1,3:4)], by="exon_id")
	genes <- unique(subset(m, type == "c")$gene_id)
	m <- subset(m, gene_id %in% genes)
	
	# filter for "expressed" exons (FPKM >= 1) 
	m$tissue <- apply(m[,2:4], 1, function(x) length(which(x >= 1)))

	ls <- split(m, m$gene_id)

	df <- rbind(df, data.frame(do.call(rbind, lapply(ls, function(x) tissue.ratio(x))), species=s))
}

setwd(paste(scratch, "plots", sep="/")) #hardcoded


# plot
##########################################################################################

pdf("constitutiveExons.pdf", width=8.27/2, height=11.69/6, useDingbats=FALSE)
	ggplot(df, aes(x=species, y=ratio, fill=type, group=interaction(type, species))) + geom_boxplot(notch=TRUE, outlier.shape=NA, lwd=0.25) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key=element_rect(colour='white'), legend.key.size=unit(0.3, "cm")) + guides(fill = guide_legend(override.aes = list(colour = NULL)), legend.margin = unit(0, "cm")) + scale_fill_manual(values=c("#663399", "#fdd0fb", "#CCCCCC"), labels=c("contained", "outside", "utr"), name="") + ylab("number of tissues") + ggtitle("D: Tissue frequency of exon types")
dev.off()










# ruby
# get bed-file from coding exons 
#File.open("codingExons.bed", 'w') do |output|
#File.open("codingExons.gtf", 'r').readlines.each do |line|
#line = line.strip.split("\t")
#exon_id = line[-1].match(/(?<=transcript_id\s")(\w*.\d+)/)
#output.puts ["chr#{line[0]}", line[3..4], exon_id, 0, line[6]].join("\t")
#end
#end

# bash
#bedtools intersect -a codingExons.bed -b ../phastcons/exons_collapsed.bed -wa -wb -f 0.9 | cut -f4,10 > mapped_exons.txt

