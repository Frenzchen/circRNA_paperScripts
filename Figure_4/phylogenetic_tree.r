#!/usr/bin/env Rscript

# R 3.1
#
# Figure 4C
#
# Prodcues phylogenetic tree for mouse TEs based on TE reference sequence.
#
# Outputs PDF


# librariers 
library(ggplot2) # v3.1.1
library(ggdendro) # v0.1.20


# run
##########################################################################################

# hard-coded, files needs to be replaced with species-specific file 
df <- read.table("mouse_distmx.txt", sep="\t", as.is=TRUE, row.names=1)
dd <- dist(df, method="euclidian")
hc <- hclust(dd, method="ward.D2")
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")


# plot
##########################################################################################

# hard-coded, file name needs to be replaced with species-specific name 
pdf("repeat_phylogeny_mouse.pdf", width=8.27*0.4, height=11.69/4, useDingbats=FALSE)
	ggdendrogram(hc, rotate=TRUE, size=2) + theme_minimal(base_size=8) + ggtitle("mouse repeats phylogeny") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10)) + xlab("") + ylab("euclidian distance")
dev.off()
