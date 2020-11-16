#!/usr/bin/env Rscript

# R 3.1
#
# Figure 3B
#
# Plot repeat cpmplementarity depending on GC content in coding genes
#
# Outputs PDF


# libraries
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3


# run
##########################################################################################

# get repeat density plot
species <- c("opossum", "mouse", "rat", "rhesus", "human")
base <- c("md", "mm", "rn", "rm", "hs")

# get megablast plot
mb <- read.table("megablastSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
m <- melt(mb)
m$orientation <- ifelse(m$orientation == "as", "antisense", "sense")
m$species <- factor(m$species, levels=species)
m$isochore <- factor(m$isochore, levels=c("L1", "L2", "H1", "H2", "H3"))

q.mb <- ggplot(m[grep("coding", m$variable), ], aes(x=isochore, y=value, group=interaction(species, variable), colour=species, linetype=variable)) + geom_line(lwd=0.5) + theme_bw(base_size=8) + facet_wrap(~orientation) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.5, "cm"), legend.title=element_blank(), strip.text.x=element_text(face="bold", size=8, hjust=0), strip.background = element_rect(colour="white", fill="white"), legend.key=element_rect(colour=NA), legend.margin = margin(t=0, unit="cm")) + scale_colour_manual(values=c("#1F39B9", "#d01c8b", "#f1b6da", "#b8e186", "#4dac26")) + scale_linetype(labels=c("non-parental", "parental")) + xlab("isochores") + ylab("%") + ggtitle("B: Complementarity in coding genes")


# plot 
#########################################################################################€

pdf("repeatComplementarity.pdf", useDingbats=FALSE, width=8.27, height=11.69*1/5)
	q.mb
dev.off()


