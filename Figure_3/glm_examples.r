#!/usr/bin/env Rscript

# R 3.1
#
# Figure 3C-D
#
# Plot 500 random p-values from GLM
#
# Outputs PDF


# librariers
library(ggplot2) # v3.1.1
library(gridExtra) # v2.3

set.seed(22121987)


# run
##########################################################################################

df <- read.table("../../GLM/geneTable_row.txt", sep="\t", as.is=TRUE, header=TRUE)
pvalues <- read.table("predictions_parentalGLM.txt", sep="\t", as.is=TRUE, header=TRUE)

df <- df[pvalues$rows,]
df$p <- pvalues$pred

p1 <- ggplot(df[sample(1:nrow(df), 500), ], aes(x=exon_count, y=percentage_gc_content, fill=p, size=p)) + geom_point(shape=21, alpha=0.5, stroke=0.25) + theme_minimal(base_size=10) + scale_size_continuous(range=c(0.5, 7), name="Probability") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm"), legend.title=element_text(size=8)) + scale_fill_gradient(low="#FFFFFF", high="#663399", name="Probability", guide="legend") + xlab("#exons") + ylab("GC content [%]") + ggtitle("C: GC content vs exon count")

p2 <- ggplot(df[sample(1:nrow(df), 500), ], aes(x=log10(as.rvc), y=phastcons, fill=p, size=p)) + geom_point(shape=21, alpha=0.5, stroke=0.25) + theme_minimal(base_size=10) + scale_size_continuous(range=c(0.5, 7), name="Probability") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm"), legend.title=element_text(size=8)) + scale_fill_gradient(low="#FFFFFF", high="#663399", name="Probability", guide="legend") + xlab("reverse-complements [log10(#)]") + ylab("phastCons score") + ggtitle("D: RVCs vs phastCons score")


# plot
##########################################################################################

pdf("GLM_visuals.pdf", width=8.27, height=11.69/5, useDingbats=FALSE)
	grid.arrange(p1, p2, ncol=2)
dev.off()