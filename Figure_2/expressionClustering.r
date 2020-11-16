#!/usr/bin/env Rscript

# R 3.1
# 
# Figure 2B 
#
# Plots
# - expression clustering of circRNAs across species and tissues
# - expression clustering of parental genes across species and tissues
#
# Outputs PDF


# libraries
library(gplots) # v3.0.1.1
library(ggplot2) # v3.1.1


# run
##########################################################################################

# cluster
my_colors <- colorRampPalette(c("black", "blue", "white", "red"))(n = 399)
col_breaks <- c(seq(-0.33,0,length=100), seq(0.01,0.33,length=100), seq(0.34,0.66,length=100), seq(0.67,1,length=100)) 

# read expression for circRNAs
df <- read.table("expression2.txt", sep="\t", as.is=TRUE, header=TRUE)
mx <- cor(scale(df[df$age=="therian",][,2:16], center=TRUE, scale=TRUE), method="spearman")

# parental genes
df2 <- read.table("expression_parentalGenes.txt", sep="\t", as.is=TRUE, header=TRUE)
mx2 <- cor(df2[,7:length(df2)], method="spearman")


# plot
##########################################################################################

pdf("exprClustering.pdf", useDingbats = FALSE)
	# circRNAs
	heatmap.2(mx, margins=c(11, 11), col=my_colors, breaks=col_breaks, dendrogram="row", density.info="none", trace="none", main="clustering of \n therian circRNAs")

	# parental genes
	heatmap.2(mx2, margins=c(11, 11), col=my_colors, breaks=col_breaks, dendrogram="row", density.info="none", trace="none", main="clustering of \n therian parental genes")
dev.off()


