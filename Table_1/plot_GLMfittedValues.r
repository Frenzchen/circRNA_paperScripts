#!/usr/bin/env Rscript

# R 3.1
#
# Plots distribution of predicted values for glms separated by parental / non-parental


library(ggplot2) # v3.1.1


# run and output
##########################################################################################


species <- c("opossum", "mouse", "rat", "rhesus", "human")
df <- data.frame(trueVal=NULL, pred=NULL, rows=NULL, species=NULL)

for (s in species) {
	setwd(paste("~/Documents/scratch_folders", s, "GLMs_paper", sep="/"))
	d <- read.table("predictions_parentalGLM.txt", header=TRUE, sep="\t", as.is=TRUE)
	df <- rbind(df, d)
}

df$species <- factor(df$species, levels=species)


pdf("Pvalues_parentalGenes.pdf", width=8.27, height=11.69/5, useDingbats=TRUE)
	ggplot(df, aes(pred, colour=factor(trueVal))) + geom_density(lwd=0.5) + theme_minimal(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold"), legend.key.size=unit(0.3, "cm"), legend.key = element_rect(colour = "black"), strip.text.x = element_text(size = 8, face="bold", hjust=0), strip.background = element_rect(colour="white", fill="white")) + xlim(0,1) + xlab("Prediction values") + ylab("Density") + scale_colour_manual(values=c("#CCCCCC", "#663399"), name="", labels=c("non-parental", "parental")) + ggtitle("Prediction values for parental gene GLM") + facet_wrap(~ species, scales="free_y", ncol=5)
dev.off()
