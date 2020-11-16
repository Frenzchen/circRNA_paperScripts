#usr/bin/env Rscript

# R 3.1
#
# Figure 1C
#
# Plots
# - Figure 1C: circRNA count and loci depending on CPM
# -- table: df, plot: q1
#
# Outputs PDF


# libraries
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3


# run
##########################################################################################

species <- c("opossum", "mouse", "rat", "rhesus", "human")

df <- read.table("cpm_hotspot_development.txt", sep="\t", as.is=TRUE, header=TRUE)

df$total_percent <- df$hs_loci/df$total_loci
df$circ_percent <- df$circ_hs/df$total_circ
df$label <- round(df$circ_hs/df$hs_loci, 1)
df$rpm <- factor(df$rpm)
df$species <- factor(df$species, levels=species)

m <- melt(df[,c(1,6:9)], id.vars=c("species", "rpm", "label"))
m$variable <- factor(m$variable, levels=c("circ_percent", "total_percent"))


# plot
##########################################################################################

q1 <- ggplot(m, aes(x=rpm, y=value*100, group=variable, fill=variable)) + geom_bar(stat="identity", position="dodge", lwd=0.25) + theme_bw(base_size=8) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.title=element_blank(), legend.key.size=unit(0.3, "cm"), strip.text.x = element_text(size = 8, face="bold", hjust=0), strip.background = element_rect(colour="white", fill="white")) + scale_fill_manual(values=c("#663399", "#CCCCCC"), labels=c("circRNAs located in hotspot", "hotspot loci")) + facet_wrap(~species, ncol=2) + ylab("percentage [%]") + xlab("cpm threshold") + geom_text(data=m[16:30, ], aes(x=rpm, y=(value*100)+5, label=label), hjust=1.25, size=2.4) + ggtitle("C: CircRNA hotspot loci by CPM")

pdf("cpm_figure_1.pdf", useDingbats=FALSE)
	q1
dev.off()