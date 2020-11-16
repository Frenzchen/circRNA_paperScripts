#!/usr/bin Rscript

# R 3.1
#
# Figure 2C
#
# Calculates the median phastCons scores for different exon types and circRNAs
# 'contained' - exon, which is contained in circRNA
# 'outside' - exon, which is outside of circRNA, but in parental gene
# 'non-parental' - exon, which is not in a parental gene
# exons overlapping with UTRs were excluded
#
# Note: phastCons scores were only available for mouse, rat and human
#
# Outputs PDF


# libraries
library(ggplot2) # v3.1.1


# run
##########################################################################################

species <- c("mouse", "rat", "human")
df <- data.frame(score=NULL, label=NULL, species=NULL)
d <- data.frame(contained=NULL, non.parental=NULL, outside=NULL, species=NULL)

for (s in species) {

	# read phastcons for all exons
	pcons <- read.table(paste("/archive/cig/kaessmann/fgruhl/", s, "/phastcons/phastcons.bed", sep=""), sep="\t", as.is=TRUE)
	pcons$gene_id <- do.call(rbind, strsplit(pcons$V1, "\\|"))[,1]
	pcons$type <- do.call(rbind, strsplit(pcons$V1, "\\|"))[,2]
	pcons <- pcons[,c(1,7,6,8)]
	colnames(pcons)[c(1,3)] <- c("exon_id", "score")

	pcons$parental <- ifelse(pcons$gene_id %in% pcons$gene_id[grep("c", pcons$type)], 1, 0)
	pcons$contained <- ifelse(pcons$exon_id %in% pcons$exon_id[grep("c", pcons$type)], 1, 0)
	pcons$label <- NA
	pcons$label[pcons$type == "c" & pcons$parental == 1] <- "contained"
	pcons$label[pcons$type == "o" & pcons$parental == 1] <- "outside"
	pcons$label[pcons$type == "o" & pcons$parental == 0] <- "non-parental"

	sb <- subset(pcons, !is.na(label))[,c(3,7)]
	sb$label <- factor(sb$label, levels=c("non-parental", "outside", "contained"))
	sb$species <- s
	df <- rbind(df, sb)
	
	k <- data.frame(t(sapply(by(pcons[,"score"], pcons[, "label"], function(x) median(x)), I)), species=s)
	
	d <- rbind(d, k)	
}

df$species <- factor(df$species, levels=species)
d <- d[,c(2:3,1,4)]
d$species <- factor(d$species, levels=species)

print(d)
print("non-parental vs contained")
print(d[,3]-d[,1])
print("non-parental vs outside")
print(d[,2]-d[,1])
print("% difference explained by parental+outside")
print(round(c(d[,2]-d[,1])/c(d[,3]-d[,1]),2))

q <- ggplot(df, aes(x=species, y=score, group=interaction(species, label), fill=factor(label))) + geom_boxplot(notch=TRUE, outlier.shape=NA, lwd=0.25) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key=element_rect(colour='white'), legend.key.size=unit(0.3, "cm")) + scale_fill_manual(values=c("#CCCCCC", "#fdd0fb", "#663399"), labels=c("non-parental", "outside", "contained"), name="") + xlab("species") + ylab("PhastCons scores") + ggtitle("C: PhastCons score by exon type")


# plot
##########################################################################################

pdf("phastcons.pdf", width=8.27/2, height=11.69/6, useDingbats=FALSE)
	q
dev.off()

