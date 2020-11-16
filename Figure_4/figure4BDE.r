#!/usr/bin/env Rscript

# R 3.1
#
# Figure 4B, D and E
#
# 4B: plots the proportion of top5-dimers to all dimers
# 4C: PCA for phylogenetic distance matrix
# 4E: PCA for binding efficiency (deltaG)
#
# Outputs figure 4B / D and Figure 4E as seperate PDFs


library(ggplot2) # v3.1.1
library(gridExtra) # v2.3

# function
##########################################################################################

# Standardisation between 0 and 1
#
# Arguments
# x: value to be normalized
# minv: lowest value oberved
# maxv: highest value observed
#
# Returns
# stdx: standardised value
std01 <- function(x, minv, maxv) {
	stdx <- (x-minv)/(maxv-minv)
	return(stdx)
}


# analyses
##########################################################################################

# figure 4b analysis
species <- c("opossum", "mouse", "rat", "rhesus", "human")
df <- data.frame(percentage=NULL, category=NULL, species=NULL)

for (s in species) {
	d <- read.table(paste("dimers_", s, ".txt", sep=""), sep="\t", as.is=TRUE, header=FALSE)
	d <- subset(d, V3 == "flanking")
	d <- d[order(d$V2, decreasing=TRUE),]
	n <- round(sum(d$V2[1:5])/sum(d$V2)*100)
	df <- rbind(df, data.frame(percentage=c(n, 100-n), category=c("top5", "rest"), species=c(s, s)))
}

df$species <- factor(df$species, levels=rev(species))

q.top5 <-	ggplot(df, aes(x=species, y=percentage, group=species, fill=category)) + geom_bar(position="stack", stat="identity", colour="black", lwd=0.25) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm"), legend.title=element_blank(), legend.text=element_text(margin=margin(l=0.05, unit="cm")), legend.position="bottom", legend.spacing.x=unit(0.2, "cm"), axis.title.y=element_blank(), axis.text.x=element_text(vjust=0.5, size=8)) + scale_fill_manual(values=c("#FFFFFF", "#663399"), labels=c("remaining dimers", "top5-dimers"), name="") + xlab("species") + ylab("percentage [%]") + ggtitle("B: top-5 dimer contribution") + coord_flip()


# figure 4d analysis
dimers <- read.table("dimers_mouse.txt", sep="\t", as.is=TRUE)
distmx <- read.table("mouse_distmx.txt", sep="\t", as.is=TRUE, row.names=1)
names <- read.table("repeat_names.txt", sep="\t", as.is=TRUE)
genes <- rownames(distmx)

distmx <- apply(distmx, 2, function(x) std01(x, min(distmx), max(distmx)))
pc <- prcomp(distmx, center=TRUE, scale.=FALSE)
d <- data.frame(pc$rotation[,1:2], samples=genes)
d2 <- merge(d, names, by.x="samples", by.y="V1")

q.pca.dist <- ggplot(d2) + theme_minimal(base_size=8) + geom_point(aes(x=PC1, y=PC2, fill=V2), shape=21, size=3) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177', '#000000'), name="TE family") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm")) + ggtitle("D: TE distance matrix")


# figure 4e analysis
# !!! hard-coded, files needs to be replaced with species file !!! 
cofold <- read.table("RNAcofold_pairwise_rodents.txt", sep="\t", as.is=TRUE, header=TRUE)
l <- length(genes)
mx <- matrix(0, ncol=l, nrow=l)

# sort cofold results into matrix
for (i in 1:l) {
	for (j in i:l) {
		s <- sort(c(genes[i], genes[j]))
		id <- paste(s[1], s[2], sep="+")
		v <- cofold$mfe[cofold$pair == id]
		if (length(v) == 0) {v <- 0} else {v <- -1*(mean(v))}
		mx[i, j] <- v
		mx[j, i] <- v
	}
}

mx <- apply(mx, 2, function(x) std01(x, min(mx), max(mx)))
pc.mx <- prcomp(mx, center=TRUE, scale.=FALSE)
d.mx <- data.frame(pc.mx$rotation[,1:2], samples=genes)
d.mx <- merge(d.mx, names, by.x="samples", by.y="V1")

q.pca.deltag <- ggplot(d.mx) + theme_minimal(base_size=8) + geom_point(aes(x=PC1, y=PC2, fill=V2), shape=21, size=3) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177', '#000000'), name="TE family") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm")) + ggtitle("E: DeltaG dimers")


# output plots
##########################################################################################

pdf("figure4bd.pdf", width=8.27*0.6, height=11.69/4, useDingbats=FALSE)
	grid.arrange(q.top5, q.pca.dist, ncol=2, widths=c(0.4, 0.6))		
dev.off()

pdf("figure4e.pdf", width=8.27*0.4, height=11.69/4, useDingbats=FALSE)
	q.pca.deltag
dev.off()
