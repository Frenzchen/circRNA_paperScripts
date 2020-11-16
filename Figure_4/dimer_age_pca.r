#!/usr/bin/env Rscript

# R 3.1
#
# Figure 3B, D, Supplementary Figure 12
#
# Creates PCA for repeat dimers (all TE families and most enriched family)
# Command argument: species to be analyzed
# 
# Outputs PDF


# libraries 
library(ggplot2) # v3.1.1
library(stringr) # v1.4.0

options(scipen = 999)

args <- commandArgs(TRUE)
species <- args[1]


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


# analysis
##########################################################################################

# read data
dimers <- read.table(paste("dimers_", species, ".txt", sep=""), sep="\t", as.is=TRUE)
df <- read.table(paste(species, "_distmx.txt", sep=""), sep="\t", as.is=TRUE, row.names=1)
names <- read.table("repeat_names.txt", sep="\t", as.is=TRUE)

# frequency analysis
nb <- table(dimers$V3)[1]
nf <- table(dimers$V3)[2]
dimers.fl <- subset(dimers, V3 == "flanking")
top5 <- dimers.fl[order(dimers.fl$V2, decreasing=TRUE),]$V1[1:5]

dm <- data.frame(dimer=NULL, freq=NULL, enr=NULL)

for (i in top5) {
	f <- dimers$V2[dimers$V1 == i & dimers$V3 == "flanking"]
	b <- dimers$V2[dimers$V1 == i & dimers$V3 == "background"]
	enr <- (f/nf) / (b/nb)
	dm <- rbind(dm, data.frame(dimer=i, f=f, enr=enr))
}


# !!! hard-coded, files needs to be replaced with species file !!! 
cofold <- read.table("RNAcofold_pairwise_rodents.txt", sep="\t", as.is=TRUE, header=TRUE)

genes <- rownames(df)
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

# sort dimer frequenvy into matrix
mx1 <- matrix(0, ncol=l, nrow=l)
for (i in 1:l) {
	for (j in i:l) {
		s <- sort(c(genes[i], genes[j]))
		id <- paste(s[1], s[2], sep="+")
		v1 <- dimers$V2[dimers$V1 == id & dimers$V3 == "flanking"]
		v2 <- dimers$V2[dimers$V1 == id & dimers$V3 == "background"]
		v <- v1/v2
		if (length(v) == 0) {v <- 0} else {v <- v}
		mx1[i, j] <- v
		mx1[j, i] <- v
	}
}

# 0.6744 -> frequency included
# 0.6746 -> enrichment included
# 0.7596 -> freq excluded

# normalize cofold and dimer frequency matrix (mx = cofold results, df = distance matrix)
mx <- apply(mx, 2, function(x) std01(x, min(mx), max(mx)))
df <- apply(df, 2, function(x) std01(x, min(df), max(df)))
mx1 <- apply(mx1, 2, function(x) std01(x, min(mx1), max(mx1)))

# PCA (don't scale, because it's already standardized)
xx <- df + mx + mx1
pc <- prcomp(xx, center=TRUE, scale.=FALSE)
d <- data.frame(pc$rotation[,1:2], samples=genes)
d2 <- merge(d, names, by.x="samples", by.y="V1")



# generate line graph
b <- subset(dimers, V3 == "background")
f <- subset(dimers, V3 == "flanking")
m <- merge(b[,1:2], f[,1:2], by="V1")
colnames(m)[2:3] <- c("b", "f")
m$rep1 <- do.call(rbind, str_split(m$V1, "\\+"))[,1]
m$rep2 <- do.call(rbind, str_split(m$V1, "\\+"))[,2]
if (species == "opossum") {
	m$rep1[m$rep1 == "SINE1_Mdo"] <- "SINE-1_MD"
	m$rep2[m$rep2 == "SINE1_Mdo"] <- "SINE-1_MD"
	m$rep1[m$rep1 == "SINE1a_Mdo"] <- "SINE-1a_MD"
	m$rep2[m$rep2 == "SINE1a_Mdo"] <- "SINE-1a_MD"
}

m <- subset(m, (rep1 %in% genes & rep2 %in% genes))
m <- m[order(m$f, decreasing = TRUE),]

lines = data.frame(x=NULL, y=NULL, repeatname=NULL, dimer=NULL, rank=NULL, selfloop=NULL)

for (i in 1:5) {
	sf <- ifelse(m$rep1[i] == m$rep2[i], 1, 0)
	if (sf == 1) {
		tmp1 <- data.frame(x=d2$PC1[d2$samples == m$rep1[i]], y=d2$PC2[d2$samples == m$rep1[i]], repeatname=m$rep1[i], dimer=m$V1[i], rank = i, selfloop = sf)
		lines <- rbind(lines, tmp1)
	} else {
		tmp1 <- data.frame(x=d2$PC1[d2$samples == m$rep1[i]], y=d2$PC2[d2$samples == m$rep1[i]], repeatname=m$rep1[i], dimer=m$V1[i], rank = i, selfloop = sf)
		tmp2 <- data.frame(x=d2$PC1[d2$samples == m$rep2[i]], y=d2$PC2[d2$samples == m$rep2[i]], repeatname=m$rep2[i], dimer=m$V1[i], rank=i, selfloop=sf)
		lines <- rbind(lines, tmp1, tmp2)
	}
}

lines$rank <- factor(lines$rank, level=unique(lines$rank))

freq <- data.frame(freq=t(do.call(rbind, list(by(m[,"f"], m[, "rep1"], sum)))))
freq$samples <- rownames(freq)
d2 <- merge(d2, freq, by="samples")


# plots
# Comment: this script doesn'r really run on it's own, each plot is generated separately
# Plots need to be commented out by hand
##########################################################################################

# mouse
if (species == "mouse") {
	pdf(paste("dimer_pca_", species, ".pdf", sep=""), width=8.27*0.6, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x, y=y+(abs(y)*0.3), colour=rev(rank)), shape=1, size=4, stroke=2) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=guide_legend(ncol = 2, title.position="top", title="Repeat class", size=2), colour=guide_legend(ncol = 2, title.position="top", title="Top5 rank", override.aes=list(size = 1)), size=guide_legend(ncol = 2, title.position="top", title="Frequency")) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177')) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x+0.025, y=y+0.05, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm")) + ggtitle(paste("B: ", species, " repeat dimers"))
	dev.off()

# opossum	 
} else if (species == "opossum") {

# zoom-in
	pdf(paste("dimer_pca_", species, ".pdf", sep=""), width=8.27/4, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x-(abs(y)*0.02), y=y+(abs(y)*0.15), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=FALSE, colour=FALSE, size=guide_legend(title.position="top", title="Frequency", size=2)) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#f768a1', '#c51b8a','#7a0177'), breaks="MIR", "SINE") + geom_text(data=unique(lines[,c(1:3)]), aes(x=x, y=y+0.04, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.position="bottom", legend.box="horizontal", legend.margin=margin(t=0, unit='cm'), legend.box.just = "left", legend.spacing.x=unit(0.1, "cm"), legend.key.width=unit(0.3, "cm")) + ggtitle(species) + xlim(0.15, 0.4) + ylim(-0.3, 0.35)
	dev.off()
	
	# full plot
	pdf(paste(species, "_repeatPhylogeny.pdf", sep=""), width=8.27*0.6, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x-(abs(y)*0.02), y=y+(abs(y)*0.15), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=guide_legend(ncol = 2, title.position="top", title="Repeat class", size=2), colour=guide_legend(ncol = 2, title.position="top", title="Top5 rank", override.aes=list(size = 1)), size=guide_legend(ncol = 2, title.position="top", title="Frequency")) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#f768a1', '#c51b8a','#7a0177'), breaks="MIR", "SINE") + geom_text(data=unique(lines[,c(1:3)]), aes(x=x-0.025, y=y+0.05, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm")) + ggtitle(paste("A: ", species, " repeat dimers", sep="")) 
	dev.off()

# rat
} else if (species == "rat") {	
	pdf(paste("dimer_pca_", species, ".pdf", sep=""), width=8.27/4, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x+(abs(x)*0.005), y=y+(abs(y)*0.004), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=FALSE, colour=FALSE, size=guide_legend(title.position="top", title="Frequency", size=2)) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177'), breaks=c("ID")) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x, y=y+0.04, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.position="bottom", legend.box="horizontal", legend.margin=margin(t=0, unit='cm'), legend.box.just = "left", legend.spacing.x=unit(0.1, "cm"), legend.key.width=unit(0.3, "cm")) + ggtitle(species) + xlim(-0.165, -0.155) + ylim(0.05, 0.35)
	dev.off()
	
	# full plot
	pdf(paste(species, "_repeatPhylogeny.pdf", sep=""), width=8.27*0.6, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x-(abs(y)*0.005), y=y+(abs(y)*0.004), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=guide_legend(ncol = 2, title.position="top", title="Repeat class", size=2), colour=guide_legend(ncol = 2, title.position="top", title="Top5 rank", override.aes=list(size = 1)), size=guide_legend(ncol = 2, title.position="top", title="Frequency")) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177')) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x+0.06, y=y, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm")) + ggtitle(paste("B': ", species, " repeat dimers", sep="")) 
	dev.off()
	
# rhesus
} else if (species == "rhesus") {	
	pdf(paste("dimer_pca_", species, ".pdf", sep=""), width=8.27/4, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x+(abs(x)*0.007), y=y+(abs(y)*0.004), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=FALSE, colour=FALSE, size=guide_legend(title.position="top", title="Frequency", size=2)) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#7a0177', '#c51b8a'), breaks=c("Alu")) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x+0.1, y=y, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.position="bottom", legend.box="horizontal", legend.margin=margin(t=0, unit='cm'), legend.spacing.x=unit(0.1, "cm"), legend.key.width=unit(0.3, "cm")) + ggtitle(species) + xlim(0, 0.5) + ylim(0.05, 0.35)
	dev.off()

	# full plot
	pdf(paste(species, "_repeatPhylogeny.pdf", sep=""), width=8.27*0.6, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x-(abs(y)*0.007), y=y+(abs(y)*0.004), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=guide_legend(ncol = 2, title.position="top", title="Repeat class", size=2), colour=guide_legend(ncol = 2, title.position="top", title="Top5 rank", override.aes=list(size = 1)), size=guide_legend(ncol = 2, title.position="top", title="Frequency")) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177')) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x+0.06, y=y, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm")) + ggtitle(paste("C': ", species, " repeat dimers", sep="")) 
	dev.off()

# human
} else if (species == "human") {	
	pdf(paste("dimer_pca_", species, ".pdf", sep=""), width=8.27/4, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x+(abs(x)*0.007), y=y+(abs(y)*0.004), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=FALSE, colour=FALSE, size=guide_legend(title.position="top", title="Frequency", size=2)) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#7a0177', '#c51b8a'), breaks=c("Alu")) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x, y=y+0.075, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.position="bottom", legend.box="horizontal", legend.margin=margin(t=0, unit='cm'), legend.spacing.x=unit(0.1, "cm"), legend.key.width=unit(0.3, "cm")) + ggtitle(species) + xlim(0.325, 0.4) + ylim(-0.5, 0.2)
	dev.off()
	
	pdf(paste(species, "_repeatPhylogeny.pdf", sep=""), width=8.27*0.6, height=11.69/4, useDingbats=FALSE)
		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x+(abs(x)*0.007), y=y+(abs(y)*0.004), colour=rev(rank)), shape=1, size=4, stroke=2) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=guide_legend(ncol = 2, title.position="top", title="Repeat class", size=2), colour=guide_legend(ncol = 2, title.position="top", title="Top5 rank", override.aes=list(size = 1)), size=guide_legend(ncol = 2, title.position="top", title="Frequency")) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#7a0177', '#c51b8a'), breaks=c("Alu")) + geom_text(data=unique(lines[,c(1:3)]), aes(x=x, y=y+0.075, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.5, "cm")) + ggtitle(paste("D': ", species, " repeat dimers", sep=""))
	dev.off()
}





#!/usr/bin/ebv Rscript

library(ggplot2)
library(reshape2)

species <- c("opossum", "mouse", "rat", "rhesus", "human")
dm <- data.frame(dimer=NULL, freq=NULL, enr=NULL, species=NULL)
df <- data.frame(species=NULL, top5=NULL, rest=NULL)

for (s in species) {
	dimers <- read.table(paste("dimers_", s, ".txt", sep=""), sep="\t", as.is=TRUE)

	# frequency analysis
	nb <- table(dimers$V3)[1]
	nf <- table(dimers$V3)[2]
	top5 <- subset(dimers, V3 == "flanking")
	n <- sum(top5$V2)
	top5 <- top5[order(top5$V2, decreasing=TRUE),][1:5,]

	df <- rbind(df, data.frame(species=s, top5=sum(top5$V2/n), rest=1-sum(top5$V2/n)))

	for (i in top5$V1) {
		f <- top5$V2[top5$V1 == i]
		b <- dimers$V2[dimers$V1 == i & dimers$V3 == "background"]
		enr <- (f/nf) / (b/nb)
		dm <- rbind(dm, data.frame(dimer=i, f=f, enr=enr, species=s))
	}
}

df <- melt(df)
df$species <- factor(df$species, levels=species)
df$variable <- factor(df$variable, levels=c("top5", "rest"))

ggplot(df, aes(x=species, y=value, fill=variable)) + geom_bar(stat="identity", position="stack", color="black", lwd=0.25) + theme_minimal(base_size=8) + scale_fill_manual(values=c("#663399", "white"), labels=c("top5", "rest")) + coord_flip()



		ggplot(d2) + theme_minimal(base_size=8) + geom_line(data=subset(lines, selfloop == 0), aes(x=x, y=y, colour=rev(rank)), lwd=1.5) + geom_point(data=subset(lines, selfloop == 1), aes(x=x-(abs(y)*0.02), y=y+(abs(y)*0.15), colour=rev(rank)), shape=1, size=4, stroke=2, show.legend=FALSE) + geom_point(aes(x=PC1, y=PC2, fill=V2, size=freq), shape=21) + guides(fill=FALSE, colour=FALSE, size=guide_legend(title.position="top", title="Frequency", size=2)) + scale_colour_manual(values=c('#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5'), labels=5:1) + xlab("PC1") + ylab("PC2") + scale_fill_manual(values=c('#f768a1', '#c51b8a','#7a0177'), breaks="MIR", "SINE") + geom_text(data=unique(lines[,c(1:3)]), aes(x=x, y=y+0.04, label=repeatname), fontface="bold", size=2.8) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.position="bottom", legend.box="horizontal", legend.margin=margin(t=0, unit='cm'), legend.box.just = "left", legend.spacing.x=unit(0.1, "cm"), legend.key.width=unit(0.3, "cm")) + ggtitle(species) + xlim(0.15, 0.4) + ylim(-0.3, 0.35)
