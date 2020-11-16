#!/usr/bin/env Rscript

# R 3.1
#
# Figure 2E, Supplementary Figure 3F and Supplementary Figure 6
# 
# Figure 2E
# Calculate correlations of GC amplitude and plot for different intron/exon boundaries
# Amplitude = ratio between the last 250 intronic bp and the first 50 exonic bp
# - 'n' - non-parental (but coding genes)
# - 'o' - parental, but outside of circRNA
# - 'c' - parental and contained in circRNA
#
# Supplementary Figure 3F
# Calculate the intron length for different intron types
# - 'np' - introns of non-parental genes (but in coding genes)
# - 'po' - introns parental genes but outside of circRNA
# - 'pf' - introns flanking circRNAs
# - 'pi' - introns inside of circRNAs
#
# Supplementary Figure 6
# Correlation = Amplitude correlation ~ GC content of intron
# - 'n' - non-parental (but coding genes)
# - 'o' - parental, but outside of circRNA
# - 'c' - parental and contained in circRNA
#
# Returns PDF


# libraries
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3
library(gridExtra) # v2.3
library(grid) # v3.6.0


# methods
##########################################################################################

# Transform character label into numeric label (used for stats.intron function)
#
# Arguments
# x: data frame with introns
#
# Returns
# numeric label
label <- function(x) {
	a <- c()
	if (x[12] == 0) {
		a <- 1
	} else if ((x[9] == "o" & x[12] == 1)) {
		a <- 2
	} else if ((x[9] == "c" & x[12] == 1)) {
		a <- 3
	} else if ((x[9] == "f" & x[12] == 1)) {
		a <- 4
	}
	return(a)
}

# Calculate median intron length for different intron types (np, po, pf, pi)
#
# Arguments
# x: list (by gene ID) of introns
#
# Returns
# data frame with gene id, "isochore of gene" and length of different intron types
stats.intron <- function(x) {
	x1 <- subset(x, label==1)
	x2 <- subset(x, label==2)
	x3 <- subset(x, label==3)
	x4 <- subset(x, label==4)
	
	x1.length <- ifelse(nrow(x1) == 0, NA, median(x1$length, na.rm=TRUE))
	x2.length <- ifelse(nrow(x2) == 0, NA, median(x2$length, na.rm=TRUE))
	x3.length <- ifelse(nrow(x3) == 0, NA, median(x3$length, na.rm=TRUE))
	x4.length <- ifelse(nrow(x4) == 0, NA, median(x4$length, na.rm=TRUE))

	d <- data.frame(gene=x$gene[1], isochore=x$isochore[1], np.length=x1.length, po.length=x2.length, pc.length=x3.length, pf.length=x4.length)
	return(d)
}

# Calculate median GC amplitude for different intron/exon boundaries types (n, o, c)
#
# Arguments
# x: list (by gene ID) of exon/intron GC amplitude (separated by acceptor and donor)
#
# Returns
# data frame amplitude, GC content of exon and intron and label (n, o, c)
stats.ampl <- function(x) {

	if (x$circle[1] == 1) {
	
		outside <- subset(x, circ_id == "0")
		ampl <- median(c(unique(outside[,c(14,16)])$amplitude_acc, unique(outside[,c(15,17)])$amplitude_don))
		igc <- median(c(unique(outside[,c(6,14)])$intron_acc, unique(outside[,c(9,14)])$intron_don))
		egc <- median(c(unique(outside[,c(7,14)])$exon_acc, unique(outside[,c(8,15)])$exon_don))
		d1 <- data.frame(amplitude=ampl, gc.intron=igc, exon.gc=egc, circle=1, label='o')	
		
		inside <- subset(x, circ_id != "0")
		ampl <- median(c(unique(inside[,c(14,16)])$amplitude_acc, unique(inside[,c(15,17)])$amplitude_don))
		igc <- median(c(unique(inside[,c(6,14)])$intron_acc, unique(inside[,c(9,14)])$intron_don))
		egc <- median(c(unique(inside[,c(7,14)])$exon_acc, unique(inside[,c(8,15)])$exon_don))
		d2 <- data.frame(amplitude=ampl, gc.intron=igc, exon.gc=egc, circle=1, label='c')	
		d <- rbind(d1, d2)
		
	} else {
	
		ampl <- median(c(unique(x[,c(14,16)])$amplitude_acc, unique(x[,c(15,17)])$amplitude_don))
		igc <- median(c(unique(x[,c(6,14)])$intron_acc, unique(x[,c(9,14)])$intron_don))
		egc <- median(c(unique(x[,c(7,14)])$exon_acc, unique(x[,c(8,15)])$exon_don))
		d <- data.frame(amplitude=ampl, gc.intron=igc, exon.gc=egc, circle=0, label='n')	
	}
	
	d$isochore <- x$isochore[1]
	d$ensembl_gene_id <- x$ensembl_gene_id[1]
	return(d)
}


# intron length
##########################################################################################

isochores <- c("L1", "L2", "H1", "H2", "H3")
species <- c("opossum", "mouse", "rat", "rhesus", "human")
base <- c("md", "mm", "rn", "rm", "hs")
path="/scratch/local/monthly/fgruhl"
df.length <- data.frame(gene=NULL, isochore=NULL, variable=NULL, value=NULL, circle=NULL, species=NULL)

for (s in species) {
	setwd(paste(path, s, "introns/introns_v2", sep="/"))
	o <- read.table("../../parentalGenes/parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
	o <- unique(o[,c(1,12)])
	d <- read.table("codingIntrons_gc.txt", sep="\t", as.is=TRUE)
	d$label <- do.call(rbind, strsplit(d$V5, "\\|"))[,3]
	d$gene <- do.call(rbind, strsplit(d$V5, "\\|"))[,1]
	d$length <- d$V3-d$V2

	d <- merge(d, o, by.x="gene", by.y="ensembl_gene_id")
	genes <- unique(d$gene[d$label != 'o'])
	d$circle <- ifelse(d$gene %in% genes, 1, 0)
	d$label <- apply(d, 1, function(x) label(x))
	
	ls <- split(d, d$gene)
	r <- do.call(rbind, lapply(ls, stats.intron))
	r <- melt(r, id.vars=c("gene", "isochore"))
	r <- subset(r, !is.na(r))
	r$circle <- ifelse(r$gene %in% genes, 1, 0)
	r$species <- s	
	df.length <- rbind(df.length, r)
}

pv.length <- data.frame(species=NULL, variable=NULL, enrichment=NULL, p=NULL)

# wilcox test for length difference
for (s in species) {
	sb <- df.length[df.length$species == s,]
	d2 <- subset(sb, (variable == "pf.length" & !is.na(value)))[c(1,4)]
	
	for (v in c("np.length", "po.length", "pc.length")) {		
		d1 <- subset(sb, (variable == v & !is.na(value)))[c(1,4)]
		
		if (v == "np.length") {
			p <- wilcox.test(d1[,2], d2[,2], alternative="less")$p.value
			enr <- median(d2[,2], na.rm=TRUE)/median(d1[,2], na.rm=TRUE)
		} else {
			colnames(d1)[2] <- v
			colnames(d2)[2] <- "pf.length"
			m <- merge(d1, d2, by="gene")
		 	p <- wilcox.test(m[,2], m[,3], paired=TRUE, alternative="less")$p.value
		 	enr <- median(m[,3], na.rm=TRUE)/median(m[,2], na.rm=TRUE)
		}
		pv.length <- rbind(pv.length, data.frame(species=s, variable=v, enrichment=enr, p=p))
	}
}

print(pv.length)

# length
df.length$species <- factor(df.length$species, levels=species)
qlength <- ggplot(df.length, aes(x=log10(value), group=variable, colour=variable)) + geom_density(lwd=0.6) + theme_bw(base_size=8) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm"), strip.text.x=element_text(face="bold", size=8, hjust=0), strip.background=element_rect(fill="white", colour="white")) + scale_colour_manual(values=c("grey90", "grey70", "grey40", "grey10"), labels=c("np", "po", "pi", "pf"), name="type") + xlab("intron length [ log10 ]") + ylab("density") + facet_wrap(~species, ncol=5) + ggtitle("A: Length of different intron types")


# GC amplitude and correlation between 
##########################################################################################

df <- data.frame(amplitude=NULL, gc.intron=NULL, exon.gc=NULL, circle=NULL, label=NULL, isochore=NULL, ensembl_gene_id=NULL, sepcies=NULL)

for (s in species) {
	wd <- paste(path, s, "parentalGenes", sep="/")
	setwd(paste(path, s, "parentalGenes", sep="/"))
	d <- read.table("exontable_gc.txt", sep="\t", as.is=TRUE, header=TRUE)
	d <- subset(d, distance == 249)
	o <- read.table("parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
	o <- unique(o[,c(1,12:13)])

	d <- merge(d, o, by="ensembl_gene_id")
	colnames(d)[c(11,13)] <- c("circ_id", "circle")
	d$acc <- do.call(rbind, strsplit(d$pos, ":"))[,2]
	d$don <- do.call(rbind, strsplit(d$pos, ":"))[,3]
	d$amplitude_acc <- d$exon_acc/d$intron_acc
	d$amplitude_don <- d$exon_don/d$intron_don
	ls <- split(d, d$ensembl_gene_id)

	r <- do.call(rbind, lapply(ls, function(x) stats.ampl(x)))
	r$species <- s
	df <- rbind(df, r)
}

pv.ampl <- data.frame(species=NULL, variable=NULL, p=NULL)

# wilcox test for amplitude difference
for (s in species) {
	sb <- df[df$species == s,]
	for (v in c("n", "o")) {
	
		if (v == "n") {
			p <- wilcox.test(sb$amplitude[sb$label == v], sb$amplitude[sb$label == "c"], alternative="less")$p.value
		} else {	
			d1 <- sb[sb$label == v,][c(1,7)]
			d2 <- sb[sb$label == "c",][c(1,7)]
			colnames(d1)[1] <- v
			colnames(d2)[1] <- "c"
			m <- merge(d1, d2, by="ensembl_gene_id")
			p <- wilcox.test(m[,2], m[,3], paired=TRUE, alternative="less")$p.value
		}	
		pv.ampl <- rbind(pv.ampl, data.frame(species=s, variable=v, p=p))
	}
}

print(pv.ampl)

df$species <- factor(df$species, levels=species)
q.ampl <- ggplot(df, aes(log2(amplitude), group=label, fill=label, alpha=label)) + geom_density(colour="black", lwd=0.25) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold"), legend.key.size=unit(0.3, "cm"), legend.key = element_rect(colour = "black"), strip.text.x = element_text(size = 8, face="bold", hjust=0), strip.background = element_rect(colour="white", fill="white")) + geom_vline(xintercep=0, color="red3") + scale_alpha_discrete(range=c(0.3, 0.6, 0.9), guide=FALSE) + guides(fill = guide_legend(override.aes = list(colour = NULL))) + scale_fill_manual(values=c("grey80", "grey60", "grey40"), labels=c("np", "po", "pi"), name="type") + xlab("amplitude [ log2 ]") + ylab("density") + ggtitle("B: Splice site amplitude") + xlim(-0.5, 1) + facet_wrap(~species, ncol=3)

# get idea of correlations
cors <- data.frame(species=NULL, isochore=NULL, label=NULL, intron=NULL, exon=NULL)

for (s in species) {
	sb <- subset(df, species == s)
	for (iso in isochores) {
		for (lb in c("n", "o", "c")) {
			cintron <- cor(sb$amplitude[sb$isochore == iso & sb$label == lb], sb$gc.intron[sb$isochore == iso & sb$label == lb], use="complete.obs", method="spearman")
			cexon <- cor(sb$amplitude[sb$isochore == iso & sb$label == lb], sb$exon.gc[sb$isochore == iso & sb$label == lb], use="complete.obs", method="spearman")
			d <- data.frame(species=s, isochore=iso, label=lb, intron=cintron, exon=cexon)
			cors <- rbind(cors, d)
		}
	}
}

m <- melt(cors)

q.cors <- ggplot(m, aes(x=label, y=value, group=interaction(variable, label), fill=variable)) + geom_boxplot(lwd=0.25, outlier.shape=NA, notch=TRUE, position="identity") + theme_bw(base_size=8) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.5, "cm"), legend.key=element_rect(colour="white")) + ylim(-0.8, 0.8) + scale_fill_manual(values=c("wheat2", "wheat4"), name="type") + scale_x_discrete(labels=c("np", "po", "pi")) + xlab("intron-/exon-type") + ylab("spearman's rho") + ggtitle("C: Amplitude")

print(mean(cors$intron[cors$label=="n"], method="spearman"))
print(mean(cors$intron[cors$label=="o"], method="spearman"))
print(mean(cors$intron[cors$label=="c"], method="spearman"))
print(mean(cors$exon[cors$label=="n"], method="spearman"))
print(mean(cors$exon[cors$label=="o"], method="spearman"))
print(mean(cors$exon[cors$label=="c"], method="spearman"))

#grobz <- lapply(list(q.ampl, q.cors), ggplotGrob)
#pl <- arrangeGrob(grobs=list(cbind(grobz[[1]], size = "last"), cbind(grobz[[2]], grobz[[3]], size = "last")), ncol=, nrow=2)


# write and plot
##########################################################################################

setwd(path)

write.table(cors, "amplitudeCorrelations.txt", sep="\t", append=FALSE, quote=FALSE)
lay <- rbind(c(1,1, 1), c(2,2, 3), c(2,2,3))

pdf("amplitude.pdf", useDingbats=FALSE, width=8.27, height=11.69/2)
	grid.arrange(qlength, q.ampl, q.cors, layout_matrix=lay)
dev.off()