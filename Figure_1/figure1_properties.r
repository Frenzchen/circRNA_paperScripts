#!/usr/bin/env Rscript

# R 3.1
#
# Figure 1B, D and E
#
# Plots:
# - Figure 1B: number of circRNAs in different species and tissues
# -- table: circ.frequency, plot: q.circ.freq
# - Figure 1D: hotspot frequency in different species and tissues
# -- table: df.tsi, plot: q.tsi
# - Figure 1E: contribution ("expression strength") of top-1 and top-2 circRNA to hotspot
# -- table: df.dom, plot: q.dom
# 
# Outputs PDF


# libraries
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3
library(grid) # v3.6.0
library(gridExtra) # v2.3


# functions
##########################################################################################


# Compare how often circRNA remains dominant between tissues
#
# Arguments
# x: data frame with hotspots and rank of each circRNA within each hotspot
#
# Returns
# data frame with circRNA and number of tissues in which it is dominant
dominant.switch <- function(x) {
	x <- subset(x, position == 1)
	ntissue <- nrow(x)
	ncirc <- length(unique(x$circ_id))
	
	return(data.frame(ntissue=ntissue, ncirc=ncirc))
}

# Caluclation of circRNA dominance per tissue and hotspot
# (no accumulation possible because CPMs were not normalized across tissues!!!)
#
# Arguments
# x: list of hotspots in species
#
# Returns
# list (one entry per hotspot) with dominant circRNAs and number of tissues
hs.dominance <- function(x) {
	lv <- data.frame(NULL)
	ce <- data.frame(NULL)
	ts <- data.frame(NULL)
	
	x[,6:8] <- apply(x[,6:8], 2, function(x) x/sum(x))
	
	if (!is.na(x[,6])) {
		lv <- x[order(x[,6], decreasing = TRUE),][,c(1,6)]
		lv$tissue <- "lv"
		lv$position <- 1:nrow(lv)
		colnames(lv)[2] <- "percentage"
	} else {
		
	}
	
	if (!is.na(x[,7])) {
		ce <- x[order(x[,7], decreasing = TRUE),][,c(1,7)]
		ce$tissue <- "ce"
		ce$position <- 1:nrow(ce)
		colnames(ce)[2] <- "percentage"
	}

	if (!is.na(x[,8])) {
		ts <- x[order(x[,8], decreasing = TRUE),][,c(1,8)]
		ts$tissue <- "ts"
		ts$position <- 1:nrow(ts)
		colnames(ts)[2] <- "percentage"
	}

	sub.ls <- list(lv, ce, ts)
	x <- do.call(rbind, Filter(function(x) is.data.frame(x), sub.ls))
	
	x.switch <- dominant.switch(x)
	return(list(x, x.switch))
}


# run
##########################################################################################

species <- c("opossum", "mouse", "rat", "rhesus", "human")
base <- c("md", "mm", "rn", "rm", "hs")

circ.frequency <- data.frame(species=NULL, lv=NULL, ce=NULL, ts=NULL)
tsi <- data.frame(species=NULL, tsi=NULL, type=NULL)
hs.dom <- data.frame(circ_id=NULL, percentage=NULL, tissue=NULL, position=NULL, species=NULL)
circ.sw <- data.frame(single=NULL, total=NULL, species=NULL)

# loop through species to collect data
for (s in species) {		
	
	# read data
	df1 <- read.table(paste("/Volumes/Data/ArchiveData/fgruhl/", s, "rpms/rpm_0.01/circRNA_rpms.txt", sep="/"), sep="\t", as.is=TRUE, header=TRUE) 
	df2 <- read.table(paste("/Volumes/Data/ArchiveData/fgruhl/", s, "properties/circRNA_rpms.txt", sep="/"), sep="\t", as.is=TRUE, header=TRUE) 
	df3 <- read.table(paste("/Volumes/Data/ArchiveData/fgruhl/", s, "rpms/rpm_0.1/circRNA_rpms_0.1.txt", sep="/"), sep="\t", as.is=TRUE, header=TRUE) 	

	
	# circRNAs per tissue and species
	freq <- table(unlist(strsplit(df2$tissue, "\\+")))
	circ.frequency <- rbind(circ.frequency, data.frame(species=s, ce=freq[1], lv=freq[2], ts=freq[3]))

	
	# in how many tissues is circRNA/hotspot active?	
	sb <- subset(df2, in_hotspot=="yes")
	ls <- split(sb, sb$hotspot_id)
	
	circ.tsi <- unlist(lapply(df2$tissue, function(x) length(unlist(strsplit(x, "\\+")))))
	circ.tsi <- length(which(circ.tsi > 1))/length(circ.tsi)

	hs.tsi <- data.frame(tsi=do.call(rbind, lapply(ls, function(x) length(unique(unlist(strsplit(x$tissue, "\\+")))))))
	hs.tsi <- length(which(hs.tsi$tsi > 1))/nrow(hs.tsi)
	
	tsi <- rbind(tsi, data.frame(species=s, tsi=c(circ.tsi, hs.tsi), type=c("circRNA", "hotspot")))


	# hotspot dominance
	results.ls <- lapply(ls, function(x) hs.dominance(x))
	dom <- data.frame(do.call(rbind, lapply(results.ls, function(x) x[[1]])), row.names=NULL)	
	dom$species <- s
	hs.dom <- rbind(hs.dom, dom)
	
	df.switch <- data.frame(do.call(rbind, lapply(results.ls, function(x) x[[2]])), species=s, row.names=NULL)
	circ.sw <- rbind(circ.sw, df.switch)
}


# plots
##########################################################################################

# circRNAs per tissue and species 
circ.frequency <- melt(circ.frequency)
circ.frequency$variable <- factor(circ.frequency$variable, levels=rev(c("lv", "ce", "ts")))
q.circ.freq <- ggplot(circ.frequency, aes(x=factor(species, levels=rev(levels(species))), y=value, group=variable, fill=variable)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25) + theme_bw(base_size=8) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm"), axis.title.y=element_blank()) + scale_fill_manual(values=rev(c("lightsalmon3", "lightseagreen", "lightskyblue")), breaks=c("lv", "ce", "ts"), labels=c("liver", "cerebellum", "testis")) + coord_flip() + scale_x_discrete(labels=rev(species)) + labs(x="species", y="circRNA count", title="B: Number of circRNAs", fill="")

guides(fill = guide_legend(override.aes = list(colour = NULL)))

# plot: tissue presence of hotspot
q.tsi <- ggplot(tsi, aes(x=species, y=tsi*100, type=type, fill=type)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25) + theme_bw(base_size=8) + theme(legend.title=element_blank(), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm")) + scale_x_discrete(labels=species) + scale_fill_manual(values=c("grey75", "grey25"), labels=c("circRNAs", "hotspots")) + labs(x="species", y="percentage [%]", title="D: % of circRNAs and hotspots in multiple tissues", fill="")

# dominance of circRNAs
hs.dom$species <- factor(hs.dom$species, levels=species)
q.dom <- ggplot(subset(hs.dom, position < 3), aes(x=as.character(position), y=percentage*100, group=interaction(position, species), fill=species)) + geom_boxplot(notch=TRUE, lwd=0.25) + theme_bw(base_size=8) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size=unit(0.3, "cm")) + scale_fill_manual(values=c("#1F39B9", "#d01c8b", "#f1b6da", "#b8e186", "#4dac26"), labels=species) + scale_x_discrete(labels=c("top-1", "top-2")) + labs(x="rank", y="percentage [%]", title="E: CircRNA expression strength in hotspots", fill="")


# plots
##########################################################################################

pdf("Figure_1.pdf", useDingbats=FALSE, width=8.27, height=11.69*2/5)
	grid.arrange(arrangeGrob(q.circ.freq, layout_matrix=rbind(c(NA, 1, NA)), widths=c(0.3, 0.3, 0.4), nrow=1), arrangeGrob(q.tsi, q.dom, widths=c(1/2, 1/2), nrow=1), nrow=2)
dev.off()




