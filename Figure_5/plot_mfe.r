#!/usr/bin/env Rscript

# R 3.1
#
# Figure 5C, Figure 5-Figure supplement 3
#
# Analyse mfe of least degraded dimers
#
# Arguments
# args[1]: species
# args[2]: dimers_rnacofold_rat.txt -> species name depending on species
# args[3]: usedDimer_repeats_v2.txt
# args[4]: reference file for dimer pairing (e.g. RNAcofold_pairwise_human.txt)
#
# Outputs PDF with final plot.


library(ggplot2) # v3.1.1
library(stringr) # v1.4.0

args <- commandArgs(TRUE)

species <- args[1]
rnacofold_file <- args[2] 
dimer_file <- args[3]
ref_file <- args[4]


# run
##########################################################################################

rnacofold <- read.table(rnacofold_file, sep="\t", as.is=TRUE, header=TRUE)
dimers <- read.table(dimer_file, sep="\t", as.is=TRUE, header=TRUE)

m <- merge(rnacofold, dimers[,1:2], by="rvc_id")
m$mean_milliDiv <- (m$milliDiv1+m$milliDiv2)/2

nspecies <- length(unique(m$gene_id[m$orthology == species]))
ntherian <- length(unique(m$gene_id[m$orthology == "therian"]))

# select least degraded dimer (based on mean milliDiv) for each gene
rows <- c()

for (gene in unique(m$gene_id)) {
	s <- subset(m, gene_id == gene)
	
	# can contain multiple dimers if they have the same mean milliDiv
	rn <- as.numeric(row.names(s[(s$mean_milliDiv == min(s$mean_milliDiv, na.rm=TRUE)),]))
	rows <- c(rows, rn)
}

rows <- rows[!is.na(rows)]
m2 <- m[rows,]

d <- data.frame(table(m2[,c(4,3)]))
ds <- d[as.numeric(rownames(d[d$orthology == species,])), c(1,3)]
dt <- d[as.numeric(rownames(d[d$orthology == "therian",])), c(1,3)]
d <- merge(ds, dt, by="dimer")
colnames(d)[2:3] <- c("species", "therian")
d$species <- d$species/nspecies
d$therian <- d$therian/ntherian
d$enrichment <- log2(d$therian/d$species)
d <- d[order(-d$therian, -d$enrichment),]

print(nspecies)
print(ntherian)
print(d)

# has been checked before manually, selection now added to make script run through
if (species == 'opossum') {
	dim <- droplevels(d$dimer[1])
} else if (species == 'mouse') {
	dim <- droplevels(d$dimer[c(1,2,4,5)])
} else if (species == 'rat') {
	dim <- droplevels(d$dimer[1:2])
} else if (species == 'rhesus') {
	dim <- droplevels(d$dimer[1])
} else if (species == 'human') {
	dim <- droplevels(d$dimer[1:2])
}

pvalues <- data.frame(dimer=NULL, p=NULL)

for (i in dim) {
	p <- t.test(m2$mfe[m2$orthology == species & m2$dimer == i], m2$mfe[m2$orthology == "therian" & m2$dimer == i])$p.value
	pvalues <- rbind(pvalues, data.frame(dimer=i, p=p))
}

# read mfe for reference sequence (need to split and order by name within dimer)
ref <- read.table(ref_file, sep="\t", as.is=TRUE, header=TRUE)
ref.tmp <- data.frame(do.call(rbind, lapply(str_split(ref$pair, "\\+"), function(x) x[order(x)])))
ref$pair <- paste(ref.tmp$X1, ref.tmp$X2, sep="+")

if (species %in% c('mouse', "rat", "human")) {
	ref <- subset(ref, pair %in% dim)
} else if (species == "opossum") {
	ref <- subset(ref, pair == "SINE-1_MD+SINE-1_MD")
	ref$pair[1] <- "SINE1_Mdo+SINE1_Mdo"
} else if (species == "rhesus") {
	ref <- subset(ref, pair == "AluMacYa3+AluY")
	ref$pair[1] <- "AluY+AluYRa1"
}

colnames(ref)[1] <- "dimer"


# output plot with dimers for which mfe is lower in therian than in species-specific circRNA loci
##########################################################################################

ymax <- max(m2$mfe[m2$dimer %in% dim], na.rm=TRUE)*0.85
ymin <- min(ref$mfe, na.rm=TRUE)*1.05
n.col <- length(dim)

plot.title <- paste("mfe_dimers_", species, ".pdf", sep="")

if (species %in% c('opossum', 'rhesus')) {
	w <- 0.3
} else if (species %in% c('rat', 'human')) {
	w <- 0.4
} else {
	w <- 0.5
}

pdf(plot.title, height=11.69/5, width=8.27*w, useDingbats=FALSE)
	ggplot(subset(m2, dimer %in% dim), aes(x=dimer, y=mfe, group=interaction(dimer, orthology), fill=orthology)) + geom_hline(data=ref, aes(yintercept=mfe), colour="#333333", lwd=1) + geom_point(aes(colour=orthology), position=position_jitterdodge(seed=1, dodge.width=0.9), size=0.1, show.legend=FALSE) + geom_violin(position=position_dodge(width=0.9), lwd=0.25, alpha=0.3, aes(fill=orthology)) + stat_summary(fun.y=median, geom="point", position=position_dodge(width=0.9), colour="#2b2d2f", size=0.75) + theme_minimal(base_size=8) + scale_fill_manual(values=c("red", "blue"), labels=c("species-specific", "shared"), name="CircRNA loci") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm"), strip.text=element_blank(), axis.text.x=element_text(hjust=1, angle=30)) + xlab("significant dimers") + ylab("mfe [kcal/mol]") + ggtitle(paste("Minimal Free Enereg (MEF)", species, sep=", ")) + geom_text(data=pvalues, aes(x=dimer, y=ymax, label=paste("p = ", round(p, 5), sep=""), group=NULL, fill=NULL, fontface="italic"), size=2.4) + facet_wrap(~dimer, scales="free_x", ncol=n.col)
dev.off()