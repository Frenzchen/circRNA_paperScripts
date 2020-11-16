#!/usr/bin/env Rscripts

# R 3.1
#
# Figure 5A, Supplementary Figure 13A
#
# Compare mean milliDiv value for top-5 dimers in shared and species-specific circRNAs
#
# Arguments
# - species name (e.g. opossum, mouse, rat, ...)
#
# Outputs pdf


# libraries
library(ggplot2)

set.seed(22121987)

# species name to be provided from command line
args <- commandArgs(TRUE)
species <- args[1]

# top-5 dimers for each species
if (species == "opossum") {
	r <- c("SINE1_Mdo+SINE1a_Mdo", "SINE1_Mdo+SINE1_Mdo", "MIR3+MIR3", "MIR3_MarsA+SINE1_Mdo", "SINE1a_Mdo+SINE1a_Mdo")
} else if (species == "mouse") {
	r <- c("B1_Mus1+B1_Mus2", "B1_Mm+B1_Mus1", "B1_Mus1+B1_Mus1", "B2_Mm2+B2_Mm2", "B1_Mm+B1_Mus2")
} else if (species == "rat") {
	r <- c("ID_Rn1+ID_Rn1", "ID_Rn1+ID_Rn2", "ID_Rn2+ID_Rn2", "BC1_Rn+ID_Rn1", "BC1_Rn+ID_Rn2")
} else if (species == "rhesus") {
	r <- c("AluSx+AluSx1", "AluSx1+AluSz", "AluSx1+AluY", "AluSx+AluSz", "AluSx+AluY")
} else if (species == "human") {
	r <- c("AluSx+AluSx1", "AluSx1+AluSz", "AluSx1+AluY", "AluSx+AluSz", "AluSx+AluY")
} else {
	print "unknown species"
}


# methods
##########################################################################################

# Determine signficance level of p-value and return corresponding label
#
# Arguments
# x: p-value, based on which the label is defined
#
# Returns
# p-value label (character) for the final plot
significance.label <- function(x) {
	if (x <= 0.001) {
		lb <- "***"
	} else if (x <= 0.01) {
		lb <- "**"
	} else if (x <= 0.05) {
		lb <- "*"
	} else {
		lb <- "ns"
	}
	return(lb)
}


# read data
##########################################################################################

df <- read.table("usedDimers_milliDiv.txt", sep="\t", as.is=TRUE, header=TRUE)
genes <- read.table(paste("../GLMs_paper/geneTable_", species, ".txt", sep=""), sep="\t", as.is=TRUE, header=TRUE)
m <- merge(df, genes[,c(1,7)], by="ensembl_gene_id")
m$dimer <- paste(m[,3], m[,4], sep="+")

m <- subset(m, dimer %in% r)
m <- subset(m, age %in% c("therian", species))

m$mean.milliDiv <- (m$milliDiv1 + m$milliDiv2)/2
pvalues <- data.frame(dimer=NULL, p=NULL, label=NULL)

for (i in r) {
	a <- subset(m, dimer == i)
	a <- subset(a, !is.na(mean.milliDiv))
	
	if (nrow(a) > 0) {
		p <- t.test(a$mean.milliDiv[a$age == "therian"], a$mean.milliDiv[a$age == species], alternative="less")$p.value
		pvalues <- rbind(pvalues, data.frame(dimer=i, p=p, label=significance.label(p)))
	} else {
		pvalues <- rbind(pvalues, data.frame(dimer=i, p=NA, label=NA))
	}
}

ymax <- max(max(m$mean.milliDiv, na.rm=TRUE)) 
pvalues$ymax <- ymax* 1.1


# plot
##########################################################################################

title <- paste("TopDimers_milliDiv_", species, ".pdf", sep="")
pvalues$p <- ifelse(!is.na(pvalues$p), round(pvalues$p, 3), NA)

pdf(title, width=8.27/2, height=11.69/5, useDingbats=FALSE)
	ggplot(m, aes(x=factor(dimer, levels=r), y=mean.milliDiv, group=interaction(dimer, age), fill=age)) + geom_boxplot(notch=TRUE, lwd=0.25, outlier.size=0.5) + theme_minimal(base_size=8) + scale_fill_manual(values=c("red", "blue"), labels=c("species-specific", "shared"), name="Overlap") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm"), axis.text.x=element_text(angle=30)) + xlab("Top5 dimers") + ylab("milliDiv") + ylim(0, ymax*1.15) + ggtitle(paste("MilliDivs", species, sep=", ")) + geom_text(data=pvalues, aes(x=factor(dimer, levels=r), y=ymax, label=paste("p = ", p, sep=""), group=NULL, fill=NULL, fontface="italic"), size=2.4)
dev.off()
