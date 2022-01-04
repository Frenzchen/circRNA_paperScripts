#!/usr/bin/env Rscript

# R 3.1
#
# Figure 5C, Figure 5-Figure supplement 3
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
	r <- c("MAR1b_Mdo+MAR1b_Mdo", "SINE1_Mdo+SINE1a_Mdo", "MAR1a_Mdo+MAR1b_Mdo", "MAR1a_Mdo+MAR1a_Mdo", "SINE1_Mdo+SINE1_Mdo")
} else if (species == "mouse") {
	r <- c("B1_Mus1+B1_Mus2", "B2_Mm2+B2_Mm2", "B1_Mus1+B1_Mus1", "B2_Mm1t+B2_Mm2", "B1_Mm+B1_Mus1")
} else if (species == "rat") {
	r <- c("BC1_Rn+ID_Rn1", "ID_Rn2+ID_Rn2", "ID_Rn1+ID_Rn1", "BC1_Rn+ID_Rn2", "ID_Rn1+ID_Rn2")
} else if (species == "rhesus") {
	r <- c("AluSx+AluSx1", "AluSx+AluSx1", "AluSx+AluYRa1", "AluY+AluYRa1", "AluSx+AluSz")
} else if (species == "human") {
	r <- c("AluSx+AluSz", "AluSx1+AluSz", "AluSx+AluY", "AluSx1+AluY", "AluSx+AluSx1")
} else {
	print("unknown species")
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


# functions
##########################################################################################

df <- read.table("usedDimers_milliDiv_v2.txt", sep="\t", as.is=TRUE, header=TRUE)
df$dimer <- paste(df[,3], df[,4], sep="+")
df <- subset(df, dimer %in% r)

genes <- read.table(paste("../GLMs_paper/geneTable_", species, ".txt", sep=""), sep="\t", as.is=TRUE, header=TRUE)
df <- merge(df, genes[,c(1,7)], by="ensembl_gene_id")
df <- subset(df, age %in% c(species, "therian"))

df <- unique(data.frame(repeat_name=c(df[,3], df[,4]), pos=c(df[,5], df[,6]), milliDiv=c(df[,7], df[,8]), age=rep(df$age, 2)))

pvalues <- data.frame(repeat_name=NULL, p=NULL, label=NULL)

for (i in unique(df$repeat_name)) {
	a <- subset(df, repeat_name == i)
	a <- subset(a, !is.na(milliDiv))
	print(i)
	print(nrow(a))
	if (nrow(a) > 0) {
		p <- t.test(a$milliDiv[a$age == "therian"], a$milliDiv[a$age == species])$p.value
		pvalues <- rbind(pvalues, data.frame(repeat_name=i, p=p, label=significance.label(p)))
	} else {
		pvalues <- rbind(pvalues, data.frame(repeat_name=i, p=NA, label=NA))
	}
}

ymax <- max(df$milliDiv, na.rm=TRUE)
pvalues$ymax <- ymax* 1.1

title <- paste("TopRepeats_milliDiv_", species, ".pdf", sep="")
pvalues$p <- ifelse(!is.na(pvalues$p), round(pvalues$p, 3), NA)

pdf(title, width=8.27/2, height=11.69/5, useDingbats=FALSE)
	ggplot(df, aes(x=repeat_name, y=milliDiv, group=interaction(repeat_name, age), fill=age)) + geom_boxplot(notch=TRUE, lwd=0.25, outlier.size=0.5) + theme_minimal(base_size=8) + scale_fill_manual(values=c("red", "blue"), labels=c("species-specific", "shared"), name="CircRNA loci") + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm"), axis.text.x=element_text(angle=30)) + xlab("repeats of top-5 dimers") + ylab("milliDiv") + ylim(0, ymax*1.15) + ggtitle(paste("MilliDivs", species, sep=", ")) + geom_text(data=pvalues, aes(x=repeat_name, y=ymax, label=paste("p = ", p, sep=""), group=NULL, fill=NULL, fontface="italic"), size=2.4)
dev.off()
