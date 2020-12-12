#!/usr/bin/env Rscript

# R 3.1
#
# Figure 5B, Supplementary Figure 13B
#
# Plot dimer enrichment in shared vs species-specific circRNAs together with dimer age
# Dimer age is based on mean lineage rank:
# - therian = rank 4
# - eutherian = rank 3
# - lineage (e.g. rodents) = rank 2
# - species (e.g. mouse) = rank 1
#
# Arguments
# - species name (e.g. opossum, mouse, rat, ...)
#
# Outputs pdf and tab-delimited table with dimer frequency, enrichment and age


# libraries
library(stringr)
library(ggplot2)

set.seed(22121987)

# species name to be provided from command line
s <- commandArgs(TRUE)

if (s %in% c("mouse", "rat")) {
	lineage <- "rodents"
} else if (s %in% c("human", "rhesus")) {
	lineage <- "primates"
} else if (s == "opossum") {
	lineage <- NA
}


# read adat and seperate into shared and species-specific repeats
##########################################################################################

# repeat overview
df <- read.table("dimers_by_sharedness.txt", sep="\t", as.is=TRUE, header=TRUE)
r <- data.frame(n=do.call(cbind, list(by(df[, "count"], df[, "dimer"], sum))))
r$dimer <- rownames(r)
r <- r[order(r$n, decreasing=TRUE),]
rownames(r) <- 1:nrow(r)

top100 <- head(r, n=100)
base100 <- tail(r, n=100)

df <- subset(df, (dimer %in% top100$dimer | dimer %in% base100$dimer))

shared <- subset(df, orthology == "therian")
shared <- data.frame(n=do.call(cbind, list(by(shared[, "count"], shared[, "dimer"], mean))))
shared$dimer <- rownames(shared)
rownames(shared) <- 1:nrow(shared)

if (s %in% c("mouse", "rat")) {
	species <- subset(df, orthology == s)
} else if (s %in% c("human", "rhesus")) {
	species <- subset(df, orthology == lineage)
} else if (s == "opossum") {
	species <- subset(df, orthology == s)
}

#species <- subset(df, orthology == s)
species <- data.frame(n=do.call(cbind, list(by(species[, "count"], species[, "dimer"], mean))))
species$dimer <- rownames(species)
rownames(species) <- 1:nrow(species)

mm <- merge(shared, species, by="dimer")
colnames(mm)[2:3] <- c("n.therian", "n.species")
mm$enrichment <- log2(mm$n.therian/mm$n.species)
mm$rep1 <- do.call(rbind, str_split(mm$dimer, "\\+"))[,1]
mm$rep2 <- do.call(rbind, str_split(mm$dimer, "\\+"))[,2]

# define mean age of dimer
age <- read.table("repeatGroups.txt", sep="\t", as.is=TRUE, header=TRUE)
age <- age[age$age %in% c("eutherian", s, lineage, "therian"),]

if (s == "opossum") {
	age$age[age$age == "therian"] <- 4
	age$age[age$age == s] <- 1

} else {
	age$age[age$age == "therian"] <- 4
	age$age[age$age == "eutherian"] <- 3
	age$age[age$age == lineage] <- 2
	age$age[age$age == s] <- 1
}

mm <- merge(mm, age, by.x="rep1", by.y="te_name")
mm <- merge(mm, age, by.x="rep2", by.y="te_name")
mm$age <- (as.numeric(mm$age.x) + as.numeric(mm$age.y))/2

# compare whether there is difference of age in therian and species-specific repeats
p <- round(t.test(mm$n.therian, mm$n.species, paired=TRUE)$p.value, 3)
xmax <- 0.9*max(mm$n.therian)
ymin <- 0.05*max(mm$n.species)
ymax <- 1.05*max(mm$n.species)


# plot
##########################################################################################

plot.title <- paste("Dimer enrichment", s, sep=" ")
table.title <- paste("DimerAge_", s, ".txt", sep="")
title <- paste("DimerAge_", s, ".pdf", sep="")

write.table(mm, file=table.title, append=FALSE, quote=FALSE, row.names=FALSE)

pdf(title, width=8.27/2, height=11.69/5, useDingbats=FALSE)
	ggplot(mm, aes(x=n.therian, y=n.species, fill=enrichment, size=age)) + geom_abline(color="#CCCCCC", linetype="dashed") + geom_point(shape=21, alpha=0.5) + theme_minimal(base_size=8) + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) + scale_size(range = c(2, 6), limits=c(1,4), breaks=c(1, 1.5, 2, 2.5, 3, 3.5, 4)) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm"), legend.margin=margin(0,0,0,0, unit="cm")) + ggtitle(plot.title) + guides(fill=guide_legend(ncol=2, title.position="top", title="log2 enrichment", size=2), size=guide_legend(ncol = 2, title.position="top", title="mean dimer age")) + xlab("Dimer count\n(overlapping circRNA loci)") + ylim(0, ymax) + ylab("Dimer count\n(species-specific circRNA loci)") + annotate(geom="text", label=paste("p =", p, paste=" "), size=2.4, x=xmax, y=ymin, fontface="italic")
dev.off()