#!/usr/bin/env Rscript

# R 3.1
#
# Figure 5A, Figure 5-Figure supplement 1
#
# Plot dimer enrichment in shared vs species-specific circRNAs together with dimer age
# Dimer age is based on mean lineage rank:
# - therian = rank 4
# - eutherian = rank 3
# - lineage (e.g. rodents) = rank 2
# - species (e.g. mouse) = rank 1
#
# Run from within folder "repeats" for each species.
#
# Arguments
# - species name (e.g. opossum, mouse, rat, ...)
#
# Outputs pdf and tab-delimited table with dimer frequency, enrichment and age


# libraries
library(stringr)
library(ggplot2)

#options("scipen"=0, "digits"=4)

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


# methods
##########################################################################################

# correction of dimers which were counted multiple times
# (see Figure 5-Figure supplement 2 for explanation)
#
# Arguments
# dim: dimer name (in script, loops through column with all dimer names)
# cor.df: data frame with correction values (overestimation_by_age.txt)
#
# Returns
# Data frame with corrected counts
correct_count <- function(dim, cor.df) {
	d1 <- strsplit(dim, "\\+")[[1]][1]
	d2 <- strsplit(dim, "\\+")[[1]][2]
	
	# therian factor for correction
	ft1 <- cor.df$therian.factor[cor.df$repeat_name == d1]
	ft2 <- cor.df$therian.factor[cor.df$repeat_name == d2]
	
	ft <- ifelse((!is.na(ft1) & !is.na(ft2)), round((ft1+ft2)/2), NA)
	if (!is.numeric(ft)) { ft <- NA }
	
	# species factor for correction
	fs1 <- cor.df$species.factor[cor.df$repeat_name == d1]
	fs2 <- cor.df$species.factor[cor.df$repeat_name == d2]
	
	fs <- ifelse((!is.na(fs1) & !is.na(fs2)), round((fs1+fs2)/2), NA)
	if (!is.numeric(fs)) { fs <- NA }
	
	return(data.frame(dimer=dim, species.factor=fs, therian.factor=ft))
}


# read adat and seperate into shared and species-specific repeats
##########################################################################################

# repeat overview
df <- read.table("dimers_by_sharedness_v2.txt", sep="\t", as.is=TRUE, header=TRUE)

r <- data.frame(n=do.call(cbind, list(by(df[, "count"], df[, "dimer"], mean))))
r$dimer <- rownames(r)
r <- r[order(r$n, decreasing=TRUE),]
rownames(r) <- 1:nrow(r)

top100 <- head(r, n=100)
base100 <- tail(r, n=100)

df <- subset(df, (dimer %in% top100$dimer | dimer %in% base100$dimer))
df <- subset(df, (orthology == s | orthology == "therian"))

# read correction table and correct counts in df
correction.df <- read.table("overestimation_by_age.txt", sep="\t", as.is=TRUE, header=TRUE)

dx <- data.frame(dimer=NULL, species.factor=NULL, therian.factor=NULL)

for (i in df$dimer) {
	dx <- rbind(dx, correct_count(i, correction.df))
}
dx <- unique(dx)

# get shared loci and associated dimers
# normalize dimer frequency by number of genes: function(x) sum(x)/nshared
shared <- subset(df, orthology == "therian")
nshared <- length(unique(shared$ensembl_gene_id))
shared <- data.frame(n=do.call(cbind, list(by(shared[, "count"], shared[, "dimer"], function(x) sum(x)/nshared))))
shared$dimer <- rownames(shared)
rownames(shared) <- 1:nrow(shared)

# correct for co-counting
for (i in 1:nrow(shared)) {
	shared$n[i] <- shared$n[i]/dx$therian.factor[dx$dimer == shared$dimer[i]]
}

# get species-specific loci and associated dimers
# normalize dimer frequency by number of genes: function(x) sum(x)/nspecies
species <- subset(df, orthology == s)
nspecies <- length(unique(species$ensembl_gene_id))
species <- data.frame(n=do.call(cbind, list(by(species[, "count"], species[, "dimer"], function(x) sum(x)/nspecies))))
species$dimer <- rownames(species)
rownames(species) <- 1:nrow(species)

# correct for co-counting
for (i in 1:nrow(species)) {
	species$n[i] <- species$n[i]/dx$species.factor[dx$dimer == species$dimer[i]]
}

# merge shared and species-specific 
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
p <- round(wilcox.test(mm$n.therian, mm$n.species, paired=TRUE)$p.value, 5)
xmax <- 0.9*max(mm$n.therian)
ymin <- 0.05*max(mm$n.species)
ymax <- 1.05*max(mm$n.species)

print(p)
print(tail(mm[order(mm$n.therian, mm$enrichment),][,c(3:6,9)], n=15))

# plot
##########################################################################################

plot.title <- paste("Dimer enrichment", s, sep=" ")
table.title <- paste("DimerAge_", s, ".txt", sep="")
title <- paste("DimerAge_", s, ".pdf", sep="")

write.table(mm, file=table.title, append=FALSE, quote=FALSE, row.names=FALSE)

pdf(title, width=8.27/2, height=11.69/5, useDingbats=FALSE)
	ggplot(mm, aes(x=n.therian, y=n.species, fill=enrichment, size=age)) + geom_abline(color="#CCCCCC", linetype="dashed") + geom_point(shape=21, alpha=0.5) + theme_minimal(base_size=8) + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) + scale_size(range = c(2, 6), limits=c(1,4), breaks=c(1, 1.5, 2, 2.5, 3, 3.5, 4)) + theme(plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.size = unit(0.5, "cm"), legend.margin=margin(0,0,0,0, unit="cm")) + ggtitle(plot.title) + guides(fill=guide_legend(ncol=2, title.position="top", title="log2 enrichment", size=2), size=guide_legend(ncol = 2, title.position="top", title="mean dimer age")) + xlab("Dimer count\n(shared circRNA loci)") + ylim(0, ymax) + ylab("Dimer count\n(species-specific circRNA loci)") + annotate(geom="text", label=paste("p =", p, paste=" "), size=2.4, x=xmax, y=ymin, fontface="italic")
dev.off()