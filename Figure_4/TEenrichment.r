#!/usr/bin/env Rscript

# R 3.1
# 
# Figure 4A, Figure 4-Figure supplement 1
#
# Compares the enrichment of different transposable elements in background and flanking introns.
#
# Outputs PDF.


# libraries
library(stringr) # v1.4.0
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3
library(gridExtra) # v2.3


# run
##########################################################################################

scratch <- "~"
species <- c("opossum", "mouse", "rat", "rhesus", "human")
base <- c("md", "mm", "rn", "rm", "hs")

df <- data.frame(Var1=NULL, flanking=NULL, background=NULL, enr=NULL, species=NULL)
plot.ls <- list()

for (i in 1:5) {
	setwd(paste(scratch, species[i], "repeats", sep="/"))

	if (i == 4) {
		bn <- nrow(read.table("backgroundIntrons_lifted.bed", sep="\t", as.is=TRUE))
		fn <- nrow(read.table("flankingIntrons_lifted.bed", sep="\t", as.is=TRUE))
	} else {
		bn <- nrow(read.table("backgroundIntrons.bed", sep="\t", as.is=TRUE))
		fn <- nrow(read.table("flankingIntrons.bed", sep="\t", as.is=TRUE))
	}

	b <- read.table("repeat_backgroundIntron.intersection", sep="\t", as.is=TRUE)
	b$class <- do.call(rbind, str_split(b$V4, "\\:"))[,3]
	b$class <- do.call(rbind, str_split(b$class, "\\_"))[,1]
	b$name <- do.call(rbind, str_split(b$V4, "\\:"))[,1]

	f <- read.table("repeat_flankingIntron.intersection", sep="\t", as.is=TRUE)
	f$class <- do.call(rbind, str_split(f$V4, "\\:"))[,3]
	f$class <- do.call(rbind, str_split(f$class, "\\_"))[,1]
	f$name <- do.call(rbind, str_split(f$V4, "\\:"))[,1]

	# class enrichment
	bclass <- data.frame(table(b$class))
	bclass$Freq <- bclass$Freq/bn
	fclass <- data.frame(table(f$class))
	fclass$Freq <- fclass$Freq/fn
	mclass <- merge(fclass, bclass, by="Var1")
	colnames(mclass)[2:3] <- c("flanking", "background")
	mclass <- mclass[order(mclass$flanking, decreasing=TRUE),][1:10,]
	print(species[i])
	print(wilcox.test(mclass$flanking, mclass$background, paired=TRUE))
	mclass$enr <- mclass$flanking/mclass$background
	mclass <- mclass[order(mclass$enr, decreasing=TRUE),]
	mclass$species <- species[i]
	
	df <- rbind(df, mclass)
}

opossum <- subset(df, species == "opossum")
opossum$Var1 <- factor(opossum$Var1, levels=rev(unique(opossum$Var1)))
rodents <- subset(df, (species == "mouse" | species == "rat"))
rodents$Var1 <- factor(rodents$Var1, levels=rev(unique(rodents$Var1)))
primates <- subset(df, (species == "rhesus" | species == "human"))
primates$Var1 <- factor(primates$Var1, levels=rev(unique(primates$Var1)))

plot.ls[[1]] <- ggplot(opossum, aes(x=Var1, y=flanking)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25, aes(fill=enr)) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.width=unit(0.6, "cm"), legend.key.height=unit(0.3, "cm"),  legend.position="bottom", axis.title.y=element_blank(), axis.text.x=element_text(vjust=0.5, size=8)) + scale_fill_gradient2(low="#CCCCCC", mid="white", high="#663399", midpoint=1, name="enrichment") + ylab("# TEs") + coord_flip() + ggtitle("opossum")

plot.ls[[2]] <- ggplot(subset(rodents, species=="mouse"), aes(x=Var1, y=flanking)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25, aes(fill=enr)) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.width=unit(0.6, "cm"), legend.key.height=unit(0.3, "cm"), legend.position="bottom", axis.title.y=element_blank(), axis.text.x=element_text(vjust=0.5, size=8)) + scale_fill_gradient2(low="#CCCCCC", mid="white", high="#663399", midpoint=1, name="enrichment") + ylab("# TEs") + coord_flip() + ggtitle("mouse")

plot.ls[[3]] <- ggplot(subset(rodents, species == "rat"), aes(x=Var1, y=flanking)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25, aes(fill=enr)) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.width=unit(0.6, "cm"), legend.key.height=unit(0.3, "cm"), legend.position="bottom", axis.title.y=element_blank(), axis.text.x=element_text(vjust=0.5, size=8)) + scale_fill_gradient2(low="#CCCCCC", mid="white", high="#663399", midpoint=1, name="enrichment") + ylab("# TEs") + coord_flip() + ggtitle("rat")

plot.ls[[4]] <- ggplot(subset(primates, species=="rhesus"), aes(x=Var1, y=flanking)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25, aes(fill=enr)) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.width=unit(0.6, "cm"), legend.key.height=unit(0.3, "cm"), legend.position="bottom", axis.title.y=element_blank(), axis.text.x=element_text(vjust=0.5, size=8)) + scale_fill_gradient2(low="#CCCCCC", mid="white", high="#663399", midpoint=1, name="enrichment") + ylab("# TEs") + coord_flip() + ggtitle("rhesus")

plot.ls[[5]] <- ggplot(subset(primates, species=="human"), aes(x=Var1, y=flanking)) + geom_bar(stat="identity", position="dodge", colour="black", lwd=0.25, aes(fill=enr)) + theme_bw(base_size=8) + theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0, vjust=1, face="bold", size=10), legend.key.width=unit(0.6, "cm"), legend.key.height=unit(0.3, "cm"), legend.position="bottom", axis.title.y=element_blank(), axis.text.x=element_text(vjust=0.5, size=8)) + scale_fill_gradient2(low="#CCCCCC", mid="white", high="#663399", midpoint=1, name="enrichment") + ylab("# TEs") + coord_flip() + ggtitle("human")


# plots
##########################################################################################

grobz <- lapply(plot.ls, ggplotGrob)
lay <- rbind(c(1,NA, NA, NA), c(2,3,4,5))

output <- paste(scratch, "plots", "TEenrichment.pdf", sep="/")

pdf(output, height=11.69*1/2, width=8.27, useDingbats=FALSE)
	grid.arrange(grobz[[1]], grobz[[2]], grobz[[3]], grobz[[4]], grobz[[5]], layout_matrix=lay)
dev.off()	
	