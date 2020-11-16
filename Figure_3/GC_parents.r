#!/usr/bin/env Rscript

# R 3.1
#
# Figure 3A
#
# Simple isochore analysis for coding egnes and parental genes
#
# Outputs PDF


# libraries
library(ggplot2) # v3.1.1
library(reshape2) # v1.4.3


# run
##########################################################################################

mm <- read.table("mouse/parentalGenes/parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
rn <- read.table("rat/parentalGenes/parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
rm <- read.table("rhesus/parentalGenes/parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
hs <- read.table("human/parentalGenes/parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)
md <- read.table("opossum/parentalGenes/parentSummary.txt", sep="\t", as.is=TRUE, header=TRUE)

z <- data.frame(rbind(prop.table(t(table(md$isochore[md$circle == 1])),1)*100, prop.table(t(table(md$isochore[md$circle == 0])),1)*100, prop.table(t(table(mm$isochore[mm$circle == 1])),1)*100, prop.table(t(table(mm$isochore[mm$circle == 0])),1)*100, prop.table(t(table(rn$isochore[rn$circle == 1])),1)*100, prop.table(t(table(rn$isochore[rn$circle == 0])),1)*100, prop.table(t(table(rm$isochore[rm$circle == 1])),1)*100, prop.table(t(table(rm$isochore[rm$circle == 0])),1)*100, prop.table(t(table(hs$isochore[hs$circle == 1])),1)*100, prop.table(t(table(hs$isochore[hs$circle == 0])),1)*100), species=c("opossum", "opossum", "mouse", "mouse", "rat", "rat", "rhesus", "rhesus", "human", "human"), circle=rep(c(1,0), 5))

z$species <- factor(z$species, levels=unique(z$species))
m <- melt(z, id.vars = c("species", "circle"))

labs <- data.frame(xval=c("L1", "L2", "L2", "L1", "L1"), yval=rep(80, 5), lab=c(z$L1[1]+z$L2[1], z$L2[3]+z$H1[3], z$L2[5]+z$H1[5], z$L1[7]+z$L2[7]+z$H1[7], z$L1[9]+z$L2[9]+z$H1[9]), species=c("opossum", "mouse", "rat", "rhesus", "human"), circle=unique(m$circle[m$circle == '1']))

df.lines <- data.frame(xval=c("L1", "L1", "L2", "L2", "L2", "L2", "H1", "H1", "L2", "L2", "H1", "H1", "L1", "L1", "H1", "H1", "L1", "L1", "H1", "H1"), yval=rep(c(73, 75, 75, 73), 5), species=c("opossum", "opossum", "opossum", "opossum", "mouse", "mouse", "mouse", "mouse", "rat", "rat", "rat", "rat", "rhesus", "rhesus", "rhesus", "rhesus", "human", "human", "human", "human"), circle=unique(m$circle[m$circle == '1']))


# plot
##########################################################################################

pdf("isochores.pdf", width=8.27, height=11.69/4)
	ggplot(m, aes(x=factor(variable, levels=c("L1", "L2", "H1", "H2", "H3")), y=value, group=circle, fill=factor(circle))) + geom_bar(stat="identity", pos="dodge", colour="black") + theme_bw(base_size=8) + guides(fill = guide_legend(keywidth=0.3, keyheight=0.3, override.aes = list(colour = NULL))) + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), plot.title=element_text(hjust=0, vjust=1, size=10, face="bold"), legend.key = element_rect(colour = "black"), legend.title=element_blank(), strip.text.x=element_text(face="bold", size=8, hjust=0), strip.background=element_rect(colour="white", fill="white")) + facet_wrap(~ species, scales="free_x", ncol=5) + scale_fill_manual(values=c("#CCCCCC", "#663399"), labels=c("non-parental", "parental")) + ylim(0,100) + xlab("isochores") + ylab("frequency (%) in isochores") + ggtitle("A: GC content of parental genes") + geom_text(data=labs, aes(x=xval, y=yval, label=paste(round(lab, 1), "%"), group=NULL), colour="#663399", size=3, hjust=0) + geom_line(data=df.lines, aes(x=factor(xval, levels=c("L1", "L2", "H1", "H2", "H3")), y=yval, group=circle), colour="#663399")
dev.off()