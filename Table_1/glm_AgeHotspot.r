#!usr/bin/env Rscript

# R 3.1
#
# Test whether hotspots are more likely to be shared than non-hotspot loci
#
# Outputs text-file (tab) with glm summary values for all species


set.seed(22121987)

library(reghelper) #v1.0.0


# run
##########################################################################################

species <- c("opossum", "mouse", "rat", "rhesus", "human")
glmvalues <- data.frame(estimate=NULL, std.error=NULL, p=NULL, species=NULL)

for (s in species) {
	
	df <- read.table(paste("~/Documents/scratch_folders/", s, "/GLMs_paper/geneTable_", s, "_glm.txt", sep=""), sep="\t", as.is=TRUE, header=TRUE)
	
	df <- subset(df, circle==1)
	print(nrow(df))
	
	# establish factor order
	if (s == "opossum") {
		f <- c(s, "therian")
	} else if (s == "mouse") {
		f <- c(s, "rodents", "eutherian", "therian")
	} else if (s == "rat") {
		f <- c(s, "rodents", "eutherian", "therian")
	} else if (s == "rhesus") {
		f <- c(s, "primates", "eutherian", "therian")
	} else if (s == "human") {
		f <- c(s, "primates", "eutherian", "therian")
	}
	
	df$age <- factor(df$age, levels=f)
	
	d <- as.data.frame(beta(glm(formula = hs ~ age, family = "binomial", data = df))$coefficients)[-1,c(1,2,4)]
	colnames(d) <- c("coefficient", "std.error", "p")
	d$lower.ci <- d$coefficient - 1.96 * d$std.error
	d$upper.ci <- d$coefficient + 1.96 * d$std.error
	d$species <- s
	
	d$predictor <- sub(".z", "", rownames(d))
	d <- d[,c(7, 1:6)]
	
	glmvalues <- rbind(glmvalues, d)
}

glmvalues[,c(2,3,5,6)] <- apply(glmvalues[,c(2,3,5,6)], 2, function(x) round(x,4))


# output
##########################################################################################

setwd("~/Documents/scratch_folders/GLMs_paper")
write.table(glmvalues, "glmValues_AgeHotspot.txt", sep="\t", quote=FALSE, append=FALSE, row.names=FALSE)

