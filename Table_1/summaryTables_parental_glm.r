#!/sur/bin/env Rscript

# R 3.1
#
# Summarizes individual GLM tables to one single table


# functions
##########################################################################################

# Determine signficance level of p-value and return corresponding label
#
# Arguments
# x: p-value (numeric)
#
# Returns
# lb: label
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


# run
##########################################################################################

# summarise the bestglm and prediction data set
species <- c("opossum", "mouse", "rat", "rhesus", "human")
df.bestglm <- data.frame(predictor=NULL, coefficient=NULL, std.error=NULL, p.value=NULL,lower.ci=NULL, upper.ci=NULL, species=NULL)
df.pred <- data.frame(predictor=NULL, beta.coefficient=NULL, std.error=NULL, p.value=NULL,lower.ci=NULL, upper.ci=NULL, species=NULL)

for (s in species) {
	setwd(paste("~/Documents/scratch_folders", s, "GLMs_paper", sep="/"))
	d1 <- read.table("bestglmSummary_parental.txt", header=TRUE, sep="\t", as.is=TRUE)
	df.bestglm <- rbind(df.bestglm, d1)
	
	d2 <- read.table("betaCoefficients_parentalGLM.txt", header=TRUE, sep="\t", as.is=TRUE)
	df.pred <- rbind(df.pred, d2)

}

df.bestglm$p.value.label <- sapply(df.bestglm$p.value, function(x) significance.label(x))
df.bestglm[,c(2,3,5,6)] <- apply(df.bestglm[,c(2,3,5,6)], 2, function(x) round(x, 4))

df.pred$p.value.label <- sapply(df.pred$p.value, function(x) significance.label(x))
df.pred[,c(2,3,5,6)] <- apply(df.pred[,c(2,3,5,6)], 2, function(x) round(x,4))

# output
##########################################################################################

setwd("~/Documents/scratch_folders/GLMs_paper")

write.table(df.bestglm, "summary_bestglmSummary_parental.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# glm summary (prediction data)
write.table(df.pred, "summary_betaCoefficients_parentalGLM.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)
