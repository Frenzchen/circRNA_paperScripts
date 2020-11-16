#!/usr/bin/env Rscript

# R 3.1
#
# Supplementary Table 8
#
# Reads data table for variables (non-normalised) for GLM of mouse parental genes
# Script performs the following tasks:
#
# 1. Normalisation of data table   
# 2. Partitioning into training (80%) and prediction set (20%)
# 3. Run bestglm to determine significance of variables
# 4. Run glm on prediction data with formula of bestglm
#
# Outputs
# 1. tab-delimited file with confusion matrix (modelPerformance_confMatrix.txt)


set.seed(22121987)

library(bestglm) # v0.37
library(caret) # v6.0.84
library(reghelper) #v1.0.0


# function
##########################################################################################

# Creates summary table for model statistics (coefficient, std error, p, CIs)
#
# Arguments
# x: glm model
#
# Returns
# cf: summary data table
summarize.model <- function(x) {
	civ <- data.frame(confint(x$BestModel))
	civ$predictor <- rownames(civ)
	cf <- data.frame(coef(summary(x$BestModel))[,c(1,2,4)])
	cf$predictor <- rownames(cf)
	cf <- merge(cf, civ, by="predictor")
	colnames(cf)[c(2:6)] <- c("coefficient", "std.error", "p.value", "lower.ci", "upper.ci")
	return(cf[-1,])
}


# run
##########################################################################################

species <- c("opossum", "mouse", "rat", "rhesus", "human")
bst <- data.frame(predictor=NULL, coefficient=NULL, std.error=NULL, p.value=NULL, lower.ci=NULL, upper.ci=NULL, species=NULL)
pred <- data.frame(predictor=NULL, beta.coefficient=NULL, std.error=NULL, p.value=NULL, lower.ci=NULL, upper.ci=NULL, species=NULL)

for (s in species) {

	if (s == "opossum") {
		v <- c(2:5,8:9,13:14,16)
	} else if (s == "mouse") {
		v <- c(2:5,8:9,13:15,17)
	} else if (s == "rat") {
		v <- c(2:5,8:12,14)
	} else if (s == "rhesus") {
		v <- c(2:5,8:9,13:14,16)
	} else if (s == "human") {
		v <- c(2:5,8:9,13:15, 17)
	}

	# 1. Normalisation
	df <- read.table(paste("~/Documents/scratch_folders/", s, "/GLMs_paper/geneTable_", s, "_glm.txt", sep=""), sep="\t", header=TRUE, as.is=TRUE)

	df <- subset(df, circle == 1)
	df <- subset(df, (age == "therian" | age == s))
	
	print("data size")
	print(nrow(df))
	
	# 2. Data partitioning
	tr <- createDataPartition(df$hs, p = .8, list = FALSE, times = 1)
	training.set <- df[tr, v]
	print(nrow(training.set))
	prediction.set <- df[-tr, v]
	print(nrow(prediction.set))

	print("VIFs")
	print(vifx(training.set))

	# 3. run bestglm on training set, extract forumla and run on prediction set
	bst.sub <- bestglm(training.set, IC="CV", t=1000, family=binomial)
	bst.summary <- summarize.model(bst.sub)
	bst.summary$species <- s
	bst <- rbind(bst, bst.summary)
	
	# 4. get variables and predict on prediction set
	variables <- names(coef(bst.sub$BestModel))
	variables <- variables[2:length(variables)]
	f <- as.formula(paste("hs~", paste(variables, collapse="+")))

	glm.pred <- glm(data=prediction.set, formula=f, family="binomial")

	glm.pred.coef <- beta(glm.pred)$coefficients
	glm.pred.coef <- data.frame(glm.pred.coef[,c(1,2,4)])
	colnames(glm.pred.coef) <- c("beta.coefficient", "std.error", "p.value")
	glm.pred.coef$lower.ci <- glm.pred.coef$beta.coefficient - 1.96 * glm.pred.coef$std.error
	glm.pred.coef$upper.ci <- glm.pred.coef$beta.coefficient + 1.96 * glm.pred.coef$std.error
	glm.pred.coef$species <- s
	glm.pred.coef$predictor <- sub(".z", "", rownames(glm.pred.coef))
	glm.pred.coef <- glm.pred.coef[,c(7, 1:6)]
	
	pred <- rbind(pred, glm.pred.coef)
}

bst[c(2:3,5:6)] <- apply(bst[c(2:3,5:6)], 2, function(x) round(x, 4))
pred <- pred[-grep("Intercept", pred$predictor),]
pred[c(2:3,5:6)] <- apply(pred[c(2:3,5:6)], 2, function(x) round(x, 4))

# output
##########################################################################################

setwd("~/Documents/scratch_folders/GLMs_paper")

# best glm summary (training data)
write.table(bst, "bestglmSummary_hotspots.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# glm summary (prediction data)
write.table(pred, "betaCoefficients_hotspotsGLM.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)
