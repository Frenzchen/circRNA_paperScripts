#!/sur/bin/env Rscript

# R 3.1
#
# Table 1
#
# Reads data table for variables (non-normalised) for GLM of human parental genes
# Script performs the following tasks:
#
# 1. Normalisation of data table   
# 2. Partitioning into training (80%) and prediction set (20%)
# 3. Run bestglm to determine significance of variables
# 4. Run glm on prediction data with formula of bestglm
# 5. Get model statistics
#
# Outputs
# 1. tab-delimited file with stas of bestglm run (bestglmSummary_parental.txt )
# 2. tab-delimited file with stas of validation run (betaCoefficients_parentalGLM.txt)
# 3. tab-delimited file with confusion matrix (modelPerformance_confMatrix.txt)
# 4. tab-delimited file with true and predicted values (predictions_parentalGLM.txt)


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
	colnames(cf)[c(2:6)] <- c("coefficient", "std.error", "p-value", "lower.ci", "upper.ci")
	return(cf[-1,])
}

# Calculates confusion matrix for model
# Based on the 0.25 quantile as prediction border for positive/negative rates
#
# Arguments
# x: data table with ture and predicted values
#
# Returns
# v: vector with accuracy, sensitivity and specificity of model
modelperformance <- function(x) {
	# accuracy = tp
	# sensitivity = tp/tp+fn
	# specificity = tn/tn+fp
	
	q25 <- quantile(x$pred[x$trueVal==1], 0.25)
	
	tp <- nrow(subset(x, (trueVal == 1 & pred >= q25)))
	fp <- nrow(subset(x, (trueVal == 0 & pred >= q25)))
	tn <- nrow(subset(x, (trueVal == 0 & pred < q25)))
	fn <- nrow(subset(x, (trueVal == 1 & pred < q25)))

	acc <- (tp + tn)/nrow(x)
	sens <- tp / (tp + fn)
	spec <- tn / (tn + fp)
	
	v <- c(acc, sens, spec)
	
	return(v)
}


# run
##########################################################################################

# 1. Normalisation
df <- read.table("geneTable_human_glm.txt", sep="\t", header=TRUE, as.is=TRUE)
#df[,c(2:5,8:9,13:15)] <- scale(df[,c(2:5,8:9,13:15)])

# 2. Data partitioning
tr <- createDataPartition(df$circle, p = .8, list = FALSE, times = 1)
training.set <- df[tr, c(2:5,8:9,13:16)]
nrow(training.set)
prediction.set <- df[-tr, c(2:5,8:9,13:16)]
nrow(prediction.set)

print("VIFs")
vifx(training.set)

# 3. run bestglm on training set, extract forumla and run on prediction set
bst <- bestglm(training.set, IC="CV", t=1000, family=binomial)
bst.summary <- summarize.model(bst)
bst.summary$species <- "human"

# 4. get variables and predict on prediction set
variables <- names(coef(bst$BestModel))
variables <- variables[2:length(variables)]
f <- as.formula(paste("circle~", paste(variables, collapse="+")))

glm.pred <- glm(data=prediction.set, formula=f, family="binomial")

glm.pred.coef <- beta(glm.pred)$coefficients
glm.pred.coef <- data.frame(glm.pred.coef[-1,c(1,2,4)])
colnames(glm.pred.coef) <- c("beta.coefficient", "std.error", "p-value")
glm.pred.coef$lower.ci <- glm.pred.coef$beta.coefficient - 1.96 * glm.pred.coef$std.error
glm.pred.coef$upper.ci <- glm.pred.coef$beta.coefficient + 1.96 * glm.pred.coef$std.error
glm.pred.coef$species <- "human"
glm.pred.coef$predictor <- sub(".z", "", rownames(glm.pred.coef))
glm.pred.coef <- glm.pred.coef[,c(7, 1:6)]

# 5. get stats on model performance
pred <- data.frame(trueVal=prediction.set$circle, pred=glm.pred$fitted.values, rows=as.numeric(row.names(prediction.set)), species="human")

# confusion matrix
d <- modelperformance(pred)
modelperf <- data.frame(acc=d[1], sens=d[2], spec=d[3])


# output
##########################################################################################

# best glm summary (training data)
write.table(bst.summary, "bestglmSummary_parental.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# glm summary (prediction data)
write.table(glm.pred.coef, "betaCoefficients_parentalGLM.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# confusion matrix of prediction data
write.table(modelperf, "modelPerformance_confMatrix.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

# predicted vs real data
write.table(pred, "predictions_parentalGLM.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)

