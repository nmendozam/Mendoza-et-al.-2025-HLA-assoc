# This script performs logistic regression analysis on immunogenicity data
# to determine if any immunogenicity score against a set of epitopes
# is associated with a negative effect on longevity in a given population.


##########################################
# Read the predicted immunogenicity scores
##########################################
drb = read.csv("DRB1.csv", row.names="ID")
drb1 = read.csv("DRB1_1.csv", row.names="ID")

# Exclude samples that are not in both files (one allele is unknown)
library(compare)
IDs <- intersect(row.names(drb),row.names(drb1))
drb <- drb[row.names(drb) %in% IDs,]
drb1 <- drb1[row.names(drb1) %in% IDs,]

# Scale the variables
drb$LLI <- drb$LLI - 1

drb[,4:ncol(drb)] <- drb[,4:ncol(drb)] + drb1 [,4:ncol(drb1)]  / 2

# Exclude DRB*15:01 to see if it affects the association downstream
# (keep this line commented for the global association)
# drb <- drb[!(drb$DRB1 == 'DRB1*15:01' | drb$DRB1_1 == 'DRB1*15:01'), ]

#########################################
# Map the epitopes to the source proteins
#########################################

# Get the file with the protein source of the epitopes
sele <- readxl::read_excel("Selection of epitopes.xlsx")
sele$Source <- tolower(sele$Source)
sele$Epitope <- sapply(strsplit(unlist(sele$Epitope), " "), head, 1)

# Merge epitopes by source protein (get min value)
by_prot <- data.frame(LLI = drb$LLI, row.names = row.names(drb))
for (prot_name in unique(sele$Source)) {
  # prot_name <- "apolipoprotein b-100"
  print(prot_name)
  prot_list <- unlist(sele[sele$Source == prot_name, "Epitope"], use.names=FALSE)
  if (length(prot_list) > 1) {
    summary_values <- apply(drb[,prot_list], 1, max)
  } else {
    summary_values <- drb[,prot_list]
  }
  by_prot[prot_name] <- summary_values
}

#######################
# Logistic regression #
#######################

## Remove colinear variables
logit <- glm(LLI ~ ., data = by_prot, family = binomial)
library(car)
vif(logit) < 12
by_prot_no_vif <- by_prot[!(colnames(by_prot) %in% c("Calreticulin"))]
logit <- glm(LLI ~ ., data = by_prot_no_vif, family = binomial)
summary(logit)

# Get principal components
pcs <- read.csv("GSAdata_raw.merge.HapMapIII_CGRCh37.eigenvec", sep=" ", header=FALSE)[,2:5]
colnames(pcs) <- c("ID", "PC1", "PC2", "PC3")
row.names(pcs) <- pcs$ID
pcs <- subset(pcs, select = -c(ID))
# Append them to the data set to correct for population stratification
byprot_with_pcs = merge(by_prot_no_vif, pcs, by="row.names", all=TRUE) 
row.names(byprot_with_pcs) <- byprot_with_pcs$Row.names
byprot_with_pcs <- subset(byprot_with_pcs, select = -c(Row.names))
byprot_with_pcs <- na.omit(byprot_with_pcs)


# Get sex from ID
byprot_with_pcs$sex <- rep(1, nrow(byprot_with_pcs))
byprot_with_pcs$sex[grep('_m_', row.names(byprot_with_pcs))] <- 0

# Run the association
logit <- glm(LLI ~ ., data = byprot_with_pcs, family = binomial) 
# This is ran to get the model evaluation metrics, however the statistics
# are taken from the robust model below

# Evaluate the model
library(lmtest)
library(performance)
check_model(logit)
bptest(logit)

################################
# MLR robust logistic regression
################################

library(lavaan)

names(byprot_with_pcs) <- make.names(names(byprot_with_pcs))
model2 = 'LLI ~ myelin.basic.protein + cofilin.1 + apolipoprotein.b.100 + PC1 + PC2 + PC3 + sex'
       
logit.fit.robust = sem(model2, data=byprot_with_pcs, estimator = "MLR")
summary(logit.fit.robust)

lavInspect(logit.fit.robust, "scaling.factor")

logit.fit <- sem(model2, data = byprot_with_pcs, estimator = "ML")  # classical SE

# Standard estimates
coef_classical <- parameterEstimates(logit.fit, standardized = TRUE)

# Robust estimates
coef_robust <- parameterEstimates(logit.fit.robust, standardized = TRUE)

# Compare
comparison <- merge(
  coef_classical[, c("lhs","op","rhs","est","se","pvalue")],
  coef_robust[, c("lhs","op","rhs","est","se","pvalue")],
  by = c("lhs","op","rhs"),
  suffixes = c("_classical", "_robust")
)

print(comparison)


library(arm)
preds <- predict(logit.fit.robust)  # same as logit.fit.robust

# Raw residuals
res <- residuals(logit.fit.robust)
binnedplot(preds, res, xlab = "Predicted Probability",
           ylab = "Residuals", main = "Binned Residuals Plot")
