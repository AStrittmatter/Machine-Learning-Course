##############################################################################
########################  Load Packages and the Data  ########################
##############################################################################

### Load the packages  
library(fBasics)     # use for descriptive statistics
library(tidyverse)   # use for handling data
library(caret)       # use for handling data
library(lmtest)      # use for heteroscedasticity robust standard errors
library(sandwich)    # use for heteroscedasticity robust standard errors
library(hdm)         # use for Lasso and Post-Double-Selection
library(glmnet)      # use for lasso and Elastic Net regularized Generalized Linear Models
options(warn=-1)     # supress warnings

print('All packages successfully installed and loaded.')

### Load the Data
set.seed(12345678) 
df <- read.csv("job_corps.csv",header=TRUE, sep=",") # load data from csv-file
df <- df[sample(c(1:nrow(df)), size=3000, replace =F),] # Select a random subsample of 3000 observations
print('Data successfully loaded.')

##############################################################################

##############################################################################
########################  Descriptive Statistics  ############################
##############################################################################

## Table with Descriptive Statistics 
desc <- fBasics::basicStats(df) %>% t() %>% as.data.frame() %>% 
          select(Mean, Stdev, Minimum, Maximum, nobs)
print(round(desc, digits=2))

##############################################################################

#########################################################################
########################  Univariate OLS Regression #####################
#########################################################################

## Univariate OLS
ols1 <- lm(EARNY4 ~ participation, data = df)
summary(ols1)

## Store results
results <- as.matrix(coef(summary(ols1))[2, c("Estimate", "Std. Error", "Pr(>|t|)")])

# Prepare matrix to store results
res <- matrix(NA,nrow=3,ncol=5)
colnames(res) <- c("Univariate OLS", "Multivariate OLS1", "Multivariate OLS2",
                   "Multivariate OLS3", "Multivariate OLS4")
rownames(res) <- rownames(results)
res[,1] <- results

print(round(res[,1], digits=2))

########################################################################

########################################################################
########################  Standardized Differences #####################
########################################################################

## Means and standard deviations for the participants (D=1)
desc_1 <- fBasics::basicStats(df[df$participation==1,]) %>% t() %>% as.data.frame() %>% select(Mean, Stdev)

## Means and standard deviations for the non-participants (D=0)
desc_0 <- fBasics::basicStats(df[df$participation==0,]) %>% t() %>% as.data.frame() %>% select(Mean, Stdev)

# Make table and add standardized differences
desc <- cbind(desc_1[-c(1:3),],desc_0[-c(1:3),], 
        100*abs(desc_1[-c(1:3),1]-desc_0[-c(1:3),1])/sqrt(0.5*(desc_1[-c(1:3),2]^2+desc_0[-c(1:3),2]^2)))
colnames(desc) <- c("D=1 Mean", "D=1 Std.Dev.", "D=0 Mean", "D=0 Std.Dev.", "Std.Diff.")
print(round(desc, digits=2))

########################################################################

#########################################################################
########################  Multivariate OLS Regression ###################
#########################################################################

## Multivariate OLS
ols2 <- lm(EARNY4 ~ participation + age_1 + age_3 + livespou + publich, data = df)
summary(ols2)
# Question: Why do we omit age_2?

## Store results
results <- as.matrix(coef(summary(ols2))[2, c("Estimate", "Std. Error", "Pr(>|t|)")])
res[,2] <- results
print(round(res[,c(1:2)], digits=2))

## Relative change in the estimated effect
print(paste0("Relative change in the estimated effect: ",round(100*(res[1,2]-res[1,1])/res[1,1], digits=1),"%"))

########################################################################

#########################################################################

## Multivariate OLS
ols3 <- lm(EARNY4 ~ ???, data = df)
summary(ols3)

## Store results
results <- as.matrix(coef(summary(ols3))[2, c("Estimate", "Std. Error", "Pr(>|t|)")])
res[,3] <- results
print(round(res[,c(1:3)], digits=2))

## Relative change in the estimated effect
print(paste0("Relative change in the estimated effect: ",round(100*(res[1,3]-res[1,2])/res[1,2], digits=1),"%"))

########################################################################

###############################################################################

## Generate first-order interactions between all control variables
interactions <- t(apply(df[,-c(1,2,3,6,11)], 1, combn, 2, prod))
colnames(interactions) <- paste("Inter.V", combn(1:ncol(df[,-c(1,2,3,6,11)]), 2, paste, collapse="V"), sep="")
print(paste0("Maximm number of interaction terms: ", ncol(interactions)))

## Merge basline characteristics with interaction terms
df_merge <- as.data.frame(cbind(df[,-c(1,2,3,6,11)], interactions))

## Eliminate collinear variables
df2 = cor(df_merge)
df2[is.na(df2)] <- 1
hc = findCorrelation(df2, cutoff=0.8) # putt any value as a "cutoff" 
hc = sort(hc)
df_int = cbind(df[,c(1,3)],df_merge[,-c(hc)])
print(paste0("Total number of control variables: ", ncol(df_int)-2))

###############################################################################

###############################################################################

## Multivariate OLS with all baseline characteristics and interaction terms
ols4 <- lm(EARNY4 ~ ., data = df_int)

## Store results
results <- as.matrix(coef(summary(ols4))[2, c("Estimate", "Std. Error", "Pr(>|t|)")])
res[,4] <- results
print(round(res[,c(1:4)], digits=2))

## Relative change in the estimated effect
print(paste0("Relative change in the estimated effect: ",round(100*(res[1,4]-res[1,3])/res[1,3], digits=1),"%"))

########################################################################

###############################################################################

# Set starting value for replicability 
set.seed(123456) 

# Specify number of random variables
cols <- 1000

# Generate random variables
redundant_x <- matrix(rnorm(nrow(df_int)*cols), nrow = nrow(df_int)) # We draw from a random standard normal distribution
colnames(redundant_x) <- paste("Rand.", 1:cols, sep="")

# Merge random variables with baseline characteritics and interaction terms
df_rand <- as.data.frame(cbind(df_int, redundant_x))
print(paste0("Total number of control variables: ", ncol(df_rand)-2))

###############################################################################

###############################################################################

## Multivariate OLS with all baseline characteristics, interaction terms, and random variables
ols5 <- lm(EARNY4 ~ ., data = df_rand)

## Store results
results <- as.matrix(coef(summary(ols5))[2, c("Estimate", "Std. Error", "Pr(>|t|)")])
res[,5] <- results
print(round(res, digits=2))

## Relative change in the estimated effect
print(paste0("Relative change in the estimated effect: ",round(100*(res[1,5]-res[1,4])/res[1,4], digits=1),"%"))

########################################################################

###############################################################################
########################### Earnings Equation #################################
###############################################################################

# Predict earnings
N <- nrow(df)
st1 <- rlasso(as.matrix(df[,c(4:ncol(df))]), as.matrix(df$EARNY4), 
              penalty = list(homoscedastic = FALSE, c= 1.1, gamma = 0.1/log(N)))
summary(st1)

# Store selected variables
n1<- names(st1$coefficients[(st1$coefficients != 0) == TRUE])[-1]

###############################################################################

###############################################################################
######################### Participation Probability ###########################
###############################################################################

# Predict participation
N <- nrow(df)
st2 <- rlasso(as.matrix(df[,c(4:ncol(df))]), as.matrix(df$participation), 
              penalty = list(homoscedastic = FALSE, c= 1.1, gamma = 0.1/log(N)))
summary(st2)

# Store selected variables
n2<- names(st2$coefficients[(st2$coefficients != 0) == TRUE])[-1]

###############################################################################

###############################################################################
################################# Post-Lasso ##################################
###############################################################################

# Take union of selected covariates
selected_covariates <- c("participation", unique(c(n1, n2)))

# Setup the formula of the linear regression model
sumx <- paste(selected_covariates, collapse = " + ")  
linear <- paste("EARNY4",paste(sumx, sep=" + "), sep=" ~ ")
linear <- as.formula(linear)

# Post-Lasso regression
ols <- lm(linear, data = df)
summary(ols)

# Heteroskedasticity robust standard errors
#coeftest(ols, vcov = vcovHC(ols, type = "HC1"))

###############################################################################

###############################################################################
################## Estimate the Treatment Effect Directly #####################
###############################################################################

# Post-Double-Selection Procedure 
dsp <- rlassoEffect(as.matrix(df[,c(4:ncol(df))]), as.matrix(df$EARNY4)
          , as.matrix(df$participation), model = TRUE, penalty = list(homoscedastic = FALSE), method = "double selection")
summary(dsp)

###############################################################################
# Earning Equation
###############################################################################

# Predict earnings

# Store selected variables

###############################################################################
# Participation Probability
###############################################################################

# Predict participation

# Store selected variables

###############################################################################
# Post-Lasso Model
###############################################################################

# Take union of selected covariates
selected_covariates <- c("participation", unique(c(n1, n2)))

# Setup the formula of the linear regression model
sumx <- paste(selected_covariates, collapse = " + ")  
linear <- paste("EARNY4",paste(sumx, sep=" + "), sep=" ~ ")
linear <- as.formula(linear)

# Post-Lasso OLS regression
ols <- lm(linear, data = df_rand)
summary(ols)

###############################################################################

####################################################################
################# Cross-Validated Lasso ############################
####################################################################

set.seed(123456789) # Starting value

# Cross-validated Lasso in earnings equation
lasso_earn <- cv.glmnet(as.matrix(df_int[,c(3:ncol(df_int))]), as.matrix(df$EARNY4), 
                        alpha=1, nfolds = 10, type.measure = 'mse', standardize = TRUE)
# alpha =1 is Lasso, alpha = 0 is Ridgde
# nfolds - number of cross-validation folds
# type.measure - measure for model accuracy

plot(lasso_earn)

####################################################################

####################################################################

# Plot Lasso coefficients
coef(lasso_earn,s = lasso_earn$lambda.1se) 
# $lambda.min - Lambda that minimizes cross-validated MSE
# $lambda.1se - Lambda of 1 standard error rule

####################################################################

####################################################################

# Select covariates with non-zero coefficients
coef <- predict(lasso_earn,s = lasso_earn$lambda.min, type = "nonzero") #
colnames <- colnames(df_int[,c(3:ncol(df_int))])
n1 <- colnames[unlist(coef)]
print(paste0("Number of Selected Variables Earnings Equation: ",length(n1)))
print("Selected Variables:")
print(n1)

####################################################################

####################################################################

set.seed(123456789) # Starting value

# Cross-validated Lasso in participation equation
lasso_part <- cv.glmnet(???, 
                        alpha=1, nfolds = 10, type.measure = 'mse', standardize = TRUE)
plot(lasso_part)

####################################################################

####################################################################

# Select covariates with non-zero coefficients
coef <- predict(???,s = ???, type = "nonzero") #
colnames <- colnames(df_int[,c(3:ncol(df_int))])
print(paste0("Number of Selected Variables Participation Equation: ",length(n2)))
print("Selected Variables:")
print(n2)

####################################################################

###############################################################################
# Post-Lasso Model
###############################################################################

# Take union of selected covariates
selected_covariates <- c(???)

# Setup the formula of the linear regression model
sumx <- paste(selected_covariates, collapse = " + ")  
linear <- paste("EARNY4",paste(sumx, sep=" + "), sep=" ~ ")
linear <- as.formula(linear)

# Post-Lasso OLS regression
ols <- lm(linear, data = df_int)
summary(ols)

###############################################################################
