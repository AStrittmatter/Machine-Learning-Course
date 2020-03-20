###############################################################################
## Course: Machine Learning for Economists and Business Analysts
## Topic: Effect Heterogeneity
###############################################################################

rm(list = ls())
set.seed(100239) 

#getwd()
#setwd("")

#  Load Packages  
library("fBasics")
library("glmnet")
library("AER")
library("grf")
library("hdm")
library("lmtest")
library("sandwich")
library("tidyverse")

# Load data 
df <- read.csv("job_corps.csv",header=TRUE, sep=",")

###########################################
### Exercise 3: Double Machine Learning ###
###########################################

######################
## Data Preparation ##
######################
set.seed(123456789)

# Generate variable with the rows in training data
size <- floor(0.5 * nrow(df))
set_A <- sample(seq_len(nrow(df)), size = size)
set_B <- seq_len(nrow(df))[-set_A] 

##  Generate Variables  
# Outcome
earnings <- as.matrix(df[,1])

# Treatment
treat = 2 #Select treatment 2= offer to participate, 3 = actual participation
treat <- as.matrix(df[,treat])

# Covariates
covariates <- as.matrix(df[,c(4:ncol(df))])

#########################
## Nuisance Parameters ##
#########################

###############################################################################

##  Conditional Potential Earnings under Non-Treatment
p = 1 # 1 for LASSO, 0 for Ridge
set.seed(100237)

## Using Sample A to Predict Sample B
# Potential Earnings under Non-Treatment
lasso_y0_A <- cv.glmnet(covariates[c(set_A,treat==0),], earnings[c(set_A,treat==0)],
                              alpha=p, type.measure = 'mse')
plot(lasso_y0_A)
fit_y0_A <- glmnet(covariates[c(set_A,treat==0),], earnings[c(set_A,treat==0)]
                        ,lambda = lasso_y0_A$lambda.min)
y0hat_B <- predict(fit_y0_A, covariates)

## Using Sample B to Predict Sample A
# Potential Earnings under Non-Treatment
lasso_y0_B <- cv.glmnet(covariates[c(set_B,treat==0),], earnings[c(set_B,treat==0)],
                              alpha=p, type.measure = 'mse')
plot(lasso_y0_B)
fit_y0_B <- glmnet(covariates[c(set_B,treat==0),], earnings[c(set_B,treat==0)]
                        ,lambda = lasso_y0_B$lambda.min)
y0hat_A <- predict(fit_y0_B, covariates)

###############################################################################

##  Conditional Potential Earnings under Treatment
p = 1 # 1 for LASSO, 0 for Ridge
set.seed(100237)

## Using Sample A to Predict Sample B
# Potential Earnings under Treatment
lasso_y1_A <- cv.glmnet(covariates[c(set_A,treat==1),], earnings[c(set_A,treat==1)],
                        alpha=p, type.measure = 'mse')
plot(lasso_y1_A)
fit_y1_A <- glmnet(covariates[c(set_A,treat==1),], earnings[c(set_A,treat==1)]
                   ,lambda = lasso_y1_A$lambda.min)
y1hat_B <- predict(fit_y1_A, covariates)

## Using Sample B to Predict Sample A
# Potential Earnings under Treatment
lasso_y1_B <- cv.glmnet(covariates[c(set_B,treat==1),], earnings[c(set_B,treat==1)],
                        alpha=p, type.measure = 'mse')
plot(lasso_y1_B)
fit_y1_B <- glmnet(covariates[c(set_B,treat==1),], earnings[c(set_B,treat==1)]
                   ,lambda = lasso_y1_B$lambda.min)
y1hat_A <- predict(fit_y1_B, covariates, type = 'response')


###############################################################################

##  Propensity Score  
p = 1 # 1 for LASSO, 0 for Ridge
set.seed(100236)

# Using Sample A to Predict Sample B
lasso_p_A <- cv.glmnet(covariates[set_A,], treat[set_A], alpha=p, type.measure = 'mse')
plot(lasso_p_A)
fit_p_A <- glmnet(covariates[set_A,], treat[set_A],lambda = lasso_p_A$lambda.min)
pscore_B <- predict(fit_p_A, covariates)

# Using Sample B to Predict Sample A
lasso_p_B <- cv.glmnet(covariates[set_B,], treat[set_B,], alpha=p, type.measure = 'mse')
plot(lasso_p_B)
fit_p_B <- glmnet(covariates[set_B,], treat[set_B,],lambda = lasso_p_B$lambda.min)
pscore_A <- predict(fit_p_B, covariates)

####################################
## Average Treatment Effect (ATE) ##
####################################

## Efficient Score
# Generate Modified Outcome in each sample
Y_star <- matrix(NA,nrow=nrow(df),ncol=1)
Y_star[set_A] <- invisible(y1hat_A[set_A] -y0hat_A[set_A] 
          + treat[set_A]*(earnings[set_A]-y1hat_A[set_A])/pscore_A[set_A]
          - (1-treat[set_A])*(earnings[set_A]-y0hat_A[set_A])/(1-pscore_A[set_A]))
Y_star[set_B] <- invisible(y1hat_B[set_B] -y0hat_B[set_B] 
          + treat[set_B]*(earnings[set_B]-y1hat_B[set_B])/pscore_B[set_B]
          - (1-treat[set_B])*(earnings[set_B]-y0hat_B[set_B])/(1-pscore_B[set_B]))

# Average Treatment Effect (ATE)
ATE <- round(mean(Y_star), digits=2)
se_ATE <- round(sd(Y_star)/sqrt(nrow(df)), digits=2)

print(paste0("Average Treatment Effect (ATE): ", ATE))
print(paste0("Standard Error for ATE: ", se_ATE))

##########################
## Effect Heterogeneity ##
##########################

set.seed(100237)

## Predict Effect Heterogeneity
lasso_A <- cv.glmnet(covariates[set_A,], Y_star[set_A],
                        alpha=p, type.measure = 'mse')
plot(lasso_A)
fit_A <- glmnet(covariates[set_A,], Y_star[set_A] ,lambda = lasso_A$lambda.min)
coef(fit_A)
# Extrapolate to sample B
het_B <- predict(fit_A, covariates)

## Predict Effect Heterogeneity
lasso_B <- cv.glmnet(covariates[set_B,], Y_star[set_B],
                   alpha=p, type.measure = 'mse')
plot(lasso_B)
fit_B <- glmnet(covariates[set_B,], Y_star[set_B],lambda = lasso_B$lambda.min)
coef(fit_B)
# Extrapolate to sample B
het_A <- predict(fit_B, covariates)

het_dml <- matrix(NA, nrow = nrow(df),ncol =1)
het_dml[set_A] <- het_A[set_A]
het_dml[set_B] <- het_B[set_B]

# Kernel Density Plot
d_dml <- density(het_dml) 
plot(d_dml) 

##################
## Post DML-OLS ##
##################

# Multivariate OLS
ols <- lm(Y_star ~ ., data = df[,-c(1,2,3,6,11)])
summary(ols)

# Robust standard errors
coeftest(ols, vcov = vcovHC(ols, type = "HC1"))

###################
## Causal Forest ##
###################

set.seed(1234567)
cf <- causal_forest(covariates, earnings, treat)
het_cf <- predict(cf,estimate.variance = TRUE)

# Kernel Density Plot
d_cf <- density(het_cf$predictions) 
plot(d_cf) 

cor(het_dml,het_cf$predictions)

## Inference

# t-Statistics
t_stat <- as.matrix(het_cf$predictions)/ as.matrix(sqrt(het_cf$variance.estimates))

sig_pos <- (t_stat>=1.96)== TRUE
sig_neg <- (t_stat<=-1.96)== TRUE
insig <- (abs(t_stat) <1.96)== TRUE

print(paste0("Share with positive effects: ", round(mean(sig_pos), digits=4)))
print(paste0("Share with negative effects: ", round(mean(sig_neg), digits=4)))
print(paste0("Share with insignificant effects: ", round(mean(insig), digits=4)))




