###############################################################################
## Course: Machine Learning for Economists and Business Analysts
## Topic: High-Dimensional Confounding
###############################################################################

#rm(list = ls())
set.seed(100239) 

#getwd()
#setwd("")

#  Load Packages  
library("fBasics")
library("glmnet")
library("AER")
library("hdm")
library("lmtest")
library("sandwich")
library("tidyverse")

# Load data 
df <- read.csv("job_corps.csv",header=TRUE, sep=",")

########################################
### Exercise 1: Conventional Methods ###
########################################

###############################################################################
# Task 1
###############################################################################

# Table with Descriptive Statistics 
desc <- fBasics::basicStats(df) %>% t() %>% as.data.frame() %>% 
  select(Mean, Stdev, Minimum, Maximum, nobs)
print(round(desc, digits=2))

###############################################################################
# Task 2
###############################################################################

# Univariate OLS
ols <- lm(EARNY4 ~ participation, data = df)
summary(ols)

# Robust standard errors
coeftest(ols, vcov = vcovHC(ols, type = "HC1"))

###############################################################################
# Task 3
###############################################################################

# IV results
iv <- ivreg(formula = EARNY4 ~ participation | assignment, data = df)
summary(iv)

# Robust standard errors
coeftest(iv, vcov = vcovHC(ols, type = "HC1"))

##############################################
### Exercise 2: Double Selection Procedure ###
##############################################

###############################################################################
# Task 1
###############################################################################

# Predict earnings
st1 <- rlasso(as.matrix(df[,c(4:ncol(df))]), as.matrix(df$EARNY4))
summary(st1)
# Store selected variables
n1<- names(st1$coefficients[(st1$coefficients != 0) == TRUE])[-1]

###############################################################################
# Task 2
###############################################################################

# Predict participation
st2 <- rlasso(as.matrix(df[,c(4:ncol(df))]), as.matrix(df$participation))
summary(st2)
# Store selected variables
n2<- names(st2$coefficients[(st2$coefficients != 0) == TRUE])[-1]

###############################################################################
# Task 3
###############################################################################

# Take uniion of selected covariates
selected_covariates <- c("participation", unique(c(n1, n2)))

# Setup the formula of the linear regression model
sumx <- paste(selected_covariates, collapse = " + ")  
linear <- paste("EARNY4",paste(sumx, sep=" + "), sep=" ~ ")
linear <- as.formula(linear)

# Post-Lasso OLS regression
ols <- lm(linear, data = df)
summary(ols)

# Robust standard errors
coeftest(ols, vcov = vcovHC(ols, type = "HC1"))

###############################################################################
# Task 4
###############################################################################

# Double Selection Procedure 
dsp <- rlassoEffect(as.matrix(df[,c(4:ncol(df))]), as.matrix(df$EARNY4)
          , as.matrix(df$participation), model = TRUE, method = "double selection")
summary(dsp)

###########################################
### Exercise 3: Double Machine Learning ###
###########################################

###############################################################################
# Task 1
###############################################################################
set.seed(1001)

# Partition Samples for Cross-Fitting
df_part <- modelr::resample_partition(df, c(obs_A = 0.5, obs_B = 0.5))
df_obs_A <- as.data.frame(df_part$obs_A) # Sample A
df_obs_B <- as.data.frame(df_part$obs_B) # Sample B

##  Generate Variables  
# Outcome
earnings_obs_A <- as.matrix(df_obs_A[,1])
earnings_obs_B <- as.matrix(df_obs_B[,1])

# Treatment
treat = 3 #Select treatment 2= offer to participate, 3 = actual participation
treat_obs_A <- as.matrix(df_obs_A[,treat])
treat_obs_B <- as.matrix(df_obs_B[,treat])

# Covariates
covariates_obs_A <- as.matrix(df_obs_A[,c(4:ncol(df_obs_A))])
covariates_obs_B <- as.matrix(df_obs_B[,c(4:ncol(df_obs_B))])

###############################################################################
# Task 2
###############################################################################

##  Conditional Potential Earnings under Non-Participation
p = 1 # 1 for LASSO, 0 for Ridge
set.seed(100237)

## Using Sample A to Predict Sample B
# Potential Earnings under Non-Treatment
lasso_y0_A <- cv.glmnet(covariates_obs_A[treat_obs_A==0,], earnings_obs_A[treat_obs_A==0,],
                              alpha=p, type.measure = 'mse')
plot(lasso_y0_A)
fit_y0_A <- glmnet(covariates_obs_A[treat_obs_A==0,], earnings_obs_A[treat_obs_A==0,]
                        ,lambda = lasso_y0_A$lambda.min)
y0hat_B <- predict(fit_y0_A, covariates_obs_B)

## Using Sample B to Predict Sample A
# Potential Earnings under Non-Treatment
lasso_y0_B <- cv.glmnet(covariates_obs_B[treat_obs_B==0,], earnings_obs_B[treat_obs_B==0,],
                              alpha=p, type.measure = 'mse')
plot(lasso_y0_B)
fit_y0_B <- glmnet(covariates_obs_B[treat_obs_B==0,], earnings_obs_B[treat_obs_B==0,]
                        ,lambda = lasso_y0_B$lambda.min)
y0hat_A <- predict(fit_y0_B, covariates_obs_A, type = 'response')

###############################################################################
# Task 3
###############################################################################

##  Propensity Score  
p = 1 # 1 for LASSO, 0 for Ridge
set.seed(100236)

# Using Sample A to Predict Sample B
lasso_p_A <- cv.glmnet(covariates_obs_A, treat_obs_A, alpha=p, type.measure = 'mse')
plot(lasso_p_A)
fit_p_A <- glmnet(covariates_obs_A, treat_obs_A,lambda = lasso_p_A$lambda.min)
pscore_B <- predict(fit_p_A, covariates_obs_B)

# Using Sample B to Predict Sample A
lasso_p_B <- cv.glmnet(covariates_obs_B, treat_obs_B, alpha=p, type.measure = 'mse')
plot(lasso_p_B)
fit_p_B <- glmnet(covariates_obs_B, treat_obs_B,lambda = lasso_p_B$lambda.min)
pscore_A <- predict(fit_p_B, covariates_obs_A)

###############################################################################
# Task 4
###############################################################################

## Efficient Score
p_A = mean(pscore_A)
p_B = mean(pscore_B)

# Generate Modified Outcome in each sample
Y_star_A = invisible(treat_obs_A*(earnings_obs_A - y0hat_A)/p_A 
            - (1-treat_obs_A)*pscore_A*(earnings_obs_A - y0hat_A)/(p_A*(1-pscore_A)))

Y_star_B = invisible(treat_obs_B*(earnings_obs_B - y0hat_B)/p_B 
            - (1-treat_obs_B)*pscore_B*(earnings_obs_B - y0hat_B)/(p_B*(1-pscore_B)))

Y_star <- 0.5*(mean(Y_star_A) + mean(Y_star_B))

###############################################################################
# Task 5
###############################################################################

# Average Treatment Effect for Treated (ATET)
ATET <- round(Y_star, digits=3)

# Estimate variance for each sample
var_A = mean((Y_star_A- treat_obs_A*Y_star/p_A)^2) 
var_B = mean((Y_star_B - treat_obs_B*Y_star/p_B)^2) 

# Split sample estimator for standard error
SD_ATET <- round(sqrt(0.5*(var_A + (mean(Y_star_A) - Y_star)^2 
                     + var_B + (mean(Y_star_B) - Y_star)^2)
                     /(length(Y_star_A) + length(Y_star_B))),digits=3)

print(paste0("Average Treatment Effect for Treated (ATET): ", ATET))
print(paste0("Standard Error for ATET: ", SD_ATET))


