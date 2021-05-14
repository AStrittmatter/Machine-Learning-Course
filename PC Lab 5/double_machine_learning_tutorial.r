##############################################################################
########################  Load Packages and the Data  ########################
##############################################################################

### Load the packages  
library(fBasics)     # use for descriptive statistics
library(tidyverse)   # use for handling data
library(DiagrammeR)  # use for plotting trees
library(lmtest)      # use for heteroscedasticity robust standard errors
library(sandwich)    # use for heteroscedasticity robust standard errors
library(grf)         # use for generalized random forest
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

###############################################################################
######################### Sample Splitting ####################################
###############################################################################

# Set starting value 
set.seed(123456789)

# Partition Samples for Cross-Fitting
df_part <- modelr::resample_partition(df, c(obs_A = 0.5, obs_B = 0.5)) # Split sample in strata of equal size
df_obs_A <- as.data.frame(df_part$obs_A) # Sample A
df_obs_B <- as.data.frame(df_part$obs_B) # Sample B

##  Generate Variables  
# Outcome variable
earnings_obs_A <- as.matrix(df_obs_A[,1])
earnings_obs_B <- as.matrix(df_obs_B[,1])

# Treatment variable
treat = 3 #Select treatment 2= offer to participate, 3 = actual participation
treat_obs_A <- as.matrix(df_obs_A[,treat])
treat_obs_B <- as.matrix(df_obs_B[,treat])

# Covariates
covariates_obs_A <- as.matrix(df_obs_A[,c(4:ncol(df_obs_A))])
covariates_obs_B <- as.matrix(df_obs_B[,c(4:ncol(df_obs_B))])

print('Sample partitioning ready.')

##############################################################################

###############################################################################
########### Conditional Potential Earnings under Non-Participation ############
###############################################################################

p = 1 # 1 for LASSO, 0 for Ridge

# Set starting value
set.seed(123456789)

# Estimate Lasso among non-participants in Sample A
# Use cross-validation to select optimal lambda value
lasso_y0_A <- cv.glmnet(covariates_obs_A[treat_obs_A==0,], earnings_obs_A[treat_obs_A==0,],
                              alpha=p, type.measure = 'mse')
# Plot the cross-validated MSE
plot(lasso_y0_A)

# Extrapolate the fitted values to Sample B
y0hat_B <- predict(lasso_y0_A, newx = covariates_obs_B, type = 'response', s = lasso_y0_A$lambda.min)

# Estimate Lasso among non-participants in Sample B
lasso_y0_B <- cv.glmnet(covariates_obs_B[treat_obs_B==0,], earnings_obs_B[treat_obs_B==0,],
                              alpha=p, type.measure = 'mse')
# Plot the cross-validated MSE
plot(lasso_y0_B) 

# Extrapolate the fitted values to Sample A
y0hat_A <- predict(lasso_y0_B, newx = covariates_obs_A, type = 'response', s= lasso_y0_B$lambda.min)

# Merge fitted values of both samples
y0hat <- rbind(y0hat_A,y0hat_B)

#################################################################################

###############################################################################
########### Conditional Potential Earnings under Participation ############
###############################################################################

p = 1 # 1 for LASSO, 0 for Ridge

# Set starting value
set.seed(123456789)

# Estimate Lasso among participants in Sample A
# Use cross-validation to select optimal lambda value
lasso_y1_A <- cv.glmnet(covariates_obs_A[treat_obs_A==1,], earnings_obs_A[treat_obs_A==1,],
                              alpha=p, type.measure = 'mse')
plot(lasso_y1_A)

# Extrapolate the fitted values to Sample B
y1hat_B <- predict(lasso_y1_A, newx = covariates_obs_B, type = 'response', s = lasso_y1_A$lambda.min)

# Estimate Lasso among participants in Sample B
lasso_y1_B <- cv.glmnet(covariates_obs_B[treat_obs_B==1,], earnings_obs_B[treat_obs_B==1,],
                              alpha=p, type.measure = 'mse')
plot(lasso_y1_B)

# Extrapolate the fitted values to Sample A
y1hat_A <- predict(lasso_y1_B, newx = covariates_obs_A, type = 'response', s= lasso_y1_B$lambda.min)

# Merge the fitted values of both samples
y1hat <- rbind(y1hat_A,y1hat_B)

#################################################################################

###############################################################################
########################### Propensity Score ##################################
###############################################################################

#  Propensity Score  
p = 1 # 1 for LASSO, 0 for Ridge

# Set starting value
set.seed(123456789)

# Estimate Logit-Lasso in Sample A
# Use cross-validation to select optimal lambda value
lasso_p_A <- cv.glmnet(covariates_obs_A, treat_obs_A, alpha=p, type.measure = 'mse', family="binomial")
plot(lasso_p_A)

# Extrapolate the fitted values to Sample B
pscore_B <- predict(lasso_p_A, newx = covariates_obs_B, type = 'response', s= lasso_p_A$lambda.min)

# Estimate Logit-Lasso in Sample B
lasso_p_B <- cv.glmnet(covariates_obs_B, treat_obs_B, alpha=p, type.measure = 'mse', family="binomial")
plot(lasso_p_B)

# Extrapolate the fitted values to Sample A
pscore_A <- predict(lasso_p_B, newx = covariates_obs_A, type = 'response', s= lasso_p_B$lambda.min)

# Merge the fitted values of both samples
pscore <- rbind(pscore_A,pscore_B)

###############################################################################

###############################################################################
################################### ATE Score #################################
###############################################################################

# Merge earnings outcome of Sample A and B
earnings_obs <- rbind(earnings_obs_A,earnings_obs_B)

# Merge treatmente of Sample A and B
treat_obs <- rbind(treat_obs_A,treat_obs_B)

# Calculate the ATE score using the formula described above
Y_ate_star = invisible(???)

# Calculate ATE
# It is the sample average of the ATE score
ate <- round(mean(Y_ate_star), digits = 2)

# Calculate the standard errors of the ATE
# Square root of the quotient of variance of the ATE score and the sample size
se_ate <- round(sqrt(var(Y_ate_star)/length(Y_ate_star)), digits = 2)


print(paste0("Average Treatment Effect (ATE): ", ate))
print(paste0("Standard Error for ATE: ", se_ate))

###############################################################################

###############################################################################
################################## ATET Score #################################
###############################################################################

## Unconditional Treatment probability
p = mean(pscore)

# Calculate the ATET score using the formula described above
Y_atet_star = invisible(???)

# Calculate ATET
# It is the sample average of the ATET score
atet <- round(mean(Y_atet_star), digits = 2)

# Calculate the standard errors of the ATET
# Square root of the quotient of variance of the ATET score and the sample size
se_atet <- round(sqrt(var(Y_atet_star)/length(Y_atet_star)), digits = 2)

print(paste0("Average Treatment Effect for Treated (ATET): ", atet))
print(paste0("Standard Error for ATET: ", se_atet))

###############################################################################

###############################################################################
##################################### CATEs ###################################
###############################################################################

# Merge covariates of Sample A and B
covariates_obs <- rbind(covariates_obs_A,covariates_obs_B)

# Generate a new data frame
# Merge the ATE score and the covariates
colnames(Y_ate_star) <- "y_star"
Y_star <- as.data.frame(cbind(Y_ate_star,covariates_obs[,-c(3,8)]))

# Estimate an OLS regression
# Regress the ATE score on the covariates
cates <- lm(y_star ~., Y_star)

# Heteroskedasticity robust standard errors
coeftest(cates, vcov = vcovHC(cates, type = "HC1"))

###############################################################################

###############################################################################

# Calculate the predicted effect size for each observation
fit <- predict(cates)

# Count the observations with positive and negative effects
print(paste0("Number of individuals with positive effects: ", length(fit[fit>=0])))
print(paste0("Number of individuals with negative effects: ", length(fit[fit<0])))

###############################################################################

###############################################################################
################ Plot Cumulative Distribution of CATEs ########################
###############################################################################

plot(ecdf(fit), col="blue", xlim = c(-100,150), xlab="Effect Size (in Dollars)",
     ylab="Cumulative Distribution", main="Cumulative Distibution of the CATEs")
abline(v=0, col="red")

###############################################################################

###############################################################################
######################## Description of CATEs #################################
###############################################################################

## Means and standard deviations for individuals with positive effects
desc_1 <- fBasics::basicStats(Y_star[fit >= 0,-1]) %>% t() %>% as.data.frame() %>% select(Mean, Stdev)

## Means and standard deviations for individuals with negative effects
desc_0 <- fBasics::basicStats(Y_star[fit < 0,-1]) %>% t() %>% as.data.frame() %>% select(Mean, Stdev)

# Make table and add standardized differences
desc <- cbind(desc_1,desc_0, 
        100*abs(desc_1[,1]-desc_0[,1])/sqrt(0.5*(desc_1[,2]^2+desc_0[,2]^2)))
colnames(desc) <- c("Mean (Pos.)", "Std.Dev. (Pos.)", "Mean (Neg.)", "Std.Dev. (Neg.)", "Std.Diff.")
print(round(desc, digits=2))

###############################################################################

###############################################################################
########### Conditional Potential Earnings under Non-Participation ############
###############################################################################

# Set starting value
set.seed(123456789)

# Tuning parameters for forest
trees = 1000 # number of trees in the forest
frac = 0.5 # share of subsample used for each tree
cov = floor(1/2*ncol(covariates_obs)) # number of covariates used for each tree
min = 10 # minimum sample size in the terminal leaves of the trees

# Estimate Random Forest among non-participants in Sample A
forest_y0_A <- regression_forest(covariates_obs_A[treat_obs_A==0,], earnings_obs_A[treat_obs_A==0,],
                              num.trees = trees, sample.fraction = frac, mtry = cov, min.node.size = min)

# Extrapolate the fitted values to Sample B
y0hat_B <- as.matrix(predict(forest_y0_A, newdata = covariates_obs_B)$predictions)

print("Random Forest for Sample A estimated.")

#################################################################################

#################################################################################

# Plot one tree from the random forest
plot(tree <- get_tree(forest_y0_A, 1))
# the last number is the tree number
# it can be varied from 1 to 1000

#################################################################################

#################################################################################

# Count the splitting frequencies for each covariate
split <- split_frequencies(forest_y0_A, max.depth = 4)
# max.depth specifies the maximum tree depth we consider

# Label the results
colnames(split) <- colnames(covariates_obs)
rownames(split) <- c("Depth 1", "Depth 2", "Depth 3", "Depth 4")

print(t(split))

#################################################################################

#################################################################################

# Estimate Random Forest among non-participants in Sample B
forest_y0_B <- regression_forest(???)

# Extrapolate the fitted values to Sample A
y0hat_A <- as.matrix(predict(forest_y0_B, newdata = covariates_obs_A)$predictions)

# Merge fitted values of both samples
y0hat <- rbind(y0hat_A,y0hat_B)

print("Random Forest for Sample B estimated.")

#################################################################################

###############################################################################
########################### Propensity Score ##################################
###############################################################################

# Set starting value
set.seed(123456789)

# Tuning parameters for forest
trees = 1000
frac = 0.5
cov = floor(1/2*ncol(covariates_obs))
min = 10

# Estimate Random Forest in Sample A
forest_p_A <- regression_forest(covariates_obs_A, treat_obs_A, 
                         num.trees = trees, sample.fraction = frac, mtry = cov, min.node.size = min)

# Extrapolate the fitted values to Sample B
pscore_B <- as.matrix(predict(forest_p_A, newdata = covariates_obs_B)$predictions)

##############

# Estimate Random Forest in Sample B
forest_p_B <- regression_forest(covariates_obs_B, treat_obs_B, 
                         num.trees = trees, sample.fraction = frac, mtry = cov, min.node.size = min)

# Extrapolate the fitted values to Sample A
pscore_A <- as.matrix(predict(forest_p_B, newdata = covariates_obs_A)$predictions)

# Merge the fitted values of both samples
pscore <- rbind(pscore_A,pscore_B)

print("Propensity score is estimated.")

###############################################################################

###############################################################################
################################## ATET Score #################################
###############################################################################

## Unconditional Treatment probability
p = mean(pscore)

# Calculate the ATET score using the formula described above
Y_atet_star = invisible(treat_obs*(earnings_obs - y0hat)/p 
            - (1-treat_obs)*pscore*(earnings_obs - y0hat)/(p*(1-pscore)))

# Calculate ATET
# It is the sample average of the ATET score
atet <- round(mean(Y_atet_star), digits = 2)

# Calculate the standard errors of the ATET
# Square root of the quotient of variance of the ATET score and the sample size
se_atet <- round(sqrt(var(Y_atet_star)/length(Y_atet_star)), digits = 2)

print(paste0("Average Treatment Effect for Treated (ATET): ", atet))
print(paste0("Standard Error for ATET: ", se_atet))

###############################################################################
