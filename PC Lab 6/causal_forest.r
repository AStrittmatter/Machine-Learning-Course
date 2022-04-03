########################  Load Packages  ########################

# List of required packages
pkgs <- c('fBasics', 'corrplot', 'tidyverse', 'grf', 'plotmo')

# Load packages
for(pkg in pkgs){
    library(pkg, character.only = TRUE)
}
options(warn=-1)     # supress warnings

print('All packages successfully installed and loaded.')

##################################################################

########################  Load Data Frame  ########################

# Load data frame
df <- read.csv("fundraising.csv",header=TRUE, sep=",")

# Outcome Variable
outcome <- c("char_giving")

# Treatment Variables
treatment <- c("treat")

# Covariates/Features
covariates <- c("amount_pre", "amount_lastpre", "amount_maxpre", "H_number_yearbefore", "H_ngifts",
                "H_littleask", "H_bigask", "H_nyears", "H_frequency", "H_medinc", "H_medinc_mdum",
                "H_Avg_years_ed", "H_Avg_years_ed_mdum")

    
all_variables <- c(outcome, treatment, covariates)

print('Data frame successfully loaded and sample selected.')

####################################################################

########################  Table with Descriptive Statistics  ########################

desc <- fBasics::basicStats(df) %>% t() %>% as.data.frame() %>% 
  select(Mean, Stdev, Minimum, Maximum, nobs)
print(round(desc, digits=2))

#####################################################################################

########################  Correlation Matrix  ########################

corr = cor(df[,-c(1:2)])
corrplot(corr, type = "upper", tl.col = "black")

######################################################################

########################  Partition the Samples  ########################
set.seed(100239) # set starting value for random number generator

# Partition Hold-Out-Sample
df_part <- modelr::resample_partition(df, c(obs = 0.8, hold_out = 0.2))
df_obs <- as.data.frame(df_part$obs) # Training and estimation sample
df_hold_out <- as.data.frame(df_part$hold_out) # Hold-out-sample

print('Samples are partitioned.')

########################  Generate Variables  ########################

# Outcome
giving_hold_out <- as.matrix(df_hold_out[,1])
giving_obs <- as.matrix(df_obs[,1])

# Treatment
treat_hold_out <- as.matrix(df_hold_out[,2])
treat_obs <- as.matrix(df_obs[,2])

# Covariates
covariates_hold_out <- as.matrix(df_hold_out[,c(3:ncol(df_hold_out))])
covariates_obs <- as.matrix(df_obs[,c(3:ncol(df_obs))])

print('The data is now ready for your analysis!')

#######################################################################

########################  Causal Forest  ######################## 
set.seed(100244)

# Tuning parameters
min_tree = 100 # Minimum size of terminal leaves
num_trees = 1000 # Number of trees in forest
cov_frac = 1/2 # Fraction of covariates in each tree
sample_part= 0.5 # Fraction of sample used for each tree (subsampling)

# Caual Forest
cates <- causal_forest(covariates_obs, giving_obs, treat_obs,
                  sample.fraction = sample_part, mtry = floor(cov_frac*ncol(covariates_obs)), 
                  num.trees = num_trees, min.node.size = min_tree,
                  honesty = TRUE, honesty.fraction = 0.5)

print('Forest is ready!')

###################################################################

#################################################################################

# Plot one tree from the random forest
plot(tree <- get_tree(cates, 1))
# the last number is the tree number
# it can be varied from 1 to 1000

#################################################################################

#################################################################################

# Count the splitting frequencies for each covariate
split <- split_frequencies(cates, max.depth = 4)
# max.depth specifies the maximum tree depth we consider

# Label the results
colnames(split) <- colnames(covariates_obs)
rownames(split) <- c("Depth 1", "Depth 2", "Depth 3", "Depth 4")

print(t(split))

#################################################################################

######################### ATE ###############################

average_treatment_effect(cates, target.sample = c("all"))

#############################################################

###############################################################################

# Calculate the predicted effect size for each observation
fit <- predict(cates, covariates_hold_out, estimate.variance = FALSE)$predictions

# Count the observations with positive and negative effects
print(paste0("Number of individuals with positive effects: ", length(fit[fit>=0])))
print(paste0("Number of individuals with negative effects: ", length(fit[fit<0])))

print(paste0("Share of individuals with positive effects: ", round(100*length(fit[fit>=0])/length(fit),digits=1), "%"))

###############################################################################

###############################################################################
################ Plot Cumulative Distribution of CATEs ########################
###############################################################################

plot(ecdf(fit), col="blue", xlim = c(-25,25), xlab="Effect Size (in Dollars)",
     ylab="Cumulative Distribution", main="Cumulative Distibution of the CATEs")
abline(v=0, col="red")

###############################################################################

###############################################################################
######################## Description of CATEs #################################
###############################################################################

## Means and standard deviations for individuals with positive effects
desc_1 <- fBasics::basicStats(covariates_hold_out[fit >= 0,]) %>% t() %>% as.data.frame() %>% select(Mean, Stdev)

## Means and standard deviations for individuals with negative effects
desc_0 <- fBasics::basicStats(covariates_hold_out[fit < 0,]) %>% t() %>% as.data.frame() %>% select(Mean, Stdev)

# Make table and add standardized differences
desc <- cbind(desc_1,desc_0, 
        100*abs(desc_1[,1]-desc_0[,1])/sqrt(0.5*(desc_1[,2]^2+desc_0[,2]^2)))
colnames(desc) <- c("Mean (Pos.)", "Std.Dev. (Pos.)", "Mean (Neg.)", "Std.Dev. (Neg.)", "Std.Diff.")
print(round(desc, digits=2))

###############################################################################


