########################  Load Packages and Data  ########################

# Load packages
library(rpart)
library(rpart.plot)
library(grf)
library(glmnet)

# Load data
juice <- read.csv("juice.csv", sep = ",")
new_grocery <- read.csv("new_grocery.csv", sep = ",")

print('Packages and data successfully loaded.')

#############################################################################

########################  Describe Old Data  ########################

# Print first few rows of old data
head(juice)

# Number of observations
print(paste0('Old data: ',nrow(juice),' observations'))

######################################################################

########################  Describe Old Data  ########################

# Print first few rows of new data
head(new_grocery)

# Number of observations
print(paste0('New data: ',nrow(new_grocery),' observations'))

######################################################################

########################  Data Preparation  ########################

# Generate dummy for missing prices
missing <- (is.na(juice$price) == TRUE)
new_missing <- (is.na(new_grocery$price) == TRUE)

# Replace missing prices with zero
juice$price[is.na(juice$price)] <-0
new_grocery$price[is.na(new_grocery$price)] <-0

# Generate Dummies for Brands
brand_1 <- (juice$brand == "minute.maid")
brand_2 <- (juice$brand == "dominicks")
brand_3 <- (juice$brand == "tropicana")

new_brand_1 <- (new_grocery$brand == "minute.maid")
new_brand_2 <- (new_grocery$brand == "dominicks")
new_brand_3 <- (new_grocery$brand == "tropicana")

# Generate outcome and control variables
y <- as.matrix(juice$sales)
colnames(y) <- c("sales")

x <- as.matrix(cbind(juice$price, missing, brand_1, brand_2, brand_3, juice$feat))
colnames(x) <- c("price", "missing", "minute.maid", "dominicks", "tropicana", "featured")

new_x <- as.matrix(cbind(new_grocery$price, new_missing, new_brand_1, new_brand_2, new_brand_3, new_grocery$feat))
colnames(new_x) <- c("price", "missing", "minute.maid", "dominicks", "tropicana", "featured")

# Descriptive statistics
summary(cbind(y,x))

print('Data is prepared.')

#############################################################################

########################  Training and Test Samples  ########################

set.seed(???)

# Generate variable with the rows in training data


print('Training and test samples created.')

#############################################################################

########################  LASSO, Ridge, Elastic Net  ##############################

set.seed(???)
penalized.cv <- ???


# Fitted values
pred_penalized <- ???

# Calculate the MSE
MSE_penalized <- mean((y[-training_set] - pred_penalized[-training_set])^2)
R2_penalized <- round(1- MSE_penalized/var(y[-training_set]), digits = 3)

print(paste0("R-squared Penalized Regression: ", R2_penalized))
                                   
################################################################

######################  Regression Tree  #######################

set.seed(???)
# Prepare data for tree estimator
outcome <- y[training_set]
tree_data <- data.frame(outcome, x[training_set,])

deep_tree <- ???

# Optimal tree size
op.index <- ???

## Select the Tree that Minimises CV-MSE
cp.vals <- ???

# Prune the deep tree
pruned_tree <- ???

## Plot tree structure
#rpart.plot(pruned_tree,digits=3)

# Fitted values
predtree <- ???

# Calculate the MSE
MSEtree <- mean((y[-training_set] - predtree[-training_set])^2)
R2tree <- round(1- MSEtree/var(y[-training_set]), digits = 3)

print(paste0("R-squared Tree: ", R2tree))

################################################################

########################  Random Forest  #######################

set.seed(???)

rep <- ??? # number of trees
cov <- ??? # share of covariates
frac <- ??? # fraction of subsample
min_obs <- ??? # max. size of terminal leaves in trees

# Build Forest
forest <- ???

# Fitted values
predforest <- ???

# Calculate MSE
MSEforest <- mean((y[-training_set] - predforest[-training_set])^2)
R2forest <- round(1- MSEforest/var(y[-training_set]), digits = 3)

print(paste0("R-squared Forest: ", R2forest))

################################################################

########################  Out-of-Sample Prediction  #######################

# Fitted values
new_prediction <- ???

print('Out-of-sample sales are predicted.')

###########################################################################

########################  Store Results  #######################

id_new <- as.matrix(new_grocery$id)

# Replace ??? with your group name
write.csv(cbind(id_new,new_prediction),"???.csv")

print('File is stored.')
print('Send your results to anthony.strittmatter@unibas.ch')

################################################################


