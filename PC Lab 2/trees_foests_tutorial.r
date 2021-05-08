########################  Load Packages and Data  ########################

# Load packages
library(rpart)
library(rpart.plot)
library(grf)
library(DiagrammeR)

# Load data
data_2006 <-read.csv("browser_2006.csv", sep = ",")
data_new <-read.csv("browser_new.csv", sep = ",")

# Data preparation
y_2006 <- as.matrix(data_2006[,2])
x_2006 <- as.matrix(data_2006[,c(3:ncol(data_2006))])
id_2006 <- as.matrix(data_2006[,1])
x_new <- as.matrix(data_new[,c(2:ncol(data_new))])
id_new <- as.matrix(data_new[,1])

print('Packages and data successfully loaded.')

#############################################################################

########################  Average Spending  ########################

spending <- round(???, digits=2)
print(paste0("In 2006, the average spending is ", spending, " US-dollars"))

####################################################################

########################  Online Time  ########################

freq <- round(x_2006[id_2006==921,x_2006[id_2006==921,] == ???], digit = 0)
page <- names(freq)

print(paste0("Household 921 is most of the time on the webpage ", page))
print(paste0(freq, "% of the online time is the household on this webpage"))

################################################################

########################  Log Transformation  ########################

log_y_2006 = as.matrix(???) # take logarithm

# Cumulative Distribution of Spending
plot(ecdf(y_2006), xlab = "Spending in US-Dollars", sub = "(Truncated at 20,000 US-Dollars)",
     ylab = "cdf", main = "Distribution of Spending", xlim= c(0,20000))

# Cumulative Distribution of Log Spendiung
plot(ecdf(log_y_2006), xlab = "log Spending", ylab = "cdf", main = "Distribution of Log Spending")

#######################################################################

########################  Training and Test Samples  ########################

set.seed(1001)
# Generate variable with the rows in training data
size <- floor(0.5 * nrow(data_2006))
training_set <- sample(seq_len(nrow(data_2006)), size = size)

print('Training and test samples created.')

#############################################################################

########################  Shallow Tree  ########################

# Prepare data for tree estimator
outcome <- log_y_2006[training_set]
tree_data_2006 <-  data.frame(outcome, x_2006[training_set,])

# Build shallow tree
set.seed(1001)
shallow_tree <- rpart(formula = outcome ~., data = tree_data_2006, method = "anova", xval = 10,
                             y = TRUE, control = rpart.control(cp = 0.00002, minbucket=150))
# Note: 'minbucket=100' imposes the restriction that each terminal leave should contain at least 100 observations. 
# The algorithm 'rpart' stops growing trees when either one leave has less than 100 observations or 
# the MSE gain of addidng one addidtional leave is below cp=0.00002.

## Plot tree structure
rpart.plot(shallow_tree,digits=3)

# bizrate.com
# fedex.com

################################################################

########################  Deep Tree  ########################
set.seed(1001)
deep_tree <- rpart(formula = outcome ~., data = tree_data_2006, ???)

print('Relative CV-MSE for different tree sizes')
print(deep_tree$cptable)

# Plot CV-MSE
plotcp(deep_tree)

#############################################################

########################  Optimal Tree Size  ########################

op.index <- which.min(deep_tree$cptable[, "xerror"])
op.size <- deep_tree$cptable[op.index, "nsplit"] +1
print(paste0("Optimal number final leaves: ", op.size))

#####################################################################

########################  Pruned Tree  ########################

# Select the Tree that Minimises CV-MSE 
# Get cp-value that corresponds to optimal tree size
cp.vals <- deep_tree$cptable[op.index, "CP"]

# Prune the deep tree
pruned_tree <- prune(???, cp = cp.vals)

## Plot tree structure
rpart.plot(pruned_tree,digits=3)

# aggregateknowledge.com

################################################################

########################  Out-of-Sample Performance  ########################

# Predict log online spending 
pred_tree <- predict(???, newdata= as.data.frame(x_2006))

# Test sample data
outcome_test <- log_y_2006[-training_set]
pred_tree_test  <- pred_tree[-training_set]

# R-squared
MSE_tree <- mean((outcome_test-pred_tree_test)^2)
r2_tree <-  round(1- MSE_tree/var(outcome_test), digits = 3) 
print(paste0("Test sample R-squared: ", r2_tree))

##############################################################################

########################  Random Forest  ########################

rep <- 1000 # number of trees
cov <- 1/3 # share of covariates
frac <- 1/2 # fraction of subsample
min_obs <- 100 # max. size of terminal leaves in trees

# Build Forest
set.seed(10001)
forest1 <- regression_forest(x_2006[training_set,],log_y_2006[training_set,], 
                            mtry = floor(cov*ncol(x_2006)), sample.fraction = frac, num.trees = rep, 
                            min.node.size = min_obs, honesty=FALSE)

print('Forest is built.')

##################################################################

########################  Plot Example Tree  ########################

# Plot a tree of the forest
# Just an illustration, overall the forest contains 1000 trees
tree <- get_tree(???,1) # here we select tree number 1
plot(tree)

#####################################################################

########################  Variable Importance  ########################

# Plot the variable importantance
# First we consider only first split
imp1 <- variable_importance(forest1, max.depth = 1)
print(cbind(colnames(x_2006[,imp1>0.02]),imp1[imp1>0.02]))

# Now we consider the first four splits
imp2 <- round(variable_importance(forest1, decay.exponent = 2, max.depth = 4), digits = 3)
print(cbind(colnames(x_2006[,imp2>0.02]),imp2[imp2>0.02]))

########################################################################

########################  Out-of-Sample Performance  ########################

# Prediction
fit <- predict(???, newdata = x_2006[-training_set,])$predictions

# R-squared
SST <- mean(((log_y_2006[-training_set,])-mean((log_y_2006[-training_set,])))^2)
MSE1 <- mean(((log_y_2006[-training_set,])-fit)^2)
r2_1 <-  round(1- MSE1/SST, digits = 3) 
print(paste0("Test sample R-squared: ", r2_1))

#############################################################################

########################  Area Under the Curve (AUC)  ########################

sizes <- c(1000,500,400,300, 200, 100, 50, 40,30,20,10, 5,4,3,2,1) # Select a grid of sample sizes
# Prepare matrix to store results
auc <- matrix(NA, nrow = length(sizes), ncol = 3)
colnames(auc) <- c("Trees", "AUC", "Marginal AUC")
auc[,1] <- sizes
# Sum of Squares Total
SST <- mean(((log_y_2006[-training_set,])-(mean(log_y_2006[-training_set,])))^2)

set.seed(10001) # set starting value
for (t in sizes){
  # Estimate Forests
  forest <- regression_forest(x_2006[training_set,],(log_y_2006[training_set,]), mtry = floor(cov*ncol(x_2006)),
                              sample.fraction = frac, num.trees = t, min.node.size = min_obs, honesty=FALSE)
  fit <- predict(forest, newdata = x_2006[-training_set,])$predictions # prediction in test sample
  auc[auc[,1]== t,2] <- 1- mean(((log_y_2006[-training_set,])-fit)^2)/SST # store R-squared
}
auc[,3] <- auc[,2] - rbind(as.matrix(auc[-1,2]),auc[nrow(auc),2])

# Marginal AUC
plot(auc[,1],auc[,2],type = "o",xlab="Trees", ylab= "R-squared", main = "AUC")
abline(a=0,b=0, col="red")

################################################################################

########################  Deep Forest  ########################

min_obs <- 5
# Build Forest
forest2 <- regression_forest(x_2006[training_set,],log_y_2006[training_set,], 
                            ???)

# Prediction
fit <- predict(forest2, newdata = x_2006[-training_set,])$predictions

# R-squared
SST <- mean((log_y_2006[-training_set,]-mean(log_y_2006[-training_set,]))^2)
MSE2 <- mean((log_y_2006[-training_set,]-fit)^2)
r2_2 <-  round(1- MSE2/SST, digits = 3)
print(cbind(r2_1,r2_2))

# Plot tree
tree <- get_tree(forest2, 34)
plot(tree)

###############################################################

########################  Store Prediction for Hold-out-Sample  ########################

# Hold-out-Sample Prediction
fit_new <- predict(???, newdata = x_new)$predictions

results <- as.matrix(cbind(id_new,fit_new)) # store ID's and predictions in oine matrix
colnames(results) <- c("id","predictions") # label columns

# Store results
write.csv(results, "predictions.csv")

print('Results for the hold-out-sample stored.')

#########################################################################################
