########################  Load Packages and Data  ########################

# Load packages
library(glmnet)
library(corrplot)

# Load data
load("student-mat-train.Rdata")
load("student-mat-test.Rdata")

# Number of observations
print(paste0('Training set: ',nrow(train),' obs'))
print(paste0('Test set: ',nrow(test),' obs'))

###########################################################################

########################  Correlation analysis  ########################
cor <- round(cor(train[,c(1:25)]),2) # Variable 26 is the depedendent variable
corrplot(cor)

########################  Estimation of the linear regression  ########################

ols <- lm(G3 ~ ., data = train)
summary(ols)

# Calculate the MSE
test$predols <- predict(ols, newdata = test)

predMSEols <- mean((test$G3 - test$predols)^2)
print(predMSEols)

########################################################################################

########################  OLS model  ########################

ols_small <- lm(??? , data = train)

# Calculate the MSE
test$predols_small <- predict(ols_small, newdata = test)

predMSEols_small <- mean((test$G3 - test$predols_small)^2)
print(predMSEols_small)

########################  Lasso Path  ########################

# We make a plot that shows how the Lasso coefficients change with lambda
# glmnet is the standard R package for Lasso, Ridge, and Elastic Net
# alpha is a parmeter that allows to specify a Lasso, Ridge, or Elastic Net model
# alpha = 1 for Lasso; alpha = 0 for Ridge, 0 < alpha < 1 for Elastic Net
# The control variables are train[,c(1:25)]
# The outcome variable is train$G3 (math grades)

# Estimate a Lasso model
lasso <- glmnet(as.matrix(train[,c(1:25)]), train$G3, alpha = 1) # We save the model under the name "lasso"
plot(lasso, xvar = "lambda", label = TRUE)

###############################################################

########################  Cross-Validaton  ########################

# Set starting value for replicability
set.seed(27112019) 

# cv.glmnet performs a cross-validation to determine the optimal lambda value
# type.measure specifies the measure we use to assess the model accuracy (here MSE)
# nfolds specifies the number of cross-validation folds we use (here 5)

# Cross-validate the Lasso
lasso.cv <- cv.glmnet(as.matrix(train[,c(1:25)]), train$G3, type.measure = "mse", nfolds = 5, alpha = 1)

# Plot the MSE for the different lambda values
plot(lasso.cv)

#####################################################################

########################  Optimal Lambda Value  ########################

# Print the optimal lambda value
print(paste0("Optimal lambda that minimizes cross-validated MSE: ", lasso.cv$lambda.min))
print(paste0("Optimal lambda using one-standard-error-rule: ", lasso.cv$lambda.1se))

#########################################################################

########################  Lasso Coefficients  ########################

# Print Lasso coefficients
print(coef(lasso.cv, s = "lambda.min"))

# Save for later comparison
coef_lasso1 <- coef(lasso.cv, s = "lambda.min") 

#######################################################################

########################  Test Sample MSE  ########################

# Estimate the fitted values of the Lasso model in the test sample
# We use the model "lasso.cv" and the lambda value which we estimated in the training sample
# The control variables "newx" are from the test sample

# Fitted values
test$predlasso <- predict(lasso.cv, newx = as.matrix(test[,c(1:25)]), s = lasso.cv$lambda.min)

# Calculate the MSE
predMSElasso <- mean((test$G3 - test$predlasso)^2)
print(paste0("MSE: ", predMSElasso))
      
#####################################################################

########################  Different Starting Value  ########################

# Change the starting value
set.seed(27112025) # 27112024

# Re-estimate the Lasso model
lasso.cv <- cv.glmnet(???)

# Store the coefficients
coef_lasso2 <- coef(lasso.cv, s = ???)
print(cbind(coef_lasso1, coef_lasso2))

# Calculate the fitted values
test$predlasso2 <- predict(lasso.cv, newx = as.matrix(test[,c(1:25)]), s = lasso.cv$lambda.min)

# Correlation between the fitted values of the two Lasso models
cor_fit <- cor(test$predlasso,test$predlasso2)
print(paste0("Correlation between fitted values: ", cor_fit))

########################  Ridge Path  ########################

# alpha = 0 specifies a Ridge model

# Estimate the Ridge
ridge <- glmnet(as.matrix(train[,c(1:25)]), train$G3, alpha = ???)

# Plot the path of the Ridge coefficients
plot(ridge, xvar = "lambda", label = TRUE)

###############################################################

########################  Cross-Validation  ########################

# Set starting value
set.seed(27112019)

# Cross-validate the Ridge model 
ridge.cv <- cv.glmnet(???)

# Plot the MSE in the cross-validation samples
plot(ridge.cv)

#####################################################################

########################  Optimal Lambda Value  ########################

# Print the optimal lambda value
print(paste0("Optimal lambda that minimizes cross-validated MSE: ", ???))
print(paste0("Optimal lambda using one-standard-error-rule: ", ???))

#########################################################################

########################  Ridge Coefficients  ########################

# Print Ridge coefficients
print(coef(ridge.cv, s = "lambda.min"))

# Save for later comparison
coef_ridge <- coef(ridge.cv, s = "lambda.min") 

#######################################################################

########################  Test Sample MSE  ########################

# Estimate fitted values in test sample
test$predridge <- predict(ridge, newx = ???, s = ???)

# Calculate the MSE
predMSEridge <- ???
print(paste0("MSE: ", predMSEridge))

###################################################################

########################  Compare Lasso and Ridge Coefficients  ########################

# Pick the coefficients of Dalc and Walc
comp <- cbind(coef(ols)[23:24], coef_lasso1[23:24], coef_lasso2[23:24], coef_ridge[23:24]) 
colnames(comp) <- c("OLS", "Lasso1", "Lasso2", "Ridge")
print(comp)

#########################################################################################

########################  Compare the MSE  ########################

# Print the MSE of the OLS, Lasso and Ridge models
print(c(predMSEols, predMSElasso, predMSEridge))

####################################################################

########################  Compare models  ########################

# Visualize the predictions (Predicted vs Actual)
plot(test$G3,test$predols,xlim=c(5,20),ylim=c(4,16), col= "darkgreen", xlab = "Actual Grades", ylab = "Predicted Grades" )
par(new=TRUE)
plot(test$G3,test$predlasso,xlim=c(5,20),ylim=c(4,16), col= "blue", xlab = "", ylab = "" )
par(new=TRUE)
plot(test$G3,test$predridge,xlim=c(5,20),ylim=c(4,16), col= "red", xlab = "", ylab = "" )
abline(a=0,b=1)
legend(16, 9, c("OLS", "Lasso", "Ridge"), col = c("darkgreen", "blue", "red"), pch = c(21, 21, 21))

####################################################################
