
# Set the working directory where the data set is located
setwd("C:\\Users\\ranar\\OneDrive\\Desktop\\R")

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read the gene_data.csv file into R
gene_data <- read.csv("gene_data.csv")

#### Task 2
# Extract columns from the data set
time <- gene_data[, 1]  # extracting first column and assigning to time 
x1 <- gene_data[, 2]    # extracting second column and assigning to x1 
x2 <- gene_data[, 3]    # extracting third column and assigning to x2
x3 <- gene_data[, 4]    # extracting fourth column and assigning to x3
x4 <- gene_data[, 5]    # extracting fifth column and assigning to x4
x5 <- gene_data[, 6]    # extracting sixth column and assigning to x5

## Task 2.1
# Model 1: Define matrix X1 for Model 1 (x4, x3^2, bias) (i.e.y=θ_1 x_4+θ_2 x_3^2+θ_bias)
X1 <- cbind(x4, I(x3^2), rep(1, length(time)))
theta1 <- solve(t(X1) %*% X1) %*% t(X1) %*% x2
cat("Model 1 parameters (θ1, θ2, θ_bias): ", theta1, "\n\n")

# Model 2: Define matrix X2 for Model 2 (x4, x3^2, x5, bias) (i.e. y=θ_1 x_4+θ_2 x_3^2+θ_3 x_5+θ_bias)
X2 <- cbind(x4, I(x3^2), x5, rep(1, length(time)))
theta2 <- solve(t(X2) %*% X2) %*% t(X2) %*% x2
cat("Model 2 parameters (θ1, θ2, θ3, θ_bias): ", theta2, "\n\n")

# Model 3: Define matrix X3 for Model 3 (x3, x4, x5^3) (i.e.y=θ_1 x_3+θ_2 x_4+θ_3 x_5^3)
X3 <- cbind(x3, x4, I(x5^3))
theta3 <- solve(t(X3) %*% X3) %*% t(X3) %*% x2
cat("Model 3 parameters (θ1, θ2, θ3): ", theta3, "\n\n")

# Model 4: Define matrix X4 for Model 4 (x4, x3^2, x5^3, bias) (i.e.y=θ_1 x_4+θ_2 x_3^2+θ_3 x_5^3+θ_bias)
X4 <- cbind(x4, I(x3^2), I(x5^3), rep(1, length(time)))
theta4 <- solve(t(X4) %*% X4) %*% t(X4) %*% x2
cat("Model 4 parameters (θ1, θ2, θ3, θ_bias): ", theta4, "\n\n")

# Model 5: Define matrix X5 for Model 5 (x4, x1^2, x3^2, bias) (i.e.y=θ_1 x_4+θ_2 x_1^2+θ_3 x_3^2+θ_bias)
X5 <- cbind(x4, I(x1^2), I(x3^2), rep(1, length(time)))
theta5 <- solve(t(X5) %*% X5) %*% t(X5) %*% x2
cat("Model 5 parameters (θ1, θ2, θ3, θ_bias): ", theta5, "\n\n")



################################################################################


## Task 2.2

# Calculating y hat for each Model or Performing predictions for each model
y_hat1 <- X1 %*% theta1  # Prediction for Model 1
y_hat2 <- X2 %*% theta2  # Prediction for Model 2
y_hat3 <- X3 %*% theta3  # Prediction for Model 3
y_hat4 <- X4 %*% theta4  # Prediction for Model 4
y_hat5 <- X5 %*% theta5  # Prediction for Model 5

# Calculates residuals for each model
residuals1 <- x2 - y_hat1  # Residuals for Model 1
residuals2 <- x2 - y_hat2  # Residuals for Model 2
residuals3 <- x2 - y_hat3  # Residuals for Model 3
residuals4 <- x2 - y_hat4  # Residuals for Model 4
residuals5 <- x2 - y_hat5  # Residuals for Model 5

# Calculates RSS for each model
RSS1 <- sum(residuals1^2)  # RSS for Model 1
RSS2 <- sum(residuals2^2)  # RSS for Model 2
RSS3 <- sum(residuals3^2)  # RSS for Model 3
RSS4 <- sum(residuals4^2)  # RSS for Model 4
RSS5 <- sum(residuals5^2)  # RSS for Model 5


# Printing Outputs RSS for each model
cat("Model 1 RSS:", RSS1, "\n")
cat("Model 2 RSS:", RSS2, "\n")
cat("Model 3 RSS:", RSS3, "\n")
cat("Model 4 RSS:", RSS4, "\n")
cat("Model 5 RSS:", RSS5, "\n")


## Task 2.3
# Calculate the number of observations (n)
n <- nrow(gene_data)

# Calculate the degrees of freedom (df) for each model
df1 <- ncol(X1) - 1
df2 <- ncol(X2) - 1
df3 <- ncol(X3) - 1
df4 <- ncol(X4) - 1
df5 <- ncol(X5) - 1

# Calculate the error variance (sigma squared) for each model
sigmasq_1 <- RSS1 / df1
sigmasq_2 <- RSS2 / df2
sigmasq_3 <- RSS3 / df3
sigmasq_4 <- RSS4 / df4
sigmasq_5 <- RSS5 / df5

# Compute the log-likelihood for each model
log_likelihood_1 <- -n/2 * log(2 * pi) - n/2 * log(sigmasq_1) - 1/(2 * sigmasq_1) * RSS1
log_likelihood_2 <- -n/2 * log(2 * pi) - n/2 * log(sigmasq_2) - 1/(2 * sigmasq_2) * RSS2
log_likelihood_3 <- -n/2 * log(2 * pi) - n/2 * log(sigmasq_3) - 1/(2 * sigmasq_3) * RSS3
log_likelihood_4 <- -n/2 * log(2 * pi) - n/2 * log(sigmasq_4) - 1/(2 * sigmasq_4) * RSS4
log_likelihood_5 <- -n/2 * log(2 * pi) - n/2 * log(sigmasq_5) - 1/(2 * sigmasq_5) * RSS5

# Print out the log-likelihood values for each model
cat("Model 1 log-likelihood:", log_likelihood_1, "\n")
cat("Model 2 log-likelihood:", log_likelihood_2, "\n")
cat("Model 3 log-likelihood:", log_likelihood_3, "\n")
cat("Model 4 log-likelihood:", log_likelihood_4, "\n")
cat("Model 5 log-likelihood:", log_likelihood_5, "\n")

#################
## Task 2.4
# Define degrees of freedom (K) for each model
K1 <- ncol(X1)  # Number of parameters for Model 1 (θ1, θ2, θ_bias)
K2 <- ncol(X2)  # Number of parameters for Model 2 (θ1, θ2, θ3, θ_bias)
K3 <- ncol(X3)  # Number of parameters for Model 3 (θ1, θ2, θ3)
K4 <- ncol(X4)  # Number of parameters for Model 4 (θ1, θ2, θ3, θ_bias)
K5 <- ncol(X5) # Number of parameters for Model 5 (θ1, θ2, θ3, θ_bias)

# Compute AIC for each model
AIC1 <- 2 * K1 - 2 * log_likelihood_1
AIC2 <- 2 * K2 - 2 * log_likelihood_2
AIC3 <- 2 * K3 - 2 * log_likelihood_3
AIC4 <- 2 * K4 - 2 * log_likelihood_4
AIC5 <- 2 * K5 - 2 * log_likelihood_5

# Compute BIC for each model
BIC1 <- K1 * log(n) - 2 * log_likelihood_1
BIC2 <- K2 * log(n) - 2 * log_likelihood_2
BIC3 <- K3 * log(n) - 2 * log_likelihood_3
BIC4 <- K4 * log(n) - 2 * log_likelihood_4
BIC5 <- K5 * log(n) - 2 * log_likelihood_5

# Print AIC and BIC values for each model on the same line
cat("AIC for Model 1:", AIC1, "  BIC for Model 1:", BIC1, "\n")
cat("AIC for Model 2:", AIC2, "  BIC for Model 2:", BIC2, "\n")
cat("AIC for Model 3:", AIC3, "  BIC for Model 3:", BIC3, "\n")
cat("AIC for Model 4:", AIC4, "  BIC for Model 4:", BIC4, "\n")
cat("AIC for Model 5:", AIC5, "  BIC for Model 5:", BIC5, "\n")


#################



##Task 2.5 
# Setting up a 2x3 grid for plotting
par(mfrow = c(2, 3))
# Plotting histograms of residuals for each model
hist(residuals1, main = "Model 1 Residuals", xlab = "Residuals")
hist(residuals2, main = "Model 2 Residuals", xlab = "Residuals")
hist(residuals3, main = "Model 3 Residuals", xlab = "Residuals")
hist(residuals4, main = "Model 4 Residuals", xlab = "Residuals")
hist(residuals5, main = "Model 5 Residuals", xlab = "Residuals")

# Plotting Q-Q plots of residuals for each model
qqnorm(residuals1, main = "Q-Q Plot: Model 1 Residuals")
qqline(residuals1, col = 2)

qqnorm(residuals2, main = "Q-Q Plot: Model 2 Residuals")
qqline(residuals2, col = 2)

qqnorm(residuals3, main = "Q-Q Plot: Model 3 Residuals")
qqline(residuals3, col = 2)

qqnorm(residuals4, main = "Q-Q Plot: Model 4 Residuals")
qqline(residuals4, col = 2)

qqnorm(residuals5, main = "Q-Q Plot: Model 5 Residuals")
qqline(residuals5, col = 2)


##############################



## Task 2.7
# Set seed for reproducibility
set.seed(123)

# Sample 70% of rows for training data indices
train_indices <- sample(1:nrow(gene_data), 0.7 * nrow(gene_data))

# Create training and testing datasets
train_data <- gene_data[train_indices, ]
test_data <- gene_data[-train_indices, ]

# Extract columns for training data
train_x1 <- train_data[, 2]
train_x3 <- train_data[, 4]
train_x4 <- train_data[, 5]
train_time <- train_data[, 1]
train_y <- train_data[, 3]

# Extract columns for testing data
test_x1 <- test_data[, 2]
test_x3 <- test_data[, 4]
test_x4 <- test_data[, 5]
test_time <- test_data[, 1]
test_y <- test_data[, 3]

# Construct design matrix for training with Model 5 (y = θ1*x4 + θ2*x1^2 + θ3*x3^2 + θ_bias)
X_train <- cbind(train_x4, I(train_x1^2), I(train_x3^2), rep(1, length(train_time)))

# Compute parameters θ using the normal equation
theta <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% train_y

# Construct design matrix for testing using the same predictors as in X_train
X_test <- cbind(test_x4, I(test_x1^2), I(test_x3^2), rep(1, length(test_time)))

# Predict y using the trained coefficients θ
y_pred <- X_test %*% theta

# Calculate residuals (difference between actual and predicted values)
residuals <- test_y - y_pred

# Calculate residual variance
residual_var <- sum(residuals^2) / (length(test_y) - length(theta))

# Calculate standard errors
se <- sqrt(diag(X_test %*% solve(t(X_train) %*% X_train) %*% t(X_test)) * residual_var)

# Calculate 95% confidence intervals for predictions
ci_upper <- y_pred + qt(0.975, df = length(test_y) - length(theta)) * se
ci_lower <- y_pred - qt(0.975, df = length(test_y) - length(theta)) * se


# Plotting results
plot(test_time, test_y, pch = 16, col = "blue", xlab = "Time", ylab = "Gene Expression (x2)", 
     main = "Model 5 Predictions with 95% CI")
lines(test_time, y_pred, col = "red", lwd = 2)  # Plot predicted values
segments(test_time, ci_lower, test_time, ci_upper, col = "green", lwd = 2)  # Plot 95% CI
legend("topright", legend = c("Actual", "Predicted", "95% CI"), 
       col = c("blue", "red", "green"), lwd = c(1, 1, 2), pch = c(16, NA, NA))

################


#Task 3.
# Parameters from Model 5
theta1 <- -0.8312983
theta2 <- 0.5385828
theta3 <- -0.1096679
theta_bias <- 1.295152
RSS5 <- 8.516473  # observed RSS from Model 5

# Define prior distributions
theta1_prior <- runif(10000, min = -2, max = 0)
theta2_prior <- runif(10000, min = 0, max = 1)

# Initialize storage for accepted samples
accepted_theta1 <- numeric()
accepted_theta2 <- numeric()

# Perform rejection ABC
for (i in 1:length(theta1_prior)) {
  # Generate RSS for current sample
  y_hat <- theta_bias + theta1_prior[i] * x4 + theta2_prior[i] * (x1^2) + theta3 * (x3^2)
  RSS <- sum((x2 - y_hat)^2)
  
  # Acceptance criterion
  if (RSS < RSS5 * 1.2) {  
    accepted_theta1 <- c(accepted_theta1, theta1_prior[i])
    accepted_theta2 <- c(accepted_theta2, theta2_prior[i])
  }
}

# Check if any samples were accepted
if (length(accepted_theta1) == 0) {
  cat("No samples accepted based on the rejection criterion.")
} else {
  # Plotting posterior distributions
  
  # Joint posterior plot
  hist(accepted_theta1, breaks=30, main="Posterior distribution of θ1", xlab="θ1")
  hist(accepted_theta2, breaks=30, main="Posterior distribution of θ2", xlab="θ2")
  
  plot(accepted_theta1, accepted_theta2, main="Joint posterior distribution", xlab="θ1", ylab="θ2")
  
  # Marginal posterior plot for θ1
  hist(accepted_theta1, breaks=30, main="Marginal posterior distribution of θ1", xlab="θ1")
  
  # Marginal posterior plot for θ2
  hist(accepted_theta2, breaks=30, main="Marginal posterior distribution of θ2", xlab="θ2")
}
