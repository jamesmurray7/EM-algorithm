#' ##############
#' Implementation of the Expectation Conditional Maximisation algorithm
#' on a linear mixed effects model with random intercept only.
#' Toy example
#' Y = XB + ui + eps
#' ##############

rm(list = ls())

# Load packages and functions housed in other scripts ---------------------
library(lme4)
source("~/Documents/PhD/EM-Algorithm/fns/simlong.R")
source("~/Documents/PhD/EM-Algorithm/fns/Classic.R")

# Simulate some data ------------------------------------------------------
long.data <- simlong(num_subj = 100, num_times = 5)
X <- long.data$X; Y <- long.data$Y; Z <- long.data$Z
init.beta <- beta <- matrix(long.data$lmer.fit@beta, ncol = 1)

# Classical EM - Single run -----------------------------------------------
# Doing this as a single iteration first, easy enough to functionise once all working.
n <- nrow(Y)
q1 <- ncol(Z)
var.0 <- var.1 <- 1
l0 <- ll(X,Y,Z,init.beta, var.0, var.1)

# E step ----
# Covariance stuff
getV <- function(var.0,var.1) c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
V <- getV(var.0,var.1)
Vinv <- solve(V)
V.sqrt <- chol(V)
# X %*% beta and resulting residuals
XtX <- crossprod(X)
Xb <- X %*% beta
resid <- Y-Xb
z <- forwardsolve(t(V.sqrt), Z)
# Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)

u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
          tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
           tr(c(var.1) * diag(q1) - c(var.1)^2 * crossprod(z))

# M-step ----
# CM step 1 ---
# Estimate \sigma^(t+1)
var.0.new <- u.0/n
var.1.new <- u.1/q1
message("CM step 1 done, new sigma = [", sqrt(var.0.new), ", ", sqrt(var.1.new),"]")
# CM step 2 ---
V <- getV(var.0.new, var.1.new)
Vinv <- solve(V)
Gamma <- X %*% beta + c(var.0.new) * Vinv %*% resid
# And estimate \beta
beta.new <- solve(XtX) %*% t(X) %*% Gamma
message("CM step 2 done, new beta = ", beta.new)

ll.new <- ll(X,Y,Z,beta.new,var.0.new,var.1.new)
var.0 <- var.0.new; var.1 <- var.1.new; beta <- beta.new



# Functionise -------------------------------------------------------------
# First, simulate some data
rm(list = ls())
source("~/Documents/PhD/EM-Algorithm/fns/simlong.R")
source("~/Documents/PhD/EM-Algorithm/fns/Classic.R")
long.data <- simlong(num_subj = 100, num_times = 5)
X <- long.data$X; Y <- long.data$Y; Z <- long.data$Z

# Computational parameters: 
diff <- 10; tol <- 1e-10
iter <- 1; maxiter <- 2e3

# Initial conditions
var.1 <- var.0 <- 0.25
init.beta <- beta <- matrix(long.data$lmer.fit@beta, ncol = 1)

# Function to compute new covariance matrix, V.
getV <- function(var.0,var.1) c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)


while(diff > tol & iter <= maxiter){
  n <- nrow(Y)
  q1 <- ncol(Z)
  l0 <- ll(X,Y,Z, beta, var.0, var.1)
  
  # E step ----
  # Covariance stuff
  
  V <- getV(var.0, var.1)
  Vinv <- solve(V)
  V.sqrt <- chol(V)
  # X %*% beta and resulting residuals
  XtX <- crossprod(X)
  Xb <- X %*% beta
  resid <- Y-Xb
  
  # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
  z <- forwardsolve(t(V.sqrt), Z)
  u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
    tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
  u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
    tr(c(var.1) * diag(q1) - c(var.1)^2 * crossprod(z))
  
  # M-step ----
  # CM step 1 ---
  # Estimate \sigma^(t+1)
  var.0.new <- u.0/n
  var.1.new <- u.1/q1
  message("CM step 1 done, new sigma = [", sqrt(var.0.new), ", ", sqrt(var.1.new),"]")
  # CM step 2 ---
  V <- getV(var.0.new, var.1.new)
  Vinv <- solve(V)
  Gamma <- X %*% beta + c(var.0.new) * Vinv %*% resid
  # And estimate \beta
  beta.new <- solve(XtX) %*% t(X) %*% Gamma
  message("CM step 2 done, new beta = ", beta.new)
  
  l1 <- ll(X,Y,Z, beta.new, var.0.new, var.1.new)
  
  # New ll
  diff <- abs(l1-l0)
  # Print out
  L0.round <- round(l0, 5); L1.round <- round(l1,5); diff.round <- round(diff, 11)
  beta.new.round <- round(beta.new,2); b1 <- beta.new.round[1]; b2 <- beta.new.round[2]; b3 <- beta.new.round[3]
  sigma0.print <- round(sqrt(var.0.new), 2); sigma1.print <- round(sqrt(var.1.new),2)
  
  message("Iteration ", iter, " diff = ", diff.round)
  message("For beta = [", b1, ", ", b2, ", ", b3, "], sigma01 = [", sigma0.print, ", ", sigma1.print,"]")
  # Set new parameters
  var.0 <- var.0.new
  var.1 <- var.1.new
  beta <- beta.new
  iter <- iter + 1
}

# Remember true values should be beta = [40, -10, 1.5]; sigma01=[1.5, 0.5]