#' ##############
#' Implementation of the 'classical' EM algorithm
#' on a linear mixed effects model with random intercept only.
#' Toy example
#' Y = XB + ui + eps
#' ##############

rm(list = ls())
library(lme4)

# Simulate some data ------------------------------------------------------

# Parameters
sigma.0 <- 1.5
sigma.1 <- 0.5
num_subj <- 100; num_times <- 6
N <- num_subj * num_times
beta <- matrix(c(40, -10, 1.5), ncol = 1)
# Covariates 
x1 <- rbinom(num_subj, 1, 0.5)
x2 <- floor(runif(num_subj, 50, 75))
X <- model.matrix(~rep(x1, each = num_times)+rep(x2, each = num_times))
# Random effects and error term
u <- rnorm(num_subj, 0, sigma.1)
eps <- rnorm(N, 0, sigma.0)
Y <- X %*% beta + rep(u, each = num_times) + eps

long.data <- data.frame(id = rep(1:num_subj, each = num_times),
                        x1 = rep(x1, each = num_times),
                        x2 = rep(x2, each = num_times),
                        Y = Y)
# Actual fits
actual.fit <- lmer(Y ~ x1 + x2 + (1|id), data = long.data, REML = F)
summary(actual.fit)
actual.ll <- summary(actual.fit)$logLik
Z <- as.matrix(getME(actual.fit, "Z"))

# Initial estimates -- Coefficients / parameters ----
var.1 <- var.0 <- 1
beta <- matrix(c(30, -5, 0.5),ncol=1)

# Computational parameters
tol <- 1e-10
maxiter <- 2e3

# Define functions we'll need

# Trace - sum of diagonals matrix
tr <- function(x) sum(diag(x))

# Calculating the log likelihood for LME models
ll <- function(Y,X,Z,beta,var.0,var.1){
  # Covariance matrix
  V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(length(Y))
  Vinv <- solve(V) # sqrt
  
  # X %*% beta and the residuals (epsilon)
  Xb <- X %*% beta
  resid <- Y-Xb
  
  # The log-likelihood 
  #   same as t(resid) %*% Vinv %*% resid --> sum(forwardsolve(l = t(V.sqrt), x = resid) ^ 2 )
  -0.5 * (t(resid) %*% Vinv %*% resid) - 0.5 * log(det(V)) # determinant(V, logarithm = T)
                                                           # or sum(log(diag(cV))) * 2 where cV cholesky(V)
}

ll.chol <- function(Y,X,Z,beta,var.0,var.1){ # Should be faster
  V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(length(Y))
  V.sqrt <- chol(V)
  detV <- sum(log(diag(V.sqrt))) * 2
  # X, beta and residuals
  Xb <- X %*% beta
  resid <- Y - Xb
  temp <- sum(forwardsolve(l = t(V.sqrt), x = resid) ^ 2)
  return(
    -0.5 * temp - 0.5 * detV
  )
}

# Values from actual fit to check ...
ll(Y,X,Z,beta = actual.fit@beta, var.0 = sigma(actual.fit) ^ 2, var.1 = 0.5 ^ 2) 
ll.chol(Y,X,Z,beta = actual.fit@beta, var.0 = sigma(actual.fit) ^ 2, var.1 = 0.5 ^ 2)
actual.ll #  Seems to be appx. double that produced by my ll function. 

# Classical EM - Single run -----------------------------------------------
# Doing this as a single iteration first, easy enough to functionise once all working.
n <- nrow(Y)
q1 <- nrow(Z) #  Unsure if this should be 100 (i.e. number of unique U's or 600, length of matrix)
l0 <- ll(Y,X,Z,beta, var.0, var.1)
i <- 0

# E step ----
# Covariance stuff
V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
Vinv <- solve(V)
# X %*% beta and resulting residuals
Xb <- X %*% beta
resid <- Y-Xb
# Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)

u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
          tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
           tr(c(var.1) * diag(100) - c(var.1)^2 * t(Z) %*% Vinv %*% Z) #  This not working for diag(600) makes me think Z may be wrong.

# Alternative, commented out atm
# temp1 <- Vinv %*% resid
# u.0 <- c(var.0)^2 * t(temp1) %*% temp1 + c(var.0) * n - c(var.0)^2 * tr(Vinv)
# u.1 <- c(var.1)^2 * t(temp1) %*% Z %*% t(Z) %*% temp1 + c(var.1) * q1 - c(var.1)^2 * tr(t(Z) %*% Vinv %*% Z)

# Calculate expectation y-ZU conditional on y, call this w
w <- Xb + c(var.0) * Vinv %*% resid

 # M-step ----
# Estimate \sigma^(t+1)
var.0.new <- u.0/600
var.1.new <- u.1/100

# And estimate \beta
beta.new <- solve(t(X) %*% X) %*% t(X) %*% w
message(
  "old sigma.0 = ", sqrt(var.0), " new sigma.0 = ", sqrt(var.0.new), 
  "\nold sigma.1 = ", sqrt(var.1), " new sigma.1 = ", sqrt(var.1.new),
  "\nold beta: [", beta[1], beta[2], beta[3],"] new beta: [", beta.new[1], beta.new[2], beta.new[3],"]"
)
ll.new <- ll(Y, X, Z, beta.new, var.0.new, var.1.new)
ll.new <- ll.chol(Y,X,Z,beta.new,var.0.new,var.1.new)
var.0 <- var.0.new; var.1 <- var.1.new; beta <- beta.new


# And then repeat the algorithm until diff < tol ! 
#    Choleskys / SVDs
#    Matrix inversions


# Functionise -------------------------------------------------------------
# Computational parameters: 
diff <- 10; tol <- 1e-10
iter <- 1; maxiter <- 2e3

n <- nrow(Y)
q1 <- nrow(Z) #  Unsure if this should be 100 (i.e. number of unique U's or 600, length of matrix)
# Initial conditions
var.1 <- var.0 <- 0.25
# beta <- matrix(c(30, -5, 0.5),ncol=1)
beta <- matrix(lm(Y~x1+x2, data = long.data)$coef, ncol = 1)

while(diff > tol & iter <= maxiter){
  L0 <- ll.chol(Y, X, Z,  beta, var.0, var.1) 
  
  # E step ----
  # Covariance stuff
  V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
  Vinv <- solve(V)
  # X %*% beta and resulting residuals
  Xb <- X %*% beta
  resid <- Y-Xb
  # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
  
  u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
    tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
  u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
    tr(c(var.1) * diag(100) - c(var.1)^2 * t(Z) %*% Vinv %*% Z) #  This not working for diag(600) makes me think Z may be wrong.
  
  w <- Xb + c(var.0) * Vinv %*% resid
  
  # M step ----
  var.0.new <- u.0/600
  var.1.new <- u.1/100
  beta.new <- solve(t(X) %*% X) %*% t(X) %*% w
  
  # New ll
  L1 <- ll.chol(Y,X,Z,beta.new,var.0.new, var.1.new)
  diff <- abs(L1-L0)
  # Print out
  L0.round <- round(L0, 5); L1.round <- round(L1,5); diff.round <- round(diff, 11)
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