#' ##############
#' Implementation of the 'classial' EM algorithm
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
num_subj <- 100
num_times <- 6
N <- num_subj * num_times
beta <- matrix(c(40, -10, 1.5), ncol = 1)
# Covariates 
x1 <- rbinom(num_subj, 1, 0.5)
x2 <- floor(runif(num_subj, 50, 75))
X <- model.matrix(~rep(x1, each = num_times)+rep(x2,each = num_times))
# Random effects and error term
u <- rnorm(num_subj, 0, sigma.1)
eps <- rnorm(N, 0, sigma.0)
Y = X %*% beta + rep(u, each = num_times) + eps

long.data <- data.frame(id = rep(1:num_subj, each = num_times),
                        x1 = rep(x1, each = num_times),
                        x2 = rep(x2, each = num_times),
                        Y = Y)
actual.fit <- lmer(Y ~ x1 + x2 + (1|id), data = long.data, REML = F)
summary(actual.fit)
actual.ll <- summary(actual.fit)$logLik
Z <- as.matrix(getME(actual.fit, "Z"))

# Initial estimates -- Coefficients and parameters ----
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
  Vinv <- solve(V)
  
  # X %*% beta and the residuals (epsilon)
  Xb <- X %*% beta
  resid <- Y-Xb
  
  # The log-likelihood
  -0.5 * (t(resid) %*% Vinv %*% resid) - 0.5 * log(det(V))
}
ll(Y,X,Z,beta = c(40.08809, -9.75489, 1.49525), var.0 = 2.4323, var.1 = 0.1775) # Values from actual fit to check ...
actual.ll
# Classical EM - Single run -----------------------------------------------

n <- nrow(Y)
q1 <- nrow(Z)
l0 <- ll(Y,X,Z,beta, var.0, var.1)
i <- 0

# E step ----
# Covariance stuff
V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
Vinv <- solve(V)
# X %*% beta and resulting residuals
Xb <- X %*% beta
resid <- Y-Xb
temp1 <- Vinv %*% resid
# Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
# Old - not working.
# u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
#           tr(c(var.0) * diag(n) - c(var.0)^2 * t(diag(n)) %*% Vinv %*% diag(n))
u.0 <- c(var.0)^2 * t(temp1) %*% temp1 + c(var.0) * n - c(var.0)^2 * tr(Vinv)
u.1 <- c(var.1)^2 * t(temp1) %*% Z %*% t(Z) %*% temp1 + c(var.1) * q1 - c(var.1)^2 * tr(t(Z) %*% Vinv %*% Z)
# Old - not working.
# u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
#           tr(c(var.1) * diag(100) - c(var.1)^2 * t(Z) %*% Vinv %*% Z)


# Calculate expectation y-ZU conditional on y, call this w
w <- Xb + c(var.0) * Vinv %*% resid

 # M-step ----
# Estimate \sigma^(t+1)
var.0.new <- u.0/600
var.1.new <- u.1/600
# And estimate \beta
beta.new <- solve(t(X) %*% X) %*% t(X) %*% w
message(
  "old sigma.0 = ", sqrt(var.0), " new sigma.0 = ", sqrt(var.0.new), 
  "\nold sigma.1 = ", sqrt(var.1), " new sigma.1 = ", sqrt(var.1.new),
  "\nold beta: [", beta[1], beta[2], beta[3],"] new beta: [", beta.new[1], beta.new[2], beta.new[3],"]"
)
ll.new <- ll(Y, X, Z, beta.new, var.0.new, var.1.new)
ll2.new <- ll2(Y,X,Z,beta.new,var.0.new,var.1.new)
var.0 <- var.0.new; var.1 <- var.1.new; beta <- beta.new
# And then repeat the algorithm until convergence ! 
