#' #########
#' Bernhardt.R
#' Proof of concept for the approach taken by P. Bernhardt (2015)
#' Applied to a simple LME model, so may not be super useful.
#' #########

rm(list=ls())
setwd("~/Documents/PhD/EM-Algorithm/fns/")
source("simlong.R")

# 1. Using all data, get initial estimates for \bm{\theta}
x <- simlong()
beta0 <- x$lmer.fit@beta
var.0.init <- sd(x$Y - x$X %*% beta0)^2 # instead of just getting out of lmer.fit
var.1.init <- as.numeric(summary(x$lmer.fit)$varcor[1])
theta0 <- c(beta0, var.0.init, var.1.init)

# 2. Maximise u0u0t and u1u1t using stats::optim
tr <- function(x) sum(diag(x))

u0u0 <- function(theta){
  n <- length(Y); q <- ncol(Z)
  beta <- matrix(theta[1:3], ncol = 1)
  V <- c(theta[5]) * tcrossprod(Z) + c(theta[4]) * diag(n)
  V.inv <- solve(V)
  # X %*% beta and resulting residuals
  Xb <- X %*% beta
  resid <- Y-Xb
  # Calculate expectations u0u0t
  
  c(theta[4])^2 * crossprod(resid, V.inv) %*% V.inv %*% resid +
    tr(c(theta[4]) * diag(n) - c(theta[4])^2 * V.inv)
}
X <- x$X; Y <- x$Y; Z <- x$Z
u0u0(theta0)

optim(par = theta0, u0u0, method = "L-BFGS-B",#, lower = c(35, -15, 1, 0.1, 0.1), upper = c(45, -5, 2, 3, 1),
      control = list(fnscale = -1, REPORT = 1))
nlm(u0u0, theta0,fscale=-1)
