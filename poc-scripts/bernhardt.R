#' #########
#' Bernhardt.R
#' Proof of concept for the approach taken by P. Bernhardt (2015)
#' Applied to a simple LME model, so may not be super useful.
#' #########

rm(list=ls())
library(lme4)
setwd("~/Documents/PhD/EM-Algorithm/fns/")
source("simlong.R")

# 1. Using all data, get initial estimates for \bm{\theta}
x <- simlong(num_subj = 100, num_times = 5)
X <- x$X; Y <- x$Y; Z <- x$Z
beta0 <- x$lmer.fit@beta
var.0.init <- as.matrix(sigma(x$lmer.fit)^2)
var.1.init <- as.matrix(as.numeric(summary(x$lmer.fit)$varcor[1])) # = D in univariate case.
theta0 <- c(t(beta0), var.0.init, var.1.init)
# Initial estimates for bi
b.inits <- ranef(x$lmer.fit)$id$`(Intercept)`
# Estimate for B, 
B <- mean(b.inits)

# 2.maximise f(bi) to obtain hat(bi) and hat(Sigmai) ----------------------

# Density function
# Written atm to take in one value of b at a time, doesn't work otherwise.
b.dens <- function(b,b.input){
  b.diff <- as.matrix(b[1]-b.input)
  exp(-0.5 * t(b.diff) %*% solve(b[2]) %*% b.diff)/sqrt(2*pi*det(matrix(b[2])))
}

b.hat <- c()
Sigma.hat <- c()
for(i in 1:length(b.inits)){
  op <- optim(c(0,1), b.dens, gr = NULL, b.inits[i],
        control = list(fnscale = -1),
        method = "L-BFGS", lower = c(-5,0.01), upper = c(5,3))
  b.hat[i] <- op$par[1]
  Sigma.hat[i] <- op$par[2]
}

# Try univarite normal instead of the MVN??
b.dens.univ <- function(b, b.input){
  b.diff <- b.input-b[1]
  -1/(sqrt(2 * pi) * b[2]) * exp(-0.5 * (b.diff/b[2])^2)
}

nlminb(c(0,1), b.dens.univ, NULL, NULL, b.inits[1], lower = c(-4,0.01), upper = c(4,2))
# Same results as initial estimates so doing something wrong here!

# Is it supposed to be log-likelihood?
b.ll <- function(b, b.input){
  # b[1] = mu; b[2] = sigma^2
  ll <- -1/2 * log(2 * pi) - 1/2 * log(b[2]) - 1/(2*b[2]) * (b.input-b[1])^2
}

b.hat <- c()
sigma.hat <- c()
for(i in 1:length(b.inits)){
  op <- optim(c(0,0.5), b.ll, NULL, b.inits2[i],
              control = list(fnscale = -1), method = "L-",
              lower = c(min(b.inits), 0.01), upper = c(max(b.inits), 2))
  b.hat[i] <- op$par[1]; sigma.hat[i] <- op$par[2]
}

# Update parameter vector \bm{\theta} via EM steps ------------------------
B.new <- mean(b.hat)
