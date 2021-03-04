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
  b.diff <- as.matrix(b.input-b[1])
  exp(-0.5 * t(b.diff) %*% solve(b[2]) %*% b.diff)/sqrt(2*pi*det(matrix(b[2])))
}

b.hat <- c()
Sigma.hat <- c()
for(i in 1:length(b.inits)){
  op <- optim(c(5,5), b.dens, gr = NULL, b.inits[i],
        control = list(fnscale = -1),
        method = "L-BFGS", lower = c(-15, 0.00001), upper = c(15, 10))
  b.hat[i] <- op$par[1]
  Sigma.hat[i] <- op$par[2]
}

# Try univarite normal instead of the MVN??
b.dens.univ <- function(b, b.input){
  b.diff <- b.input-b[1]
  -1/(sqrt(2 * pi) * b[2]) * exp(-0.5 * (b.diff/b[2])^2)
}

nlminb(c(0,1), b.dens.univ, NULL, NULL, b.inits[5], lower = c(-4,0.01), upper = c(4,2))
# Same results as initial estimates so doing something wrong here!

# Is it supposed to be log-likelihood?
b.ll <- function(b, b.input){
  # b[1] = mu; b[2] = sigma^2
  ll <- -1/2 * log(2 * pi) - 1/2 * log(b[2]) - 1/(2*b[2]) * (b.input-b[1])^2
}

b.ll2 <- function(b, b.input){
  b1 <- b[1:(length(b)-1)]
  sigma <- b[length(b)]
  n <- length(b.inits)
  sum(-n/2 * log(2*pi) - 1/2 * log(sigma) - 1/(2*sigma) * (b.input-b1)^2)
}

b.ll3 <- function(b, b.input){ # Maybe this
  sigma <- (b.input-b) * t(b.input-b)/5; 
  sum(-5/2 * log(2*pi) - 1/2 * log(sigma) - 1/(2*sigma) * (b.input-b)^2)
}
b0 <- rep(0,100)

b.ll3(b0, rep(1,100))
op <- optim(b0, b.ll3, NULL, b.inits, method = "L-BFG",
            control = list(fnscale = -1),
            lower = c(rep(-5,100)),
            upper = c(rep(5,100)), hessian = T)
op$hessian

Sigma.i <- diag((b.inits-new.b) %*% t(b.inits-new.b)/5)

b0 <- c(b.inits,var.0.init)
optim(b0, b.ll2, NULL, b.inits, method = "L-BFG",
      lower = c(rep(min(b.inits), 100, 0)),
      upper = c(rep(max(b.inits), 100, 0)))

op <- optim(b0, b.ll2, NULL, rep(0,100), method = "L-BFG",
      control = list(fnscale = -1),
      lower = c(rep(min(b.inits), 100), 0.01),
      upper = c(rep(max(b.inits), 100), 5))

b.hat <- c()
sigma.hat <- c()
for(i in 1:length(b.inits)){
  op <- optim(c(0,0.5), b.ll, NULL, b.inits[i],
              control = list(fnscale = -1), method = "L-",
              lower = c(min(b.inits), 0.01), upper = c(max(b.inits), 2))
  b.hat[i] <- op$par[1]; sigma.hat[i] <- op$par[2]
}

# Is it supposed to be parameter estimate of D from hat(bi)?
laird <- function(d){
  d <- as.matrix(d)
  Sigma <-  Zi %*% d %*% t(Zi) + c(var.0.init) * diag(5)
  bi.hat <- t(Zi %*% d) %*% solve(Sigma) %*% (Yi - Xi%*%beta0)
  bi.hat
}

uids <- unique(x$long.data$id)
Zi <- rep(1,5)
dat <- x$long.data
for(i in uids){
  i.dat <- subset(dat, id == i)
  Yi <- i.dat$Y
  Xi <- model.matrix(~x1+x2, data = i.dat)
  op <- optim(var.1.init, laird,
              lower=0.01, upper = 2,
              method = "Brent", control = list(fnscale=-1))
  print(op)
}

# Update parameter vector \bm{\theta} via EM steps ------------------------
B.new <- mean(b.hat)
