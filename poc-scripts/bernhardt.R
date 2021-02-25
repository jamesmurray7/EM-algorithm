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
var.1.init <- as.matrix(as.numeric(summary(x$lmer.fit)$varcor[1]))
theta0 <- c(t(beta0), var.0.init, var.1.init)
# Initial estimates for bi
b.inits <- c()
uids <- unique(x$long.data$id)
for(i in uids){
  i.data <- subset(x$long.data, id == i)[1,]
  b.inits[i] <- i.data$Y - cbind(1,i.data$x1,i.data$x2) %*% beta0
}
# Estimate for B, 
B <- mean(b.inits)

# 2.maximise f(bi) to obtain hat(bi) and hat(Sigmai) ----------------------

# Density function
# Written atm to take in one value of b at a time, doesn't work otherwise.
b.dens <- function(b,b.input){
  b.diff <- as.matrix(b[1]-b.input)
  exp(-0.5 * t(b.diff) %*% solve(b[2]) %*% b.diff)/sqrt(2*pi*c(b[2]))
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

# Update parameter vector \bm{\theta} via EM steps ------------------------


B.new <- 
