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
var.0.init <- sd(x$Y - x$X %*% beta0)^2 # instead of just getting out of lmer.fit
var.1.init <- as.matrix(as.numeric(summary(x$lmer.fit)$varcor[1]))
theta0 <- c(t(beta0), var.0.init, var.1.init)

b.dens <- function(i){
  exp()/(s*pi*det)
}