#' #########
#' Doing N simulations of EM under different parameters
#' and initial conditions
#' dependent upon fns/classic.R and fns/simlong.R
#' #########

rm(list=ls())
library(lme4)
# Load scripts from function folder
setwd("~/Documents/PhD/EM-Algorithm")
source("./fns/Classic.R")
source("./fns/simlong.R")

# Single run - is it working?
set.seed(1995)
long <- simlong()
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)
# Run classic.em
classic.em(X, Y, Z, init.beta, init.var.0 = 1, init.var.1 = 1)
