#' ##########
#' Benchmarking different approaches to computation
#' ##########

rm(list = ls())
library(rbenchmark)
library(lme4)
source("~/Documents/PhD/EM-Algorithm/fns/simlong.R")
long.data <- simlong()
X <- long.data$X; Y <- long.data$Y; Z <- long.data$Z
beta <- matrix(long.data$lmer.fit@beta, ncol = 1)

V <- c(1) * Z %*% t(Z) + c(1) * diag(length(Y))

# Quickest way to find inverse V^-1
R <- chol(V)
benchmark(
  "Solve" = {
    Vinv <- solve(V)
  },
  "Solve-Crossprod" = {
    RtR.inv <- solve(crossprod(R))
  },
  "RinvRinvT" = {
    Rinv <- solve(R)
    RinvRinvT <- Rinv %*% t(Rinv)
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self"))
# Base appears to be fastest (?)
# test replications elapsed relative user.self sys.self
# 3       RinvRinvT           10  22.131    4.328    21.977    0.112
# 1           Solve           10   5.113    1.000     5.054    0.045
# 2 Solve-Crossprod           10  16.244    3.177    16.033    0.104


# Is XtX slower than crossprod(X)? ----
benchmark(
  "XtX" = {
    XtX.inv <- solve(t(X) %*% X)
  },
  "crossprodX" = {
    XtX.inv2 <- solve(crossprod(X))
  },
  replications = 1000,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)
# Crossprod marginally faster
# test replications elapsed relative user.self sys.self
# 2 crossprodX         1000    0.05      1.0     0.049    0.000
# 1        XtX         1000    0.07      1.4     0.062    0.008

# t(Z) %*% Vinv %*% Z vs forwardsolve ----
Vinv <- solve(V)
R <- chol(V)
benchmark(
  "base" = {
    f0 <- t(Z) %*% Vinv %*% Z
  },
  "fwdsolve" = {
    z <- forwardsolve(t(R),Z)
    f1 <- crossprod(z)
  },
  replications = 100,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)
# forwardsolve appears to be much faster.
# test replications elapsed relative user.self sys.self
# 1     base          100  37.849    6.949    37.337    0.364
# 2 fwdsolve          100   5.447    1.000     5.306    0.125
