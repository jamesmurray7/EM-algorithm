#' ##########
#' Benchmarking different approaches to computation
#' ##########

rm(list = ls())
library(rbenchmark)
source("~/Documents/PhD/EM-Algorithm/fns/simlong.R")
long.data <- simlong()
X <- long.data$X; Y <- long.data$Y; Z <- long.data$Z
beta <- matrix(long.data$lmer.fit@beta, ncol = 1)

V <- c(1) * Z %*% t(Z) + c(1) * diag(length(Y))

# Quickest way to find inverse V^-1
benchmark(
  "Solve" = {
    Vinv <- solve(V)
  },
  "Solve-Crossprod" = {
    R <- chol(V)
    RtR.inv <- solve(crossprod(R))
  },
  "RinvRinvT" = {
    R <- chol(V)
    Rinv <- solve(R)
    RinvRinvT <- Rinv %*% t(Rinv)
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self"))
# Base appears to be fastest (?)

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
