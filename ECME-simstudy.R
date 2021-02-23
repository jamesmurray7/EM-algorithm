#' #########
#' Doing N simulations of ECME under different parameters
#' and seeing if it's (broadly) faster than EM and ECM.
#' Dependencies:
#'  fns/Classic-paramtol.R 
#'  fns/simlong.R
#'  fns/ECM-paramtol.R
#'  fns/ECME.R
#' #########

rm(list=ls())
library(lme4)
library(tidyverse)
library(rbenchmark)
# Load scripts from function folder ----
setwd("~/Documents/PhD/EM-Algorithm/fns")
# To simulate longitudinal data
source("simlong.R") 
source("./ECMp2.R")
source("./Classic-paramtol.R")
source("./ECME.R")

# Small sample ------------------------------------------------------------
long <- simlong(num_subj = 100, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)

benchmark(
  "EM" = {
    em.res <- em(X,Y,Z,init.beta,1,1)
  },
  "ECM" = {
    ecm.res <- ecm(X,Y,Z,init.beta,1,1)
  },
  "ECME" = {
    ecme.res <- ecme(X,Y,Z,init.beta,1,1)
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)


# Medium sample
rm(long, em.res, ecm.res, ecme.res)
long <- simlong(num_subj = 300, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)
benchmark(
  "EM" = {
    em.res <- em(X,Y,Z,init.beta,1,1)
  },
  "ECM" = {
    ecm.res <- ecm(X,Y,Z,init.beta,1,1)
  },
  "ECME" = {
    ecme.res <- ecme(X,Y,Z,init.beta,1,1)
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# test replications elapsed relative user.self sys.self
# 2  ECM           10 523.567    1.000   509.836    5.441
# 3 ECME           10 607.633    1.161   577.381    7.137
# 1   EM           10 538.084    1.028   506.021    6.584

