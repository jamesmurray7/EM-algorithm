#' ###########
#' Bern-ECM-simstudy.R
#' --
#' Simulation study on ECM vs Bernhardt's approach on intercept + slope model.
#' > Benchmarking
#' > Different starting conditions
#' > Different N, n_i
#' ###########

rm(list=ls())
library(lme4)
library(rbenchmark)

# Load in sources
source("~/Documents/PhD/EM-Algorithm/fns/simlong2.R")
source("~/Documents/PhD/EM-Algorithm/fns/bernhardt2d.R")
source("~/Documents/PhD/EM-Algorithm/fns/ECM2d.R")


# Benchmarking ------------------------------------------------------------

# Test data - small n
set.seed(1995)
x <- simlong(100, 5, Sigma = SigmaGen(2.5, 1.5))
dat <- x$long.data
beta <- lm(Y~x1+x2,dat)$coef

benchmark(
  "ECM" = {
    ecm <- em(data = dat, init.beta = beta, init.var.0 = 1, tol = 1e-3)
  },
  "Bernhardt" = {
    bern <- em.bern(data = dat, beta.init = beta, var.0.init = 1, tol = 1e-3)
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# Bernhardt slower on smaller samples - i.e. no material gain!
#         test replications elapsed relative user.self sys.self
# 2 Bernhardt           10  41.351    1.586    40.938    0.261
# 1       ECM           10  26.070    1.000    25.817    0.178

# Larger sample
set.seed(1995)
x <- simlong(500, 10, Sigma = SigmaGen(2.5, 1.5))
dat <- x$long.data
beta <- lm(Y~x1+x2,dat)$coef

benchmark(
  "ECM" = {
    ecm <- em(data = dat, init.beta = beta, init.var.0 = 1, tol = 1e-3)
  },
  "Bernhardt" = {
    bern <- em.bern(data = dat, beta.init = beta, var.0.init = 1, tol = 1e-3)
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# Poor starting values
set.seed(1995)
x <- simlong(100,5, Sigma = SigmaGen(2.5,1.5))
dat <- x$long.data
beta <- c(22, -30, 0)

benchmark(
  "ECM" = {
    ecm <- em(data = dat, init.beta = beta, init.var.0 = 1, tol = 1e-3, init.D = SigmaGen(2.5,1.5))
  },
  "Bernhardt" = {
    bern <- em.bern(data = dat, beta.init = beta, var.0.init = 1, tol = 1e-3, D.init = SigmaGen(2.5,1.5))
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# Focus in on Bernhardt - change optimisers -------------------------------
rm(list=ls())
source("~/Documents/PhD/EM-Algorithm/fns/bernhardt2d.R")
source("~/Documents/PhD/EM-Algorithm/fns/simlong2.R")
x <- simlong(100,5,Sigma = SigmaGen(2.5,1.5))
beta <- c(40, -10, 5)
benchmark(
  "Bernhardt - default" = {
    bern1 <- em.bern(data = x$long.data, beta.init = beta, var.0.init = 1)
  },
  "Bernhardt - good b.init" = {
    bern2 <- em.bern(data = x$long.data, beta.init = beta,
                     var.0.init = 1, 
                     b.init = rbind(
                       ranef(x$lmer.fit)$id$`(Intercept)`,
                       ranef(x$lmer.fit)$id$`time`
                     ))
  },
  "Bernhardt - good D.init" = {
    bern3 <- em.bern(data = x$long.data,
                     beta.init = beta,
                     var.0.init = 1,
                     D.init = SigmaGen(2.5, 1.5))
  },
  "Bernhardt - all good inits" = {
    bern4 <- em.bern(data = x$long.data,
                     beta.init = beta, var.0.init = 1,
                     b.init = rbind(
                       ranef(x$lmer.fit)$id$`(Intercept)`,
                       ranef(x$lmer.fit)$id$`time`
                     ), D.init = SigmaGen(2.5, 1.5)
                     )
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)
