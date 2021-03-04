#' #########
#' Bernhardt.R
#' Proof of concept for the approach taken by P. Bernhardt (2015)
#' Applied to a random intercept model.
#' #########

rm(list=ls())
library(lme4)
setwd("~/Documents/PhD/EM-Algorithm/fns/")
source("simlong.R")

# 1. Using all data, get initial estimates for \bm{\theta}
x <- simlong(num_subj = 100, num_times = 5)
X <- x$X; Y <- x$Y; Z <- x$Z; dat <- x$long.data
beta0 <- matrix(lm(Y~x1+x2,data=x$long.data)$coef,nc=1)
var.0.init <- as.matrix(sigma(x$lmer.fit)^2)
var.1.init <- as.matrix(as.numeric(summary(x$lmer.fit)$varcor[1])) # = D in univariate case.
theta0 <- c(t(beta0), var.0.init, var.1.init)

# Initial estimates for bi
b.inits <- ranef(x$lmer.fit)$id$`(Intercept)`
bi.hat <- rep(0, length(b.inits))

# 2. Maximise f(bi) to obtain hat{bi} and hat{Sigmai} ---------------------

b.ll <- function(b, b.input){ # Note this is strictly 1-D.
  sigma <- (b.input-b) * t(b.input-b)/5; # May have to set-up n_i here in the future.
  sum(-5/2 * log(2*pi) - 1/2 * log(sigma) - 1/(2*sigma) * (b.input-b)^2)
}

op <- optim(bi.hat, b.ll, NULL, b.inits, method = "L-BFG",
            control = list(fnscale = -1),
            lower = c(rep(-5,100)),
            upper = c(rep(5,100)))

b.new <- op$par
sigma.i <- (b.inits-b.)

# 3. Continue with EM

var.1.new <- 
