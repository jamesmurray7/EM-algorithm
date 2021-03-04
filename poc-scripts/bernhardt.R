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
uids <- unique(dat$id); n <- length(uids)
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
sigma.i <- (b.inits-b.new) * (b.inits-b.new)/5

# 3. Continue with EM
# Importing helper functions from intslope.R
getZ <- function(data=dat, idx){
  i.dat <- subset(data, id == idx)
  Z <- cbind(1, i.dat$time)
  Z
}

getV <- function(var.1, Z, var.0, ni)  c(var.1) * tcrossprod(Z) + c(var.0) * diag(ni) # ni = nrow(getZ(dat, idx))

getResid <- function(data=dat, idx, beta){
  i.dat <- subset(data, id == idx)
  Xi <- model.matrix(~i.dat$x1 + i.dat$x2)
  Yi <- i.dat$Y
  Yi - Xi %*% beta
}

getX <- function(data = dat, idx){
  i.dat <- subset(data, id == idx)
  model.matrix(~i.dat$x1 + i.dat$x2)
}

getY <- function(data = dat, idx){
  i.dat <- subset(data, id == idx)
  matrix(i.dat$Y, nc = 1)
}

tr <- function(x) sum(diag(x))

B <- sum(b.new)/length(b.new)
var.1.new <- crossprod(b.new-B)/n
var.0 <- var.0.init
beta <- beta0

# Rest of the E-step ----
e <- nis <-  c()
for(i in uids){
  Zi <- matrix(rep(1,5),nc = 1); ni <- nrow(Zi)
  V <- solve(getV(var.1.new, Zi, var.0, ni))
  residi <- getResid(dat, idx = i, beta = beta0)
  e[i] <- tr(
    c(var.0)^2 * V %*% residi %*% t(residi) %*% V + c(var.0) * (diag(ni) - c(var.0) * V)
  )
  nis[i] <- ni
}

# M step ----
# CM step 1 //
var.0.new <- sum(e)/sum(nis)

# CM step 2 //
beta.p1 <- beta.p2 <- list()
for(i in uids){
  Zi <- matrix(rep(1,5), nc = 1); ni <- nrow(Zi)
  Yi <- getY(dat, i); Xi <- getX(dat, i);
  V <- solve(getV(var.1.new, Zi, var.0.new, ni))
  beta.p1[[i]] <- t(Xi) %*% V %*% Xi
  beta.p2[[i]] <- t(Xi) %*% V %*% Yi
}

beta.new <- solve(Reduce('+', beta.p1)) %*% Reduce('+', beta.p2)
# and then loop

