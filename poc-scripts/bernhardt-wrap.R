#' #########
#' Bernhardt.R
#' Proof of concept for the approach taken by P. Bernhardt (2015) on RI model
#' Aiming to functionise the wrap, specifically for 1D problem such as this.
#' #########

rm(list=ls())
library(lme4)
setwd("~/Documents/PhD/EM-Algorithm/fns/")
source("simlong.R")


# Setting out useful functions --------------------------------------------
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

# Maximisation function and its wrapper
b.ll <- function(b, b.input){ # Note this is strictly 1-D.
  sigma <- (b.input-b) * t(b.input-b)/5
  sum(-5/2 * log(2*pi) - 1/2 * log(sigma) - 1/(2*sigma) * (b.input-b)^2)
}

optim.b <- function(current.b){
  op <- optim(rep(0,100), b.ll, NULL, current.b, method = "L-BFG",
              control = list(fnscale = -1),
              lower = c(rep(-5,100)),
              upper = c(rep(5,100)), hessian = T)
  return(list(
    par = op$par, hess = op$hessian
  ))
}

em.bern <- function(data,
                    beta.init, var.0.init, var.1.init,
                    b.init,
                    tol = 1e-3, maxiter = 200,
                    logLik = T){
  
  diff <- 100; iter <- 0
  
  # Data-specific
  uids <- unique(data$id)
  n <- length(uids)
  Zi <- matrix(rep(1, 5), nc = 1) # We can cheat a little as no missingness.
  
  # Sort parameter estimates out and intialise the parameter vector
  beta <- beta.init;
  if(!"matrix"%in%class(beta)){
    beta <- matrix(beta, nc = 1)
  }
  
  var.0 <- var.0.init; var.1 <- var.1.init
  
  b <- b.init
  b0 <- rep(0,length(b))
  B <- sum(b)/n
  
  params <- c(t(beta), var.0, B, var.1)
  
  # Begin while loop
  t0 <- proc.time()[3]
  while(diff > tol & iter <= maxiter){
    # E-step ----
    op <- optim.b(b)
    b.hat <- op$par; hess <- solve(op$hess)
    
    # Updates for B and var.1
    B.new <- sum(b.hat)/n
    D.cont <- c()
    for(i in 1:length(b.hat)){
      D.cont[i] <- hess[i,i] + tcrossprod(b.hat[i]-B.new)
    }
    var.1.new <- sum(D.cont)/n
    
    # Rest of E-step, calculating E[eiei^T] to obtain var.0.new
    e <- nis <- c()
    for(i in uids){
      ni <- nrow(Zi); nis[i] <- ni
      V <- solve(getV(var.1.new, Zi, var.0, ni))
      residi <- getResid(data, i, beta)
      
      e[i] <- tr(
        c(var.0)^2 * V %*% tcrossprod(residi) %*% V + c(var.0) * (diag(ni) - c(var.0) * V)
      )
    }
    
    # M step ----
    # CM step 1
    var.0.new <- sum(e)/sum(nis)
    
    # CM step 2
    beta.p1 <- beta.p2 <- list()
    for(i in uids){
      ni <- nrow(Zi)
      Yi <- getY(dat, i); Xi <- getX(dat, i);
      V <- solve(getV(var.1.new, Zi, var.0.new, ni))
      beta.p1[[i]] <- t(Xi) %*% V %*% Xi
      beta.p2[[i]] <- t(Xi) %*% V %*% Yi
    }
    
    beta.new <- solve(Reduce('+', beta.p1)) %*% Reduce('+', beta.p2)
    
    # Collecting parameters and finding max absolute difference
    params.new <- c(t(beta.new), var.0.new, B.new, var.1.new)
    diffs <- abs(params-params.new)
    diff <- max(diffs)
    message("Iteration ", iter, " max difference: ", diff)
    
    # Set next parameter values
    beta <- beta.new
    var.0 <- var.0.new
    var.1 <- var.1.new
    B <- B.new
    # b0 <- b
    b <- b.hat
    params <- params.new
    
    iter <- iter + 1
    
  }
  t1 <- proc.time()[3]
  message("EM Complete after ", iter, " iterations, this took ", round(t1-t0, 2), " seconds.")
  
  return(list(
    beta = beta, var.0 = var.0, var.1 = var.1, B = B, REs = b, iter= iter, time=t1-t0
  ))
}

# Test run on small set of simulated data
source("~/Documents/PhD/EM-Algorithm/fns/simlong.R")
x <- simlong(num_subj = 100, num_times = 5)
dat <- x$long.data
beta0 <- matrix(lm(Y~x1+x2,data=x$long.data)$coef,nc=1)
var0 <- as.matrix(sigma(x$lmer.fit)^2)
var1 <- as.matrix(as.numeric(summary(x$lmer.fit)$varcor[1])) # = D in univariate case.
b.inits <- ranef(x$lmer.fit)$id$`(Intercept)`

test.bern <- em.bern(data = dat,
        beta.init = beta0, 
        var.0.init = 1, var.1.init = 1, b.init = b.inits)

sqrt(test.bern$var.0)
sqrt(test.bern$var.1)

plot(test.bern$REs,b.inits)
