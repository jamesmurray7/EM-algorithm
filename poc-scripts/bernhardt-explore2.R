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

ll <- function(X, Y, Z, beta, var.0, var.1, ni){ 
  V <- getV(var.1, Z, var.0, ni) 
  V.sqrt <- chol(V)
  detV <- sum(log(diag(V.sqrt))) * 2
  # X, beta and residuals
  Xb <- X %*% beta
  resid <- Y - Xb
  temp <- sum(forwardsolve(l = t(V.sqrt), x = resid) ^ 2)
  return(
    -0.5*length(Y)*log(2*pi) -0.5 * temp - 0.5 * detV
  )
}

# Maximisation function and its wrapper
b.ll <- function(b, Y, X, Z, beta, n){ # Note this is strictly 1-D.
  sigma <- b[2]
  -n/2 * log(2 * pi) - n/2 * log(sigma) - 1/(2*sigma) * crossprod(Y-X%*%beta-Z%*%b[1])
}

optim.b <- function(b0, Y, X, Z, beta, n){
  op <- optim(b0, b.ll, NULL, Y, X, Z, beta, n, 
              method = "L-BFG",
              control = list(fnscale = -1),
              lower = c(-5,0.01),
              upper = c(5,2))
  op
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
  Zi <- matrix(rep(1, 5), nc = 1) # We can cheat a little with Zi and ni as no missingness!
  ni <- 5
  
  # Sort parameter estimates out and intialise the parameter vector
  beta <- beta.init;
  if(!"matrix"%in%class(beta)){
    beta <- matrix(beta, nc = 1)
  }
  
  var.0 <- var.0.init; var.1 <- var.1.init
  
  b0 <- b.init
  sigma0 <- rep(1, length(uids))
  B <- sum(b0)/n
  
  params <- c(t(beta), var.0, B, var.1)
  
  # Begin while loop
  t0 <- proc.time()[3]
  while(diff > tol & iter <= maxiter){
    # E-step ----
    b.hat <- sigma.hat <- D.cont <- c()
    for(i in uids){
      start.val <- c(b0[i], sigma0[i])
      Yi <- getY(dat, i)
      Xi <- getX(dat, i)
      op <- optim.b(start.val, Yi, Xi, Zi, beta, ni)
      b.hat[i] <- op$par[1]
      sigma.hat[i] <- op$par[2]
    }
    
    # Rest of E-step, calculating E[eiei^T] to obtain var.0.new
    e <- nis <- c()
    for(i in uids){
      ni <- nrow(Zi); nis[i] <- ni
      V <- solve(getV(var.1, Zi, var.0, ni))
      residi <- getResid(data, i, beta)
      
      e[i] <- tr(
        c(var.0)^2 * V %*% tcrossprod(residi) %*% V + c(var.0) * (diag(ni) - c(var.0) * V)
      )
    }
    
    # M step ----
    # CM step 1
    # Updates for B and var.1
    B.new <- sum(b.hat)/n
    for(i in 1:length(b.hat)){
      #D.cont[i] <- sigma.hat[i] + tcrossprod(b.hat[i]-B.new)
      #D.cont[i] <- var(b.hat) + tcrossprod(b.hat[i]-B.new)
      #D.cont[i] <- 1/sigma.hat[i] * crossprod(Zi) + tcrossprod(b.hat[i] - B.new)
      D.cont[i] <- solve(1/sigma.hat[i] * crossprod(Zi) * -1) + tcrossprod(b.hat[i] - B.new)
    }
    var.1.new <- sum(D.cont)/n
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
    b0 <- b.hat
    sigma0 <- sigma.hat
    params <- params.new
    
    iter <- iter + 1
    
  }
  t1 <- proc.time()[3]
  message("EM Complete after ", iter, " iterations, this took ", round(t1-t0, 2), " seconds.")
  
  message("\nStarting post-hoc calculations of REs")
  REs <- c()
  for(i in uids){
    V <- solve(getV(var.1, Zi, var.0, ni))
    residi <- getResid(data, i, beta)
    REs[i] <- t(var.1 * Zi) %*% V %*% residi
  }
  message("done")
  
  if(logLik){
    message("\nCalculating log-likelihood")
    ll.temp <- c()
    for(i in uids){
      Yi <- getY(data, i)
      Xi <- getX(data, i)
      ll.temp[i] <-  ll(Xi, Yi, Zi, beta, var.0, var.1, ni)
    }
    message("done")
  }
  
  rtn <- list(
    beta = beta, var.0 = var.0, var.1 = var.1, B = B, last.b = b.hat,
    REs = REs, iter = iter, time=t1-t0
  )
  if(logLik) rtn <- c(rtn, logLik = sum(ll.temp))
  
  rtn
}

test <- em.bern(data = dat,
        beta.init = beta, 
        var.0.init = 1, var.1.init = 1, b.init = b.inits, tol = 1e-5)
plot(test$REs, ranef(x$lmer.fit)$id$`(Intercept)`); lines(-1:1,-1:1) # This is probably expected
sqrt(test$var.0) # target: 1.5
sqrt(test$var.1) # target: 0.5
summary(x$lmer.fit)$logLik; test$logLik # Broadly good agreement.
