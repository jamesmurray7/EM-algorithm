#' #########
#' Bernhardt2D.R
#' Extending bernhardt-explore to model on random intercept and slope.
#' #########

rm(list=ls())
library(lme4)

SigmaGen <- function(int.sd, slp.sd){
  Sigma <- matrix(
    c(int.sd * int.sd, int.sd * slp.sd,
      slp.sd * int.sd, slp.sd * slp.sd),
    2, 2, byrow = T)
  rho <-  matrix(c(1,0.2,0.2,1), 2, 2, byrow = T)
  
  Sigma * rho
}

simlong <- function(num_subj = 250, num_times = 5, 
                    beta = c(40, -10, 1.5), sigma.0 = 1.5, Sigma){
  
  N <- num_subj * num_times
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }
  # Covariates 
  x1 <- rbinom(num_subj, 1, 0.5)
  x2 <- floor(runif(num_subj, 50, 75))
  X <- model.matrix(~rep(x1, each = num_times)+rep(x2, each = num_times))
  # Set up time 
  max.time <- num_times - 1
  t <- 0:max.time
  # Random effects and error term
  U <- MASS::mvrnorm(num_subj, c(0, 0), Sigma = Sigma)
  u1 <- U[,1]; u2 <- U[,2]
  # Epsilon term
  eps <- rnorm(N, 0, sigma.0)
  
  Y <- X %*% beta + rep(u1, each = num_times) + rep(t, num_subj) * rep(u2, each = num_times) + eps
  
  long.data <- data.frame(id = rep(1:num_subj, each = num_times),
                          time = rep(t, num_subj),
                          x1 = rep(x1, each = num_times),
                          x2 = rep(x2, each = num_times),
                          Y = Y)
  # Actual fit using lme4::lmer
  actual.fit <- lmer(Y ~ x1 + x2 + (1 + time|id), data = long.data, REML = F)
  Z <- as.matrix(getME(actual.fit, "Z"))
  return(
    list(X = X, Y = Y, Z = Z, long.data = long.data, lmer.fit = actual.fit)
  )
}

x <- simlong(100,5, Sigma = SigmaGen(2.5, 1.5))
dat <- x$long.data
beta <- lm(Y~x1+x2,dat)$coef
b.int.inits <- ranef(x$lmer.fit)$id$`(Intercept)`
b.slp.inits <- ranef(x$lmer.fit)$id$`time`
b.inits <- cbind(b.int.inits, b.slp.inits)

# Setting out useful functions --------------------------------------------
getV <- function(D, Z, var.0, ni) t(Z) %*% D %*% Z + c(var.0) * diag(ni) # ni = nrow(getZ(dat, idx))

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

getZ <- function(data, idx){
  i.dat <- subset(data, id == idx)
  cbind(1, i.dat$time)
}

tr <- function(x) sum(diag(x))

ll <- function(X, Y, Z, beta, var.0, D, ni){ 
  V <- getV(D, Z, var.0, ni) 
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

# Dummy turn
Z1 <- getZ(dat, 1); X1 <- getX(dat, 1); Y1 <- getY(dat ,1)
b1 <- b.inits[1,]

# Maximisation function and its wrapper
# Random intercept and slope - 

b.ll <- function(b, bi, Y, X, beta, Z, var.e, n){ 
  # Just the square D matrix?
  D <- matrix(c(b[3], 0, 0, b[4]), 2, 2)
  V <- Z %*% D %*% t(Z) + var.e * diag(n)
  Sigma <- D - D %*% t(Z) %*% solve(V) %*% Z1 %*% D
  mu.b <- D %*% t(Z) %*% solve(V) %*% (Y - X %*% beta)
  
  -n/2 * log(2 * pi) -n/2 * log(det(Sigma)) + 
    t(c(b[1]-b[2]) - mu.b) %*% solve(Sigma) %*% (c(b[1],b[2]) - mu.b)
  
  # Say b0 and b1 are iid so off-diagonals just zero.
  #sigma <- Z %*% diag(c(b[3], b[4])) %*% t(Z) + var.e * diag(n)
  #sigma <- b[3]
  
  #-n/2 * log(2 * pi) - n/2 * log(det(sigma)) - 1/2 * t(Y - X %*% beta - Z %*% b[1:2]) %*% solve(sigma) %*% (Y-X %*% beta-Z %*% b[1:2])
}

mu.b <- matrix(NA, nr = nrow(b.inits), nc = 2)
d00 <- d11 <- c()
for(i in 1:nrow(b.inits)){
  Zi <- getZ(dat, i); Yi <- getY(dat, i); Xi <- getX(dat, i)
  op <- optim(c(1,1,1,1), b.ll, NULL, 
              b.inits[i,], Yi, Xi, beta, Zi, var.e, nrow(Zi),
              control = list(fnscale=-1),
              method = "L-",
              lower = c(-5, -5, 1e-3, 1e-3),
              upper = c(5, 5, 2, 2))
  mu.b[i,] <- op$par[1:2]
  d00[i] <- op$par[3]; d11[i] <- op$par[4]
}

optim.b <- function(b0, Y, X, Z, beta, var.e, n){
  op <- optim(b0, b.ll, NULL, Y, X, Z, beta, var.e, n, 
              control = list(fnscale = -1),
              method = "L-BFG",
              lower = c(-5,-5, 0.01),
              upper = c(5, 5, Inf))
  op
}

b1 <- b2 <- d1 <- d2 <- c()
for(i in 1:nrow(b)){
  Yi <- getY(dat, i); Xi <- getX(dat, i); Zi <- getZ(dat,i); ni <- nrow(Zi)
  op <- optim.b(c(0, 0, 1), Yi, Xi, Zi, beta, 1, ni)
  b1[i] <- op$par[1]
  b2[i] <- op$par[2]
  d1[i] <- op$par[3]
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
  ni <- nrow(data)/n
  Zi <- matrix(rep(1, ni), nc = 1) # We can cheat a little with Zi and ni as no missingness!
  
  
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

# Do some comparison with lmer...
test <- em.bern(data = dat,
                beta.init = beta, 
                var.0.init = 1, var.1.init = 1, b.init = rep(0,250), tol = 1e-5)
plot(test$REs, ranef(x$lmer.fit)$id$`(Intercept)`); lines(-1:1,-1:1) # This is probably expected
sqrt(test$var.0) # target: 1.5
sqrt(test$var.1) # target: 0.5
summary(x$lmer.fit)$logLik; test$logLik # Broadly good agreement.
