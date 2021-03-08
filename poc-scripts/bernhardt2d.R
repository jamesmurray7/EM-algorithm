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

# Initial conditions from the simulated data
x <- simlong(100,5, Sigma = SigmaGen(2.5, 1.5))
dat <- x$long.data
beta <- lm(Y~x1+x2,dat)$coef
b.int.inits <- ranef(x$lmer.fit)$id$`(Intercept)`
b.slp.inits <- ranef(x$lmer.fit)$id$`time`
b.inits <- cbind(b.int.inits, b.slp.inits)
var.e <- 1.5
D <- SigmaGen(2.5,1.5)


# Setting out useful functions --------------------------------------------
getV <- function(D, Z, var.0, ni) Z %*% D %*% t(Z) + c(var.0) * diag(ni) # ni = nrow(getZ(dat, idx))

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

# Maximisation function and its wrapper
# Random intercept and slope - 

b.ll <- function(b.hat, Y, X, Z, n, D, var.e, beta){
  V <- getV(D, Z, var.e, n)
  bi <- matrix(c(b.hat[1], b.hat[2]), nc = 1)
  -n/2 * log(2 * pi) - n/2 * log(det(V)) - 1/2 * 
    t(Y - X %*% beta - Z %*% bi) %*% solve(V) %*% (Y - X %*% beta - Z %*% bi)
}

getd2 = function(Z, D, var.e){ # Calculate 2nd derivative of f(y|..., theta)
  ni <- nrow(Z)
  V <- solve(getV(D, Z, var.e, ni))
  -(0.5 * t(Z) %*% V %*% Z + 0.5 * t(V %*% Z) %*% Z)
}

optim.b <- function(b0, Y, X, Z, n, D, var.e, beta){
  op <- optim(b0, b.ll, NULL, Y, X, Z, n, D, var.e, beta, 
              control = list(fnscale = -1),
              method = "L-BFG",
              lower = c(-Inf,-Inf),
              upper = c(Inf, Inf))
  op
}

b1 <- b2 <- c() # This looks quite promising!
for(i in 1:100){
  Yi <- getY(dat, i); Xi <- getX(dat, i); Zi <- getZ(dat, i)
  ni <- nrow(Zi)
  op <- optim.b(c(0,0), Yi, Xi, Zi, ni, D, var.e, beta)
  b1[i] <- op$par[1]; b2[i] <- op$par[2]
}
# Quick plot to show how iteration zero shapes b
par(mfrow=c(2,1))
plot(b1, b.int.inits); lines(-10:10,-10:10) 
plot(b2, b.slp.inits); lines(-10:10,-10:10) 
par(mfrow=c(1,1))

getd1 <- function(X,Y,Z,D,var.e,beta,ni,b.hat){
  V <- solve(getV(D, Z, var.e, ni))
  temp <- Y - X%*% beta - Z %*% b.hat
  0.5 * t(Z) %*% V %*% temp + 0.5 * t(Z) %*% t(V) %*% temp
}

test <- list()
for(i in 1:100){
  print(i)
  Xi <- getX(dat, i); Yi <- getY(dat, i); Zi <- getZ(dat, i)
  ni <-  nrow(Zi)
  
  sigmai = solve(tcrossprod(getd1(Xi,Yi,Zi,D,var.e,beta,ni,b.hat[,i])))
  test[[i]] <- sigmai + tcrossprod(b.hat[,i] - B.mn)
}


em.bern <- function(data,
                    beta.init, var.0.init, 
                    D.init = NULL, b.init = NULL,
                    tol = 1e-3, maxiter = 200,
                    logLik = T){
  
  diff <- 100; iter <- 0
  
  # Data-specific
  uids <- unique(data$id)
  n <- length(uids)
  
  # Sort parameter estimates out and intialise the parameter vector
  beta <- beta.init;
  if(!"matrix"%in%class(beta)){
    beta <- matrix(beta, nc = 1)
  }
  
  var.0 <- var.0.init; 
  # Initial D-matrix
  if(is.null(D.init)) D <- diag(2) else D <- D.init
  if(!isSymmetric(D) | dim(D)[1] != dim(D)[2]) stop("D must be symmetric square matrix")
  
  if(is.null(b.init)) b0 <- matrix(rep(0, n * 2), nr = 2) else b0 <- b.init
  
  B <- rowMeans(b0)
  
  params <- c(t(beta), var.0, as.vector(B), as.vector(D))
  
  # Begin while loop
  t0 <- proc.time()[3]
  while(diff > tol & iter <= maxiter){
    # E-step ----
    b1 <- b2 <- c()
    for(i in uids){
      Yi <- getY(data, i); Xi <- getX(data, i); Zi <- getZ(data, i)
      ni <- nrow(Zi)
      op <- optim.b(b0[,i], Yi, Xi, Zi, ni, D, var.0, beta)
      b1[i] <- op$par[1]; b2[i] <- op$par[2]
    }
    b.hat <- rbind(b1, b2)
    
    # Rest of E-step, calculating E[eiei^T] 
    e <- nis <- c()
    for(i in uids){
      Zi <- getZ(data, i); ni <- nrow(Zi)
      V <- solve(getV(D, Zi ,var.0, ni))
      residi <- getResid(data, i, beta)
      
      e[i] <- tr(
        c(var.0)^2 * V %*% tcrossprod(residi) %*% V + c(var.0) * (diag(ni) - c(var.0) * V)
      )
      nis[i] <- ni
    }
    # M step ----
    
    # CM step 1
    # Updates for B and D
    B.new <- rowMeans(b.hat)
    
    D.store <- list()
    for(i in uids){
      Zi <- getZ(data, i); ni <- nrow(Zi)
      sigmai <- solve(-1 * getd2(Zi, D, var.0))
      # D.store[[i]] <- sigmai + tcrossprod(b.hat[,i] - B.new)
      D.store[[i]] <- tcrossprod(b.hat[,i]-B.new)
    }
    D.new <- Reduce('+', D.store)/n ; 
    var.0.new <- sum(e)/sum(nis);
    
    # CM step 2
    beta.p1 <- beta.p2 <- list()
    for(i in uids){
      Zi <- getZ(data, i); ni <- nrow(Zi)
      Yi <- getY(data, i); Xi <- getX(data, i);
      V <- solve(getV(D.new, Zi, var.0.new, ni))
      beta.p1[[i]] <- t(Xi) %*% V %*% Xi
      beta.p2[[i]] <- t(Xi) %*% V %*% Yi
    }
    
    beta.new <- solve(Reduce('+', beta.p1)) %*% Reduce('+', beta.p2)
    # Collecting parameters and finding max absolute difference
    params.new <- c(t(beta.new), var.0.new, B.new, as.vector(D.new))
    diffs <- abs(params-params.new)
    diff <- max(diffs)
    
    print(rbind(params,params.new))
    message("Iteration ", iter, " max difference: ", diff)
    
    # Set next parameter values
    beta <- beta.new
    var.0 <- var.0.new
    D <- D.new
    B <- B.new
    b0 <- b.hat
    params <- params.new
    
    iter <- iter + 1
    
  }
  t1 <- proc.time()[3]
  message("EM Complete after ", iter, " iterations, this took ", round(t1-t0, 2), " seconds.")
  
  message("\nStarting post-hoc calculations of REs")
  REs <- matrix(NA, nc = 2, nr = n)
  for(i in uids){
    Zi <- getZ(data, i); ni <- nrow(Zi)
    V <- solve(getV(D, Zi, var.0, ni))
    residi <- getResid(data, i, beta)
    REs[i,] <- t(Zi %*% D) %*% V %*% residi
  }
  message("done")
  
  if(logLik){
    message("\nCalculating log-likelihood at final parameter values\n")
    ll.temp <- c()
    for(i in uids){
      Zi <- getZ(data, i); ni <- nrow(Zi)
      Yi <- getY(data, i)
      Xi <- getX(data, i)
      ll.temp[i] <-  ll(Xi, Yi, Zi, beta, var.0, D, ni)
    }
    message("Log-likelihood calculated")
  }
  
  rtn <- list(
    beta = beta, var.0 = var.0, D = D, B = B, last.b = b.hat,
    REs = REs, iter = iter, time=t1-t0
  )
  if(logLik) rtn <- c(rtn, logLik = sum(ll.temp))
  
  rtn
}

# Do some comparison with lmer...
test <- em.bern(data = dat,
                beta.init = beta, 
                var.0.init = 1, D.init = NULL, 
                b.init = matrix(c(b.int.inits, b.slp.inits), nr = 2), tol = 1e-3)


par(mfrow=c(2,1))
plot(test$REs[, 1], ranef(x$lmer.fit)$id$`(Intercept)`); 
lines(-10:10,-10:10) # This is probably expected
plot(test$REs[, 2], ranef(x$lmer.fit)$id$time); lines(-10:10, -10:10)
par(mfrow=c(1,1))
sqrt(test$var.0) # target: 1.5
sqrt(as.vector(diag(test$D))) # Target: 2.5, 1.5
summary(x$lmer.fit)$logLik; test$logLik # Broadly good agreement.



