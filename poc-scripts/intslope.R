#' ###############
#' Run through of ECM on a random intercept and slope.
#' ###############

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
# Define functions used throughout
getZ <- function(data=dat, idx){
  i.dat <- subset(data, id == idx)
  Z <- cbind(1, i.dat$time)
  Z
}

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

source("~/Documents/PhD/EM-Algorithm/fns/simlong2.R")

x <- simlong(250, 10, Sigma = SigmaGen(1.5,0.5))
X <- x$X; Y <- x$Y; Z <- x$Z
dat <- x$long.data
# Starting values
beta.inits <- lm(Y~x1+x2,dat)$coef
var.0 <- matrix(sigma(x$lmer.fit)^2)
var.u <- attr(VarCorr(x$lmer.fit)$id, "stddev")^2
sd.u <- sqrt(var.u)
corr <- attr(VarCorr(x$lmer.fit)$id, "correlation")
D.init <- matrix(c(var.u[1], sd.u[1] * sd.u[2],
                   sd.u[1] * sd.u[2], var.u[2]),2,2,byrow = T) * corr

# Define EM function ------------------------------------------------------
em <- function(data,
               init.beta = beta0, init.var.0 = var.0, init.D = D.init,
               tol = 1e-3, maxiter = 2000, history = F, logLik = T){
  message("EM Algorithm ----\n")
  
  diff <- 100
  iter <- 0
  
  var.0 <- init.var.0; 
  beta <-  init.beta
  D <- init.D
  
  if(det(D) == 0){
    stop("D.init not positive definite")
  }
  
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }
  
  params <- c(t(beta), var.0, as.vector(D))
  if(history) iter.hist <- data.frame(iter, t(params))
  
  uids <- unique(data$id)

  n <- length(Y); q1 <- length(uids)
  
  t0 <- proc.time()[3]
  # Begin while loop
  while(diff > tol & iter <= maxiter){
    
    # E-step ----
    b <- list()
    e <- c()
    nis <- c()
    for(i in uids){
      Zi <- getZ(data, idx = i)
      ni <- nrow(Zi)
      V <- solve(getV(D, Zi, var.0, ni))
      residi <- getResid(data, idx = i, beta)
      
      bb <- tcrossprod(D, Zi) %*% V %*% residi %*% t(residi) %*% V %*% Zi %*% D + 
        D - tcrossprod(D, Zi) %*% V %*% Zi %*% D
      b[[i]] <- tcrossprod(D, Zi) %*% V %*% residi %*% t(residi) %*% V %*% Zi %*% D + 
        D - tcrossprod(D, Zi) %*% V %*% Zi %*% D
      
      e[i] <- tr(
        c(var.0)^2 * V %*% residi %*% t(residi) %*% V + 
          c(var.0) * (diag(ni) - c(var.0) * V)
      )

      nis[i] <- ni
    }
    
    # M step ----
    # CM-step 1
    D.new <- Reduce('+', b)/q1
    var.0.new <- sum(e)/sum(nis)
    
    # CM-step 2
    betas.p1 <- list()
    betas.p2 <- list()
    for(i in uids){
      Zi <- getZ(data, idx = i)
      ni <- nrow(Zi)
      V <- solve(getV(D.new, Zi, var.0.new, ni))
      Xi <- getX(data, idx = i)
      Yi <- getY(data, idx = i)
      betas.p1[[i]] <- t(Xi) %*% V %*% Xi
      betas.p2[[i]] <- t(Xi) %*% V %*% Yi
    }
      
    beta.new <- solve(Reduce('+', betas.p1)) %*% Reduce('+', betas.p2)
    
    # New parameter set
    params.new <- c(t(beta.new), var.0.new, as.vector(D.new))
    diffs <- abs(params.new-params)
    diff <- max(diffs)
    
    message("Iteration ", iter, " largest difference = ", diff)
    
    # Set new parameters
    var.0 <- var.0.new
    D <- D.new
    beta <- beta.new
    iter <- iter + 1
    params <- params.new 
    
    # Record history if needed
    if(history) iter.hist <- rbind(iter.hist, c(iter, t(params.new)))
  }
  t1 <- proc.time()[3]
  message("EM Complete after ", iter, " iterations, this took ", round(t1-t0,2), " seconds.")
  
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
  
  message("\nStarting post-hoc calculation of random effects\n")
  bi <- matrix(NA, nc = 2, nr = q1)
  for(i in uids){
    Zi <- getZ(data, i); ni <- nrow(Zi)
    V <- solve(getV(D, Zi, var.0, ni))
    residi <- getResid(data, i, beta)
    bi[i,] <- t(Zi %*% D) %*% V %*% residi
  }
  message("Done")
  
  # Collecting things to output
  rtn <- list(data = data, beta = beta, var.0 = var.0, D = D, REs = bi, time = t1-t0, iter = iter)
  if(history) rtn <- c(rtn, iter.hist = iter.hist)
  if(logLik) rtn <- c(rtn, logLik = sum(ll.temp))
  
  rtn
}

em(data = dat, init.beta = beta.inits, init.var.0 = 1, init.D = D.init, tol=1e-3)
