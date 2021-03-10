#' ###############
#' Expectation conditional-maximisation
#' ###############

# Helper functions
tr <- function(x) sum(diag(x))
ll <- function(X, Y, Z, beta, var.0, var.1){ 
  V <- c(var.1) * tcrossprod(Z) + c(var.0) * diag(nrow(Z)) # sigma.u * ZZt + sigma.eIn
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

getX <- function(data, idx) model.matrix(~x1+x2, data = data[data$id == idx,])
getY <- function(data, idx) matrix(data[data$id == idx, ]$Y, nc=1)
getResid <- function(data, idx, beta) getY(data, idx) - getX(data, idx) %*% beta
getZ <- function(data, idx) matrix(rep(1, nrow(data[data$id == idx, ])), nc=1)
getV <- function(Z, D, var.0) D * tcrossprod(Z) + var.0 * diag(nrow(Z))

# ECM function 
ecm <- function(data, beta.init, var.0.init, D.init,
                tol = 1e-3, maxiter = 2e3, history = F, logLik = T){
  message("ECM Algorithm ----\n")
  
  diff <- 100
  iter <- 0
  uids <- unique(data$id)
  n <- length(uids)
  
  var.0 <- var.0.init
  D <- D.init
  beta <- beta.init
  
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }
  
  params <- c(t(beta), var.0, D)
  if(history) iter.hist <- data.frame(iter, t(params))
  
  t0 <- proc.time()[3]
  # Begin while loop
  while(diff > tol & iter <= maxiter){
    
    # E step ----
    b <- e <- nis <- c()
    for(i in uids){
      Zi <- getZ(data, i); ni <- nrow(Zi)
      residi <- getResid(data, i, beta)
      V <- solve(getV(Zi, D, var.0))
      b[i] <- D * crossprod(Zi, V) %*% tcrossprod(residi) %*% V %*% Zi * D + D - D * crossprod(Zi, V) %*% Zi * D
      e[i] <- tr(
        c(var.0) ^ 2 * V %*% tcrossprod(residi) %*% V + c(var.0) * (diag(ni) - c(var.0) * V)
      )
      nis[i] <- ni
    }
    
    # M step ----
    D.new <- sum(b)/n
    var.0.new <- sum(e)/sum(nis)
    beta.1 <- beta.2 <- list()
    for(i in uids){
      Xi <- getX(data, i); Yi <- getY(data, i); Zi <- getZ(data, i)
      V <- solve(getV(Zi, D.new, var.0.new))  
      beta.1[[i]] <- t(Xi) %*% V %*% Xi
      beta.2[[i]] <- t(Xi) %*% V %*% Yi
    }
    beta.new <- solve(Reduce('+', beta.1)) %*% Reduce('+', beta.2)
    
    # New parameter set
    params.new <- c(t(beta.new), var.0.new, D.new)
    diffs <- abs(params.new-params)
    diff <- max(diffs)
    
    message("Iteration ", iter, " largest difference = ", diff)
    
    # Set new parameters
    var.0 <- var.0.new
    D <- D.new
    beta <- beta.new
    iter <- iter + 1
    params <- c(t(beta), var.0, D) 
    
    # Record history if needed
    if(history) iter.hist <- rbind(iter.hist, c(iter, t(params.new)))
  }
  t1 <- proc.time()[3]
  message("ECM Complete after ", iter, " iterations, this took ", round(t1-t0,2), " seconds.")
  
  message("\nStarting post-hoc calculations of REs")
  REs <- c()
  for(i in uids){
    Zi <- getZ(data, i)
    V <- solve(getV(Zi, D, var.0))
    residi <- getResid(data, i, beta)
    REs[i] <- D * t(Zi) %*% V %*% residi
  }
  message("done")
  
  if(logLik){
    message("\nCalculating log-likelihood at final parameter values\n")
    ll.temp <- c()
    for(i in uids){
      Zi <- getZ(data, i); Yi <- getY(data, i); Xi <- getX(data, i)
      ll.temp[i] <-  ll(Xi, Yi, Zi, beta, var.0, D)
    }
    message("Log-likelihood calculated")
  }
  
  rtn <- list(
    beta = beta, var.0 = var.0, D = D, REs = REs, iter = iter, time=t1-t0
  )
  if(logLik) rtn <- c(rtn, logLik = sum(ll.temp))
  if(history) rtn <- c(rtn, iter.hist = iter.hist)
  
  rtn
}
