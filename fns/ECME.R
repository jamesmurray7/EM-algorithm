#' ###############
#' Expectation conditional-maximisation either (ECME)
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
ecme <- function(data,
                 beta.init, var.0.init, D.init,
                tol = 1e-3, maxiter = 2e3, history = F, logLik = T){
  diff <- 100; iter <- 0
  
  uids <- unique(data$id)
  n <- length(uids)
  
  # Initial conditions
  var.0 <- var.0.init; D <- D.init; beta <- beta.init
  
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }

  # First parameter vector
  params <- c(t(beta), var.0, D)
  
  if(history) iter.hist <- data.frame(iter, t(beta), var.0 = var.0, D = D)
  
  t0 <- proc.time()[3]
  while(diff > tol & iter <= maxiter){
    
    # E step ----
    
    bi <- Di <- e <- nis <- c()
    for(i in uids){
      Zi <- getZ(data, i)
      residi <- getResid(data, i, beta)
      bi[i] <- solve(var.0 * solve(D) + crossprod(Zi)) %*% crossprod(Zi, residi) 
      Di[i] <- tcrossprod(bi[i]) + solve(solve(D) + crossprod(Zi) * 1/var.0)
      e[i] <- crossprod(residi - Zi %*% bi[i]) + 
        tr(crossprod(Zi) %*% solve(solve(D) + crossprod(Zi) * 1/var.0))
      nis[i] <- nrow(Zi)
    }
    
    # M-step ----
    
    # CM step 1 ---
    # Calculate D.new and var.0.new as in EM 
    var.0.new <- sum(e)/sum(nis)
    D.new <- sum(Di)/n
    
    # CM step 2 ---
    # Estimate new beta
    beta.p1 <- beta.p2 <- list()
    for(i in uids){
      Xi <- getX(data, i); Zi <- getZ(data, i); Yi <- getY(data, i); ni <- nrow(Zi)
      beta.p1[[i]] <- t(Xi) %*% solve(Zi %*% D.new %*% t(Zi) + var.0.new * diag(ni)) %*% Xi
      beta.p2[[i]] <- t(Xi) %*% solve(Zi %*% D.new %*% t(Zi) + var.0.new * diag(ni)) %*% Yi
    }
    beta.new <- solve(Reduce('+', beta.p1)) %*% Reduce('+', beta.p2)
    
    # New parameter vector, and differences
    params.new <- c(t(beta.new), var.0.new, D.new)
    diffs <- abs(params.new-params)
    diff <- max(diffs)
    message("iteration ", iter, " largest difference = ", diff)
    
    # Set new parameters
    var.0 <- var.0.new
    D <- D.new
    beta <- beta.new
    params <- params.new
    iter <- iter + 1

    # Record history if required
    if(history) iter.hist <- rbind(iter.hist, c(iter, t(beta), var.0, D))
    
  }
  t1 <- proc.time()[3]
  message("ECME Complete after ", iter, " iterations, this took ", round(t1-t0,2), " seconds.")
  
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


