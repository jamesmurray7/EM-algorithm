#' ###############
#' Defining functions for Expectation Conditional Maximisation
#' with and without iteration history
#' Based on poc-scripts/ECM and fns/Classical
#' ###############

rm(list = ls())

# Define necessary functions
# Trace
tr <- function(x) sum(diag(x))
# Log likelihood
ll <- function(X, Y, Z, beta, var.0, var.1){ 
  V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(length(Y))
  V.sqrt <- chol(V)
  detV <- sum(log(diag(V.sqrt))) * 2
  # X, beta and residuals
  Xb <- X %*% beta
  resid <- Y - Xb
  temp <- sum(forwardsolve(l = t(V.sqrt), x = resid) ^ 2)
  return(
    -0.5 * temp - 0.5 * detV
  )
}
# Computation of covariance matrix
getV <- function(var.0,var.1,n) c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)


# ECM function 

ecm <- function(X, Y, Z, init.beta, init.var.0, init.var.1,
                tol = 1e-10, maxiter = 2e3, history = F){
  n <- length(Y)
  q1 <- ncol(Z)
  diff <- 10
  iter <- 0
  var.0 <- init.var.0; var.1 <- init.var.1; beta <- init.beta
  
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }

  if(history) iter.hist <- data.frame(iter, t(beta), var.0=sqrt(var.0), var.1=sqrt(var.1), 
                                      ll(X,Y,Z,beta,var.0,var.1))
  
  t0 <- proc.time()[3]
  while(diff > tol & iter <= maxiter){
    l0 <- ll(X,Y,Z, beta, var.0, var.1)
    
    # E step ----
    # Covariance stuff
    
    V <- getV(var.0, var.1, n)
    Vinv <- solve(V)
    V.sqrt <- chol(V)
    # X %*% beta and resulting residuals
    XtX <- crossprod(X)
    Xb <- X %*% beta
    resid <- Y-Xb
    
    # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
    z <- forwardsolve(t(V.sqrt), Z)
    u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
      tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
    u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
      tr(c(var.1) * diag(q1) - c(var.1)^2 * crossprod(z))
    
    # M-step ----
    # CM step 1 ---
    # Estimate \sigma^(t+1)
    var.0.new <- u.0/n
    var.1.new <- u.1/q1
    message("CM step 1 done, new sigma = [", sqrt(var.0.new), ", ", sqrt(var.1.new),"]")
    # CM step 2 ---
    V <- getV(var.0.new, var.1.new, n)
    Vinv <- solve(V)
    Gamma <- X %*% beta + c(var.0.new) * Vinv %*% resid
    # And estimate \beta
    beta.new <- solve(XtX) %*% t(X) %*% Gamma
    message("CM step 2 done, new beta = ", beta.new)
    
    l1 <- ll(X,Y,Z, beta.new, var.0.new, var.1.new)
    
    # New ll
    diff <- abs(l1-l0)
    # Print out
    L0.round <- round(l0, 5); L1.round <- round(l1,5); diff.round <- round(diff, 11)
    message("Iteration ", iter, " diff = ", diff.round)
    # Set new parameters
    var.0 <- var.0.new
    var.1 <- var.1.new
    beta <- beta.new
    iter <- iter + 1
    
    # Record history if required
    if(history) iter.hist <- rbind(iter.hist, c(iter, t(beta), sqrt(var.0), sqrt(var.1), l1))
    
  }
  t1 <- proc.time()[3]
  message("ECM Complete after ", iter, " iterations, this took ", round(t1-t0,2), " seconds.")
  if(history){
    rtn <- list(beta = beta, sigma.0 = sqrt(var.0), sigma.1 = sqrt(var.1), iter.hist = iter.hist,
                time = t1-t0, iter = iter)
  }else{
    rtn <- list(beta = beta, sigma.0 = sqrt(var.0), sigma.1 = sqrt(var.1), time = t1-t0, iter = iter)
  }
  return(rtn)
}


