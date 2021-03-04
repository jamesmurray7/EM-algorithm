#' ###############
#' Defining functions for classical EM 
#' This operates based on differences on the likelihood
#' ###############

tr <- function(x) sum(diag(x))
ll <- function(X, Y, Z, beta, var.0, var.1){ 
  V <- c(var.1) * tcrossprod(Z) + c(var.0) * diag(length(Y)) # tcrossprod
  V.sqrt <- chol(V)
  detV <- sum(log(diag(V.sqrt))) * 2
  # X, beta and residuals
  Xb <- X %*% beta
  resid <- Y - Xb
  temp <- sum(forwardsolve(l = t(V.sqrt), x = resid) ^ 2)
  return(
    -0.5*length(Y)*log(2*pi) + -0.5 * temp - 0.5 * detV
  )
}

# Define EM function ------------------------------------------------------
em <- function(X,Y,Z, init.beta = beta0, init.var.0 = var.0, init.var.1 = var.1,
                       tol = 1e-5, maxiter = 2000, history = F){
  message("EM Algorithm ----\nDimensions: ", nrow(Y), " x ", ncol(X)-1, 
          "\nInitial estimates: beta =  ", init.beta, "\nsigma0 = ", sqrt(init.var.0),"; sigma1 = ", sqrt(init.var.1),
          "tol = ", tol, " maximum iterations = ", maxiter)
  
  diff <- 100
  iter <- 0
  
  var.0 <- init.var.0; var.1 <- init.var.1; beta = init.beta
  
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }
  
  # Store iteration history if wanted.
  if(history) iter.hist <- data.frame(iter, t(beta), var.0=sqrt(var.0), var.1=sqrt(var.1), 
                                      ll(X,Y,Z,beta,var.0,var.1))
  
  n <- length(Y); q1 <- ncol(Z)
  
  t0 <- proc.time()[3]
  while(diff > tol & iter <= maxiter){
    L0 <- ll(X, Y, Z,  beta, var.0, var.1) 
    
    # E step ----
    # Covariance stuff
    V <- c(var.1) * tcrossprod(Z) + c(var.0) * diag(n)
    Vinv <- solve(V)
    V.sqrt <- chol(V)
    # X %*% beta and resulting residuals
    XtX <- crossprod(X)
    Xb <- X %*% beta
    resid <- Y-Xb
    
    # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
    z <- forwardsolve(t(V.sqrt), Z)
    u.0 <- c(var.0)^2 * crossprod(resid, Vinv) %*% Vinv %*% resid +
      tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
    u.1 <- c(var.1)^2 * crossprod(resid, Vinv) %*% tcrossprod(Z) %*% Vinv %*% resid +
      tr(c(var.1) * diag(q1) - c(var.1)^2 * crossprod(z))
    
    w <- Xb + c(var.0) * Vinv %*% resid
    
    # M step ----
    var.0.new <- u.0/n
    var.1.new <- u.1/q1
    beta.new <- solve(XtX) %*% t(X) %*% w
    
    # New ll
    L1 <- ll(X, Y, Z, beta.new, var.0.new, var.1.new)
    diff <- abs(L1-L0)
    # Print out
    L0.round <- round(L0, 5); L1.round <- round(L1,5); diff.round <- round(diff, 11)
    
    message("Iteration ", iter, " diff = ", diff.round)
    # Set new parameters
    var.0 <- var.0.new
    var.1 <- var.1.new
    beta <- beta.new
    iter <- iter + 1
    
    if(history) iter.hist <- rbind(iter.hist, c(iter, t(beta), sqrt(var.0), sqrt(var.1), l1))
  }
  t1 <- proc.time()[3]
  message("EM converged after ", iter - 1, " iterations, this took ", round(t1-t0, 2), " seconds.")
  
  if(history){
    rtn <- list(beta = beta, sigma.0 = sqrt(var.0), sigma.1 = sqrt(var.1), iter.hist = iter.hist,
                time = t1-t0, iter = iter)
  }else{
    rtn <- list(beta = beta, sigma.0 = sqrt(var.0), sigma.1 = sqrt(var.1), time = t1-t0, iter = iter)
  }
  
  return(rtn)
}
