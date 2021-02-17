#' ###############
#' Defining functions for classical EM with and wihtout history
#' based on the scripts Classical_EM.R
#' ###############

tr <- function(x) sum(diag(x))
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

# Without history ---------------------------------------------------------
classic.em <- function(X,Y,Z, init.beta = beta0, init.var.0 = var.0, init.var.1 = var.1,
                       tol = 1e-10, maxiter = 2000){
  message("EM Algorithm ----\nDimensions: ", nrow(Y), " x ", ncol(X)-1, 
          "\nInitial estimates: beta =  ", beta0, "\nsigma0 = ", sqrt(init.var.0),"; sigma1 = ", sqrt(init.var.1),
          "tol = ", tol, " maximum iterations = ", maxiter)
  
  diff <- 100
  iter <- 0
  var.0 <- init.var.0; var.1 <- init.var.1; beta = init.beta
  
  t0 <- proc.time()
  while(diff > tol & iter <= maxiter){
    L0 <- ll(Y, X, Z,  beta, var.0, var.1) 
    
    # E step ----
    # Covariance stuff
    V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
    Vinv <- solve(V)
    # X %*% beta and resulting residuals
    Xb <- X %*% beta
    resid <- Y-Xb
    # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
    
    u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
      tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
    u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
      tr(c(var.1) * diag(100) - c(var.1)^2 * t(Z) %*% Vinv %*% Z) #  This not working for diag(600) makes me think Z may be wrong.
    
    w <- Xb + c(var.0) * Vinv %*% resid
    
    # M step ----
    var.0.new <- u.0/600
    var.1.new <- u.1/100
    beta.new <- solve(t(X) %*% X) %*% t(X) %*% w
    
    # New ll
    L1 <- ll(Y,X,Z,beta.new,var.0.new, var.1.new)
    diff <- abs(L1-L0)
    # Print out
    L0.round <- round(L0, 5); L1.round <- round(L1,5); diff.round <- round(diff, 11)
    beta.new.round <- round(beta.new,2); b1 <- beta.new.round[1]; b2 <- beta.new.round[2]; b3 <- beta.new.round[3]
    sigma0.print <- round(sqrt(var.0.new), 2); sigma1.print <- round(sqrt(var.1.new),2)
    
    message("Iteration ", iter, " diff = ", diff.round)
    message("For beta = [", b1, ", ", b2, ", ", b3, "], sigma01 = [", sigma0.print, ", ", sigma1.print,"]")
    # Set new parameters
    var.0 <- var.0.new
    var.1 <- var.1.new
    beta <- beta.new
    iter <- iter + 1
  }
  t1 <- proc.time()
  message("EM converged after ", iter - 1, " iterations, this took ", as.numeric(t1-t0)[3], " seconds.")
  return(list(beta = beta, sigma0 = sqrt(var.0), sigma1 = sqrt(var.1)))
}


# With history ------------------------------------------------------------

classicH.em <- function(X,Y,Z, init.beta = beta0, init.var.0 = var.0, init.var.1 = var.1,
                       tol = 1e-10, maxiter = 2000){
  message("EM Algorithm ----\nDimensions: ", nrow(Y), " x ", ncol(X)-1, 
          "\nInitial estimates: beta =  ", beta0, "\nsigma0 = ", sqrt(init.var.0),"; sigma1 = ", sqrt(init.var.1),
          "tol = ", tol, " maximum iterations = ", maxiter)

  
  diff <- 100
  iter <- 0
  var.0 <- init.var.0; var.1 <- init.var.1; beta <- init.beta
  
  iter.hist <- data.frame(iteration = iter, t(beta), sigma0=sqrt(var.0), sigma1=sqrt(var.1),
                          loglik = ll(Y,X,Z, beta = beta, var.0 = var.0, var.1 = var.1))
  
  t0 <- proc.time()
  while(diff > tol & iter <= maxiter){
    L0 <- ll(Y, X, Z,  beta, var.0, var.1) 
    
    # E step ----
    # Covariance stuff
    V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
    Vinv <- solve(V)
    # X %*% beta and resulting residuals
    Xb <- X %*% beta
    resid <- Y-Xb
    # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
    
    u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
      tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
    u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
      tr(c(var.1) * diag(100) - c(var.1)^2 * t(Z) %*% Vinv %*% Z) #  This not working for diag(600) makes me think Z may be wrong.
    
    w <- Xb + c(var.0) * Vinv %*% resid
    
    # M step ----
    var.0.new <- u.0/600
    var.1.new <- u.1/100
    beta.new <- solve(t(X) %*% X) %*% t(X) %*% w
    
    # New ll
    L1 <- ll(Y,X,Z,beta.new,var.0.new, var.1.new)
    diff <- abs(L1-L0)
    # Print out
    L0.round <- round(L0, 5); L1.round <- round(L1,5); diff.round <- round(diff, 11)
    beta.new.round <- round(beta.new,2); b1 <- beta.new.round[1]; b2 <- beta.new.round[2]; b3 <- beta.new.round[3]
    sigma0.print <- round(sqrt(var.0.new), 2); sigma1.print <- round(sqrt(var.1.new),2)
    
    message("Iteration ", iter, " diff = ", diff.round)
    message("For beta = [", b1, ", ", b2, ", ", b3, "], sigma01 = [", sigma0.print, ", ", sigma1.print,"]")
    
    # Store iteration
    iter.hist <- rbind(iter.hist, c(iter, b1, b2, b3, sqrt(var.0.new), sqrt(var.1.new), L1))
    # Set new parameters
    var.0 <- var.0.new
    var.1 <- var.1.new
    beta <- beta.new
    iter <- iter + 1
  }
  t1 <- proc.time()
  message("EM converged after ", iter - 1, " iterations, this took ", as.numeric(t1-t0)[3], " seconds.")
  return(list(beta = beta, sigma0 = sqrt(var.0), sigma1 = sqrt(var.1),
              iter.hist = iter.hist))
}
