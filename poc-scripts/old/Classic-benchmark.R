#' ###############
#' Defining functions for classical EM with and wihtout history
#' based on the scripts Classical_EM.R
#' ---
#' This version full of print statements to see where the slow-down is.
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

# Primitive function to measure time difference
pt <- function() as.numeric(proc.time()[3])
pt.diff <- function(old.time){
  new.time <- pt()
  round(new.time - old.time, 2)
}

classic.em <- function(X,Y,Z, init.beta = beta0, init.var.0 = var.0, init.var.1 = var.1,
                       tol = 1e-10, maxiter = 2000){
  message("EM Algorithm ----\nDimensions: ", nrow(Y), " x ", ncol(X)-1, 
          "\nInitial estimates: beta =  ", init.beta, "\nsigma0 = ", sqrt(init.var.0),"; sigma1 = ", sqrt(init.var.1),
          "tol = ", tol, " maximum iterations = ", maxiter)
  
  diff <- 100
  iter <- 0
  var.0 <- init.var.0; var.1 <- init.var.1; beta = init.beta
  
  n <- length(Y); q1 <- ncol(Z)
  
  t0 <- proc.time()
  while(diff > tol & iter <= maxiter){
    p.whileloop <- pt()
    p.ll <- pt()
    L0 <- ll(X, Y, Z,  beta, var.0, var.1) 
    message("L0 evaluation took ", pt.diff(p.ll), "seconds.")
    
    # E step ----
    p.estep <- pt()
    # Covariance stuff
    V <- c(var.1) * Z %*% t(Z) + c(var.0) * diag(n)
    p.vinv <- pt()
    Vinv <- solve(V)
    message("Inverting V took ", pt.diff(p.vinv), " seconds.")
    # X %*% beta and resulting residuals
    Xb <- X %*% beta
    resid <- Y-Xb
    # Calculate expectations uiuiT (call this u.0 (ie eps) and u.1)
    
    p.u0 <- pt()
    u.0 <- c(var.0)^2 * t(resid) %*% Vinv %*% Vinv %*% resid + 
      tr(c(var.0) * diag(n) - c(var.0)^2 * Vinv)
    message("Calculation of u0u0T took ", pt.diff(p.u0), " seconds.")
    p.u1 <- pt()
    u.1 <- c(var.1)^2 * t(resid) %*% Vinv %*% Z %*% t(Z) %*% Vinv %*% resid + 
      tr(c(var.1) * diag(q1) - c(var.1)^2 * t(Z) %*% Vinv %*% Z)
    message("Calculation of u1u1T took ", pt.diff(p.u1), " seconds.")
    
    w <- Xb + c(var.0) * Vinv %*% resid
    
    message("Whole E-step took ", pt.diff(p.estep), " seconds.")
    
    # M step ----
    p.mstep <- pt()
    var.0.new <- u.0/n
    var.1.new <- u.1/q1
    p.beta <- pt()
    beta.new <- solve(t(X) %*% X) %*% t(X) %*% w
    message("Evaluation of new beta took ", pt.diff(p.beta), " seconds.")
    
    message("Whole M-step took ", pt.diff(p.mstep), " seconds.")
    # New ll
    p.l1 <- pt()
    L1 <- ll(X, Y, Z, beta.new, var.0.new, var.1.new)
    message("Evaluation of new log-likelihood took ", pt.diff(p.l1), " seconds.")
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
    message("One loop took ", pt.diff(p.whileloop), " seconds.")
  }
  t1 <- proc.time()
  message("EM converged after ", iter - 1, " iterations, this took ", as.numeric(t1-t0)[3], " seconds.")
  return(list(beta = beta, sigma0 = sqrt(var.0), sigma1 = sqrt(var.1)))
}

# Simulate some data
source("~/Documents/PhD/EM-Algorithm/fns/simlong.R")

# 250 x 5 ----
long.data <- simlong(num_subj = 250, num_times = 5)
X <- long.data$X; Y <- long.data$Y; Z <- long.data$Z
beta <- matrix(long.data$lmer.fit@beta, ncol = 1)

classic.em(X, Y, Z, beta, 1, 1, maxiter = 0)
# L0 evaluation took 0.7seconds.
# Inverting V took 0.58 seconds.
# Calculation of u0u0T took 0.02 seconds.
# Calculation of u1u1T took 0.4 seconds.
# Whole E-step took 1.46 seconds.
# Evaluation of new beta took 0 seconds.
# Whole M-step took 0 seconds.
# Evaluation of new log-likelihood took 0.72 seconds.
# Iteration 0 diff = 291.96213642644
# For beta = [40.5, -10.01, 1.49], sigma01 = [1.48, 0.81]
# One loop took 2.9 seconds.

# 500 x 6 ---- 
long.data <- simlong(num_subj = 500, num_times = 6)
X <- long.data$X; Y <- long.data$Y; Z <- long.data$Z
beta <- matrix(long.data$lmer.fit@beta, ncol = 1)
classic.em(X, Y, Z, beta, 1, 1, maxiter = 0)

# L0 evaluation took 9.35seconds.
# Inverting V took 6.9 seconds.
# Calculation of u0u0T took 0.15 seconds.
# Calculation of u1u1T took 4.73 seconds.
# Whole E-step took 16.11 seconds.
# Evaluation of new beta took 0 seconds.
# Whole E-step took 0 seconds.
# Evaluation of new log-likelihood took 9.12 seconds.
# Iteration 0 diff = 590.98993869196
# For beta = [39.42, -9.95, 1.51], sigma01 = [1.44, 0.77]
# One loop took 34.58 seconds.

