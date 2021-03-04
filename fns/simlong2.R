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