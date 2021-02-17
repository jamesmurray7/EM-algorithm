simlong <- function(num_subj = 250, num_times = 5, 
                    beta = c(40, -10, 1.5), sigma.0 = 1.5, sigma.1 = 0.5){

  N <- num_subj * num_times
  if(!"matrix" %in% class(beta)){
    beta <- matrix(beta, ncol = 1)
  }
  # Covariates 
  x1 <- rbinom(num_subj, 1, 0.5)
  x2 <- floor(runif(num_subj, 50, 75))
  X <- model.matrix(~rep(x1, each = num_times)+rep(x2, each = num_times))
  # Random effects and error term
  u <- rnorm(num_subj, 0, sigma.1)
  eps <- rnorm(N, 0, sigma.0)
  Y <- X %*% beta + rep(u, each = num_times) + eps
  
  long.data <- data.frame(id = rep(1:num_subj, each = num_times),
                          x1 = rep(x1, each = num_times),
                          x2 = rep(x2, each = num_times),
                          Y = Y)
  # Actual fit using lme4::lmer
  actual.fit <- lmer(Y ~ x1 + x2 + (1|id), data = long.data, REML = F)
  Z <- as.matrix(getME(actual.fit, "Z"))
  return(
    list(X = X, Y = Y, Z = Z, long.data = long.data, lmer.fit = actual.fit)
  )
}