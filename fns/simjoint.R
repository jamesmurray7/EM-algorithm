#' ##########
#' simjoint.R
#' Simulates data under conditions from Wulfsohn '97,
#' with intercept and slope
#' ##########

SigmaGen <- function(int.sd, slp.sd){
  Sigma <- matrix(
    c(int.sd * int.sd, int.sd * slp.sd,
      slp.sd * int.sd, slp.sd * slp.sd),
    2, 2, byrow = T)
  rho <-  matrix(c(1,0.5,0.5,1), 2, 2, byrow = T)
  
  Sigma * rho
}

simjoint <- function(num_subj = 100, num_times = 6,
                     sigma.epsilon = 1.5, D = SigmaGen(3, 1), 
                     gamma = 1, 
                     theta = c(-3, 1), censoring = T, censrate = exp(-3)){
  
  N <- num_subj * num_times
  # Set up time and id
  id <- 1:num_subj
  tau <- num_times - 1
  t <- 0:tau
  # Random effects and error term
  U <- MASS::mvrnorm(num_subj, c(0, 0), Sigma = D)
  b1 <- U[,1]; b2 <- U[,2]
  eps <- rnorm(N, 0, sigma.epsilon)
  
  # Longitudinal //
  
  Y <- rep(b1, each = num_times) + rep(t, num_subj) * rep(b2, each = num_times) + eps
  
  long.data <- data.frame(id = rep(1:num_subj, each = num_times),
                          time = rep(t, num_subj),
                          Y = Y)
  
  # Survival //
  uu <- runif(num_subj)
  num <- (theta[2] + gamma * b2) * log(uu)
  den <- exp(theta[1] + gamma * b1)
  tt <- suppressWarnings((1/theta[2] + gamma * b2) * log(1 - num/den))
  survtime <- rep(tau, num_subj)
  survtime[!is.nan(tt)] <- tt[!is.nan(tt)]
  
  status <- rep(1, num_subj)
  if(censoring){
    censtime <- rexp(num_subj, censrate)
    status[censtime < survtime] <- 0
    survtime[censtime < survtime] <- censtime[censtime < survtime]
  }
  status[survtime >= tau] <- 0
  survtime[survtime >= tau] <- tau
  
  surv.data <- data.frame(
    id = 1:num_subj, survtime = survtime, status
  )
  
  jd <- dplyr::left_join(long.data, surv.data, "id")
  jd2 <- jd[jd$time < jd$survtime,]
  
  message(round(sum(status)/num_subj * 100, 2), "% experienced event")
  
  rtn <- list(long.data = long.data, surv.data = surv.data,
              joint.data = jd2)
  rtn
}
