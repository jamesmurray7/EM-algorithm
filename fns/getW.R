#' ##########
#' getW.R
#' Obtains the W matrices outlined in Wuflsohn 1997,
#' which form the parameters in the MVN distribution of f(b|z,t,omega)
#' ##########

getW <- function(data, D, var.e){
  # Work out dimensionality of data ----
  mi <- nrow(data)
  ti <- data$time
  dimD <- dim(D)[1]
  # W11 ----
  W11 <- matrix(NA, mi, mi)
  for(i in 1:mi){
    for(j in 1:mi){
      W11[i,j] <- cbind(1, ti[i]) %*% D %*% rbind(1,ti[j]) 
    }
  }
  W11 <- W11 + var.e * diag(mi)
  
  # W21 & W12 ----
  W21 <- matrix(NA, nr = dimD, nc = mi)
  for(i in 1:dimD){ # rows
    for(j in 1:mi){ # cols
      W21[i,j] <- D[i, ] %*% c(1, ti[j])
    }
  }
  
  W12 <- t(W21)
  
  # W22 ----
  W22 <- D
  
  return(list(W11=W11, W12 = W12, W21 = W21, W22 = W22))
}


