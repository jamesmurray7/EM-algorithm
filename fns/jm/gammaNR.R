#' #########
#' gammaNR.R // Newton-Raphson update for gamma::
#' Evaluation of score and information in the M-step update for the association parameter, gamma
#' #########

# Consists of evaluation of score and information at current estimates
# And gamma.new = gamma + score/information.

gammaScore <- function(dat, Deltas, 
                       b, gamma, haz, 
                       Tis, # All failure times, including censored times.
                       fail.times){
  S <- c()
  uids <- unique(dat$id)
  for(i in uids){
    i.dat <- subset(dat, id == i)
    Ti <- unique(i.dat$survtime)
    Di <- unique(i.dat$status)
    Zi <- cbind(1, i.dat$time)
    idx <- which(Tis == Ti)
    Tj <- Tis[1:idx]
    EZibi <- mean(Zi %*% b[i,])
    intpart <- c(); p <- 1
    for(u in Tj){
      ii <- which(Tis == u)
      intpart[p] <- haz[ii] * exp(gamma * cbind(1, u) %*% b[i,])
      p <- p+1
    }
    S[i] <- Di * EZibi - sum(intpart)
  }
  sum(S)
}
# gammaScore(dat.ord, Deltas, Ebi.mat, 1, l0, Tis, ft)

gammaInformation <- function(dat, Deltas, 
                             b, gamma, haz, 
                             Tis, # All failure times, including censored times.
                             fail.times){
  # Split information into two parts.
  uids <- unique(dat$id)
  # first part ////
  lhs <- rhs <- c()
  for(i in uids){
    p1 <- c(); p <- 1
    i.dat <- subset(dat, id == i)
    Ti <- unique(i.dat$survtime)
    idx <- which(Tis == Ti)
    Tj <- Tis[1:idx]
    # First leg
    for(u in Tj){
      ii <- which(Tis == u)
      p1[p] <- haz[ii] * tcrossprod(cbind(1, u) %*% b[i,]) %*% exp(cbind(1, u) %*% b[i,])
      p <- p + 1
    }
    lhs[i] <- sum(p1)
  }
  # Second pat ////
  for(i in uids){
    p2 <- c(); p <- 1
    i.dat <- subset(dat, id == i)
    Ti <- unique(i.dat$survtime)
    idx <- which(Tis == Ti)
    Tj <- Tis[1:idx]
    for(u in Tj){
      ii <- which(Tis == u)
      if(u %in% fail.times){
        sv.uids <- uids[which(Tis >= u)]
        haz2 <- haz[ii]^2
        Gammau <- c(); pp <- 1
        for(j in sv.uids){
          Gammau[pp] <- cbind(1, u) %*% b[j, ] %*% exp(gamma * cbind(1, u) %*% b[j, ])
          pp <- pp + 1
        }
        p2[p] <- haz2 * tcrossprod(sum(Gammau))
      }else{
        p2[p] <- 0
      }
      p <- p+1
    }
    rhs[i] <- sum(p2)
  }
  sum(lhs-rhs)
}

gammaUpdate <- function(dat, Deltas, 
                        b, gamma, haz, 
                        Tis, # All failure times, including censored times.
                        fail.times){    
  S <- gammaScore(dat, Deltas, b, gamma, haz, Tis, fail.times)
  I <- gammaInformation(dat, Deltas, b, gamma, haz, Tis, fail.times)
  gamma + S/I
}
