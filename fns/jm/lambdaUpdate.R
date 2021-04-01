#' #########
#' lambdaUpdate.R // Update for the hazard, lambda
#' #########

lambdaUpdate <- function(data, bh, sf, Tis, Eexpbu){
  if(class(sf) != "summary.survfit") stop("sf needs to be summary(sf), change later")
  bh.new <- bh 
  # Censored observations have zero point-mass.
  censored <- which(!bh.new[,1]%in%sf$time)
  bh.new[censored,2] <- 0
  # Get number of events from survfit object.
  nev <- sf$n.event
  # Loop over failure times
  for(u in sf$time){
    # Who survives time u?
    survived.ids <- unique(data[which(data$survtime >= u), "id"])
    # index of u
    sf.idx <- which(Tis == u)
    bh.idx <- which(bh.new[,1] == u)
    # Numerator (likely = 1 in simulation)
    num <- sum(nev[sf.idx])
    # Denominator, loop through subjects who survive u and sum E[exp(\gamma(b0+b1u))]
    denom <- c(); p <- 1
    for(i in survived.ids){
      denom[p] <- Eexpbu[[i]][sf.idx]
      p <- p + 1
    }
    denom <- sum(denom)
    bh.new[bh.idx,2] <- num/denom
  }
  bh.new  
}
