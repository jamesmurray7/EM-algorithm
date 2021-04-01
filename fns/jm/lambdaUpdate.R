#' #########
#' lambdaUpdate.R // Update for the hazard, lambda
#' #########

lambdaUpdate <- function(data, bh, sf, Eexpbu){
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
    bh.idx <- which(bh.new[,1] == u)
    sf.idx <- which(sf$time == u)
    # Numerator (likely = 1 in simulation)
    num <- sum(nev[sf.idx])
    # Denominator, loop through subjects who survive u and sum E[exp(\gamma(b0+b1u))]
    denom <- c(); p <- 1
    for(i in survived.ids){
      denom[p] <- Eexpbu[[i]][bh.idx]
      p <- p + 1
    }
    denom <- sum(denom)
    bh.new[bh.idx,2] <- num/denom
  }
  bh.new  
}

bh.new <- bh
censored <- which(!bh.new[,1]%in%sf$time)

for(u in sf$time){
  survived.ids <- unique(dat.ord[which(dat.ord$survtime >= u), "id"])
  # print(length(survived.ids))
  bh.idx <- which(bh.new[,1]==u)
  num <- sum(nev[])
}


