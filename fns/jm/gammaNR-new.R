#' #########
#' gammaNR.R // Newton-Raphson update for gamma::
#' Evaluation of score and information in the M-step update for the association parameter, gamma
#' #########

# Consists of evaluation of score and information at current estimates
# And gamma.new = gamma + score/information.

gammaUpdate <- function(data, sf, l0, Tis,
                        gamma,
                        Ebi, Eexpbu, Ebexpbu, Ebbexpbu){
  fail.times <- sf$time
  filter.data <- subset(data, survtime %in% fail.times) # Only interested in failure times due to delta Pre-multiply
  uids <- unique(filter.data$id)
  Score.cont <- Information.cont <- c()
  p <- 1
  for(i in uids){
    Ti <- unique(filter.data[filter.data$id == i, "survtime"])
    ZiTi <- cbind(1, Ti)
    Ti.idx <- which(Tis == Ti)
    Eexpbuh <- Eexpbu[[i]] * l0[1:Ti.idx]
    # Score and information
    Score.cont[p] <- ZiTi %*% t(Ebi[[i]]) - sum(Ebexpbu[[i]])/sum(Eexpbuh)
    Information.cont[p] <- (sum(Ebbexpbu[[i]]) * sum(Eexpbuh) - sum(Ebexpbu[[i]])^2)/(sum(Eexpbuh)^2)
    p <- p + 1
  }
  return(list(gamma.new = gamma + sum(Score.cont)/sum(Information.cont),
              S = Score.cont, I = Information.cont))
  # gamma + sum(Score.cont)/sum(Information.cont)
}
