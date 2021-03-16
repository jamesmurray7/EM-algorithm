rm(list=ls())
source("~/Documents/PhD/EM-Algorithm/fns/jm/simjoint.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/getW.R")

gammaScore <- function(delta, Ti, bi, gamma, lambda){
  idx <- which(Tis == Ti)
  # Score contribution from subject i
  subj.cont <- delta * bi %*% rbind(1, Ti)
  # Contributed information of those already failed at Ti
  Jpart <- c()
  for(i in 1:idx){
    bTi <- bi %*% rbind(1, Tis[i])
    Jpart[i] <- lambda[i] * bTi %*% exp(bTi * gamma)
  }
  subj.cont - sum(Jpart)
}

gammaInformation <- function(delta, Ti, bi, gamma, lambda,
                             all.times){
  idx <- which(Tis == Ti)
  p1 <- p2 <- Gamma <- c()
  for(i in 1:idx){
    bTi <- bi %*% rbind(1, Tis[i])
    Gamma[i] <- tcrossprod(bTi %*% exp(bTi * gamma))
    p2[i] <- lambda[i]^2 * Gamma[i]/length(which(all.times == Ti))
    p1[i] <- lambda[i] * tcrossprod(bTi) %*% exp(bTi * gamma)
  }
  return(list(p1 = p1, p2 = p2, Gamma = Gamma))
}

x <-simjoint()
dat <- x$joint.data
uids <- unique(dat$id); n <- length(uids)

# Initial parameters
D <- SigmaGen(3, 1); var.e <- 1; gamma <- 1
b <- cbind(rep(1,n), rep(1,n))

# Actual fits, and starting baseline hazard
long.fit <- lmer(Y~1 + (1+time|id), dat)
surv.fit <- coxph(Surv(survtime, status) ~ Y, dat[dat$time == floor(dat$survtime),])
bh <- basehaz(surv.fit)
Tis <- bh$time # ordered times
l0 <- bh$hazard # intial baseline hazard lambda0... 
dat.ord <- dat[order(dat$survtime),]
dat.ord.uniq <- dplyr::distinct(dat.ord, id, survtime, status)

# bi|yi ----
Ebi <- EbibiT <- Igamma <- D.cont <- list()
Sgamma <- e <- mis <- c()
for(i in uids){
  i.dat <- subset(dat.ord, id == i)
  residi <-i.dat$Y - cbind(1, i.dat$time) %*% b[i,]
  W <- getW(i.dat, D, var.e)
  Mui <- W$W21 %*% solve(W$W11) %*% residi
  Vi <- D - W$W21 %*% solve(W$W11) %*% W$W12
  b.mc <- MASS::mvrnorm(100, Mui, Vi)
  
  # Survival stuff
  Ti <- unique(i.dat$survtime); idx <- which(Tis == Ti)
  Deltai <- unique(i.dat$status)
  fti <- c()
  for(j in 1:nrow(b.mc)){
    lft <- Deltai * log(l0[idx]) + Deltai * gamma * (b.mc[j,] %*% rbind(1, Ti)) - sum(l0[1:idx])
    fti[j] <- exp(lft)
  }
  denom <- mean(fti)
  rhs <- fti/denom
  
  # EM
  # E-step ----
  # Random effects
  Ebi[[i]] <- colMeans(b.mc * rhs)
  EbibiT[[i]] <- 1/100 * crossprod(b.mc, b.mc * rhs) 
  # Evaluation of score and information of association parameter.
  Sgamma[i] <- gammaScore(Deltai, Ti, Ebi[[i]], gamma, l0)
  Igamma[[i]] <- gammaInformation(Deltai, Ti, Ebi[[i]], gamma, l0, all.times = dat.ord$survtime)
  
  # M-step ----
  D.cont[[i]] <- tcrossprod(Ebi[[i]])
  e[i] <- crossprod(i.dat$Y - cbind(1, i.dat$time) %*% Ebi[[i]])
  mis[i] <- nrow(i.dat)
}

D.new <- Reduce('+', D.cont)/100
var.e.new <- sum(e)/sum(mis)

Sgamma <- c() # Tomorrow - look at gamma update!
for(i in uids){
  Sgamma <- 
}




