rm(list=ls())
source("~/Documents/PhD/EM-Algorithm/fns/jm/simjoint.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/getW.R")

x <-simjoint()
dat <- x$joint.data
uids <- unique(dat$id); n <- length(uids)

# Initial parameters
D <- SigmaGen(3, 1); B <- rbind(0,0); var.e <- 1; gamma <- 1
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
gh <- statmod::gauss.quad.prob(2, 
                              dist = "normal",
                              mu = 0,
                              sigma = sqrt(0.5))
ab <- gh$nodes
w <- gh$weights * sqrt(pi)
gmat <- matrix(0, 4, 2)
gmat[, 1] <- rep(ab, each = 2)
gmat[, 2] <- rep(ab, 2)
w <- as.vector(w %x% w)

b.new <- matrix(NA, n, 2)
for(i in uids){
  i.dat <- subset(dat.ord, id == i)
  residi <-i.dat$Y - cbind(1, i.dat$time) %*% b[i,]
  W <- getW(i.dat, D, var.e)
  Wzi <- chol((W$W22 - W$W21 %*% solve(W$W11) %*% W$W12) * 2)
  bzi <- matrix(W$W21 %*% solve(W$W11) %*% residi, 4, 2, byrow =  T)
  b.new <- gmat %*% Wzi + bzi
  b.new <- cbind(b.new ^ 2, b.new[,1] * b.new[,2])
}

# logf(Xi, Delta_i | bi)
lfti <- fti <- c()
for(i in uids){
  i.dat <- dat.ord[dat.ord$id == i, ]
  Di <- unique(i.dat$status)
  Ti <- unique(i.dat$survtime)
  haz.idx <- which(Tis == Ti)
  lfti[i] <- Di * log(l0[haz.idx]) + Di * gamma * (b[i,] %*% rbind(1, Ti)) - sum(l0[1:haz.idx])
  fti[i] <- exp(lfti[i])
}

den <- mean(fti)
fti/den

chol((sigma.u - W21 %*% W3) * 2)


