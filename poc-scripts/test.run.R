rm(list=ls())
source("~/Documents/PhD/EM-Algorithm/fns/jm/simjoint.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/getW.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/gammaNR-new.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/lambdaUpdate.R")

x <-simjoint(D=SigmaGen(1,0.5))
dat <- x$joint.data

long.fit <- lmer(Y~1 + (1+time|id), dat)
surv.fit <- coxph(Surv(survtime, status) ~ Y, dat[dat$time == floor(dat$survtime),])
# Survival times
sf <- summary(survfit(surv.fit))
Zif <- cbind(1, sf$time); Zift <- t(Zif)
nf <- length(sf$time)# num fails
bh <- basehaz(surv.fit)
cphd <- coxph.detail(surv.fit)
bh <- dplyr::left_join(bh, data.frame(time = cphd$time, haz.cont = cphd$hazard), "time")
bh[is.na(bh$haz.cont),3] <- 0
l0 <- bh$haz.cont; Tis <- bh$time; rm(cphd)
# Re-ordering data, get ids
dat.ord <- dat[order(dat$survtime),]
uids <- unique(dat.ord$id); n <- length(uids)
# Take unique times
dat.ord.uniq <- dplyr::distinct(dat.ord, id, survtime, status)
raw.Ti <- matrix(dat.ord.uniq$survtime, nc = 1)
Deltas <- matrix(dat.ord.uniq$status, nc = 1)
# Initial parameters
D <- SigmaGen(1, 1); var.e <- 1; gamma <- 0
b <- rbind(rep(1,n), rep(1,n))
# Parameters for quadrature 
gh <- statmod::gauss.quad.prob(2, 
                               dist = "normal",
                               mu = 0,
                               sigma = sqrt(0.5))
ab <- gh$nodes
w <- gh$weights * sqrt(pi)

# Begin loop over subjects ----
# Initialise empty lists
Efti.mat <- matrix(NA, nr = n, nc = nf)
Ebi <- EbibiT <- Ebi2 <- Eexpbu <- Ebexpbu <- Ebbexpbu <- list()
for(i in uids){
  i.dat <- subset(dat, id == i)
  # Data specific //
  Zi <- cbind(1, i.dat$time); Zit <- t(Zi)
  Yi <- matrix(i.dat$Y, nc = 1)
  Di <- unique(i.dat$status)
  Ti <- unique(i.dat$survtime)
  idx <- which(Tis == Ti)
  idx.sf <- which(sf$time <= Ti)
  Zisf <- cbind(1, sf$time[idx.sf]); Zisft <- t(Zisf) # Design matrix (1, t I(t <= Ti))
  Zis <- cbind(1, Tis[1:idx]); Zist <- t(Zis)
  residi <- Yi - Zi %*% b[,i]
  # Mean and variance for GH-Quad transformation
  W <- getW(i.dat, D, var.e)
  Wzi <- chol((W$W22 - W$W21 %*% solve(W$W11) %*% W$W12) * 2)
  bzi <- W$W21 %*% solve(W$W11) %*% residi
  # loop over m absicca an weights
  # Set up empty matrices and vectors
  fti.store <- c()
  #fti.store <- matrix(NA, nr = 4, nc = nf)
  b.abs.store <- matrix(NA, nr = 4, nc = 2)
  bb2.store <- matrix(NA, nr = 4, nc = 3)
  expbu.store <- matrix(NA, nr = 4, nc = nrow(Zis))
  bexpbu.store <- bbexpbu.store <- expbu2.store <- expbu.store
  p <- 1
  for(ii in 1:2){ # Loop over m GPT quadrature points, taking the expectations we need...
    for(jj in 1:2){
      G <- matrix(c(ab[ii],ab[jj]), nc = 1)
      b.abs <- Wzi %*% G + bzi
      ww <- w[ii] * w[jj]
      
      # Collect quadrature points
      b.abs.store[p,] <- b.abs
      bb2.store[p,] <- cbind(t(b.abs^2), t(b.abs)[,1] * t(b.abs)[,2])
      
      # f(T_i,\Delta_i|b_i,Omega)
      expbtD <- 1 # f part, unity if censored
      if(Di == 1) expbtD <- l0[idx] * exp(gamma * cbind(1, Ti) %*% b.abs)
      expSt <- exp(l0[1:idx] %*% -exp(gamma * Zis %*% b.abs)) # Survival function
      fti.store[p] <- expbtD * expSt * ww
      
      # exp(b0+b1u)
      expbu.store[p,] <- exp(gamma * Zis %*% b.abs) # `C` in Pete's code
      # (b0 + b1u) * exp(b0+b1u)
      bexpbu.store[p,1:2] <- colSums(
        tcrossprod(Zis %*% b.abs, exp(gamma * Zis %*% b.abs)) %*% Zis * cbind(l0[1:idx], l0[1:idx])
      )
      # (b0 + b1u) * (b0 + b1u) * exp(b0 + b1u)
      bbexpbu.store[p,] <- crossprod(Zis %*% bb2.store[,1:2], exp(gamma * Zis %*% b.abs) %*% l0[1:idx])
      
      
      p <- p + 1
    }
  }
  
  denom <- sum(fti.store)
  # Take expectations ----
  # E[b], E[b^2], E[b0b1]
  Ebi[[i]] <- crossprod(fti.store, b.abs.store) / denom
  Ebi2[[i]] <- crossprod(fti.store, bb2.store) / denom
  # E[exp(\gamma(b0 + b1u))]
  Eexpbu[[i]] <- crossprod(fti.store, expbu.store) / denom
  # E[(b0 + b1u)exp(\gamma(b0+b1u))]
  Ebexpbu[[i]] <- crossprod(fti.store, bexpbu.store) / denom
  # E[(b0 + b1u)(b0 + b1u)exp(\gamma(b0+b1u))]
  Ebbexpbu[[i]] <- crossprod(fti.store, bbexpbu.store) / denom
}
# Exit loop and M-step ...

Ebi.mat <- matrix(unlist(Ebi), nr = n, nc = 2, byrow = T)
Ebi2.mat <- matrix(unlist(Ebi2), nr = n, nc = 3, byrow = T)# Update for D and var.epsilon ...
D.cont <- list()
mi <- c(); e <- c()
for(i in uids){
  i.dat <- subset(dat.ord, id == i)
  Yi <- matrix(i.dat$Y, nc = 1); Yi2 <- Yi^2
  t <- i.dat$time; ones <- matrix(1, nr = length(t))
  mi[i] <- nrow(Zi)
  # D
  D.cont[[i]] <- tcrossprod(Ebi.mat[i,])
  # e
  Eb0 <- Ebi.mat[i,1]; Eb1 <- Ebi.mat[i,2]
  Eb02 <- Ebi2.mat[i,1]; Eb12 <- Ebi2.mat[i,2]
  Eb0b1 <- Ebi2.mat[i,3]
  r <- Yi2 - Yi %*% Eb0 - Yi %*% Eb1 * t
  e[i] <- sum(r + Eb02 * ones + Eb12 * t + 2 * Eb0b1 * t)
}
D.new <- Reduce('+', D.cont)/n
var.e.new <- (sum(e)/sum(mi))

# Update for lambda - Not currently at all working!
bh.new <- lambdaUpdate(dat.ord, bh, sf, Eexpbu)
bh.new <- cbind(cumsum(bh.new[,1]), bh.new[,2])
# plot against one from coxph
plot(bh[,1]~bh.new[,1],pch=19,cex=.35); lines(-10:10,-10:10,lty=3,col="grey")

gamma.new <- gammaUpdate(dat.ord, sf, gamma, Ebi, Eexpbu, Ebexpbu, Ebbexpbu)
            
# Update parameters, and re-do loop (manually)
D <- D.new; var.e <- var.e.new; gamma <- gamma.new
b <- t(Ebi.mat); l0 <- bh.new[,1]

# Try and put this in a function ! ----------------------------------------
rm(list=setdiff(ls(), "dat.ord"))
source("~/Documents/PhD/EM-Algorithm/fns/jm/simjoint.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/getW.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/gammaNR-new.R")
source("~/Documents/PhD/EM-Algorithm/fns/jm/lambdaUpdate.R")

vech <- function(x) x[lower.tri(x, diag = T)]

jm <- function(data,
               b.init = NULL, D.init = NULL, var.e.init = 1, gamma.init = 0,
               bh, sf,
               tol = 1e-3, maxiter = 200){
  diff <- 100; iter <- 1
  # Data-specific
  uids <- unique(data$id); n <- length(uids)
  Tis <- bh$time # Distinct failure times
  Deltas <- dplyr::distinct(data, id, survtime, status)$status
  # sf
  if(!"summary.survfit"%in%class(sf)){
    if(!"survfitcox"%in%class(sf)) stop("argument `sf` must be survift(coxph), or its summary") else sf <- summary(sf)
  } 
  nf <- length(sf$time) # num fails
  Zif <- cbind(1, sf$time); Zift <- t(Zif)
  
  # Initial parameter estimates //
  # D
  if(is.null(D.init)) D <- diag(1,2) else D <- D.init
  if(!isSymmetric(D) | det(D) == 0) stop("D must be symmetric, positive-definite matrix")
  # b
  if(is.null(b.init)) b <- rbind(rep(1,n), rep(1,n)) else b <- b.init
  if(dim(b)[1] > dim(b)[2]) b <- t(b)
  if(dim(b)[2] != n)  stop("each subject must have initial RE estimate (dim(b)[2] != n)")
  B <- rowMeans(b)
  # l0
  l0 <- bh$hazard
  # gamma
  gamma <- gamma.init
  # var.epsilon
  var.e <- var.e.init
  # Set-up parameter vector
  params <- c(gamma, vech(D), var.e)
  
  # Quadrature - this doesn't change between iterations 
  gh <- statmod::gauss.quad.prob(2, 
                                 dist = "normal",
                                 mu = 0,
                                 sigma = sqrt(0.5))
  ab <- gh$nodes
  w <- gh$weights * sqrt(pi)
  
  while(diff > tol & iter < maxiter){
    # Begin loop over subjects ----
    # Initialise empty lists
    Efti.mat <- matrix(NA, nr = n, nc = nf)
    Ebi <- EbibiT <- Ebi2 <- Eexpbu <- Ebexpbu <- Ebbexpbu <- list()
    for(i in uids){
      i.dat <- subset(data, id == i)
      # Data specific //
      Zi <- cbind(1, i.dat$time); Zit <- t(Zi)
      Yi <- matrix(i.dat$Y, nc = 1)
      Di <- unique(i.dat$status)
      Ti <- unique(i.dat$survtime)
      idx <- which(Tis == Ti)
      idx.sf <- which(sf$time <= Ti)
      Zisf <- cbind(1, sf$time[idx.sf]); Zisft <- t(Zisf) # Design matrix (1, t I(t <= Ti))
      Zis <- cbind(1, Tis[1:idx]); Zist <- t(Zis)
      residi <- Yi - Zi %*% b[,i]
      # Mean and variance for GH-Quad transformation
      W <- getW(i.dat, D, var.e)
      Wzi <- chol((W$W22 - W$W21 %*% solve(W$W11) %*% W$W12) * 2)
      bzi <- W$W21 %*% solve(W$W11) %*% residi
      # loop over m absicca an weights
      # Set up empty matrices and vectors
      fti <- c()
      fti.store <- matrix(NA, nr = 4, nc = nf)
      b.abs.store <- matrix(NA, nr = 4, nc = 2)
      bb2.store <- matrix(NA, nr = 4, nc = 3)
      expbu.store <- matrix(NA, nr = 4, nc = length(idx.sf))
      bexpbu.store <- bbexpbu.store <- expbu2.store <- expbu.store
      p <- 1
      for(ii in 1:2){ # Loop over m GPT quadrature points, taking the expectations we need...
        for(jj in 1:2){
          G <- matrix(c(ab[ii],ab[jj]), nc = 1)
          b.abs <- Wzi %*% G + bzi
          ww <- w[ii] * w[jj]
          # Integral of all survival times' contribution up to this subject's
          temp <- c()
          for(k in 1:idx){ 
            temp[k] <- l0[k] * exp(cbind(1, Tis[k]) %*% b.abs)
          }
          intpart <- sum(temp)
          
          # E[\gamma(b_0 + b_1iu)] used in hazard and 
          # E[(b_0 + b_1u) * exp{\gamma(b_0 + b_1u)}] used in \gamma update.
          if(length(idx.sf) > 0){
            for(k in 1:length(idx.sf)){ 
              expbu.store[p,k] <- exp(gamma * Zisf[k,] %*% b.abs)
              bexpbu.store[p,k] <- Zisf[k,] %*% b.abs %*% exp(gamma * Zisf[k,] %*% b.abs)
              bbexpbu.store[p,k] <- tcrossprod(Zisf[k,] %*% b.abs) %*% exp(gamma * Zisf[k, ] %*% b.abs)
            }
          }
          
          # f(T_i|b_i) - Work out 4x 
          temp <- c()
          for(k in 1:nf){
            temp[k] <- exp(Zif[k, ] %*% b.abs + (Di * cbind(1, Ti) %*% b.abs) - intpart) * ww
          }
          
          fti.store[p,] <- temp
          b.abs.store[p,] <- b.abs
          bb2.store[p,] <- cbind(t(b.abs^2), t(b.abs)[,1] * t(b.abs)[,2]) # bibiT and b0b1 cross-multiplied, used in sigma update
          p <- p + 1
        }
      }
      fti <- rowMeans(fti.store)
      denom <- sum(fti)
      Efti.mat[i,] <- t(fti) %*% fti.store / denom #colMeans(fti.store) / denom
      
      # Expectations ----
      # E[b]
      Ebi[[i]] <- t(fti) %*% b.abs.store / denom
      # E[bi^2] and E[b0b1]
      Ebi2[[i]] <- t(fti) %*% bb2.store / denom
      # E[exp(\gamma(b0 + b1u))]
      Eexpbu[[i]] <- t(fti) %*% expbu.store / denom
      # E[(b0 + b1u)exp(\gamma(b0+b1u))]
      Ebexpbu[[i]] <- t(fti) %*% bexpbu.store / denom
      # E[(b0 + b1u)(b0 + b1u)exp(\gamma(b0+b1u))]
      Ebbexpbu[[i]] <- t(fti) %*% bbexpbu.store / denom
    }
    
    Ebi.mat <- matrix(unlist(Ebi), nr = n, nc = 2, byrow = T)
    Ebi2.mat <- matrix(unlist(Ebi2), nr = n, nc = 3, byrow = T)
    
    # M-step ----
    # Update for D and var.epsilon ...
    D.cont <- list()
    mi <- c(); e <- c()
    for(i in uids){
      i.dat <- subset(data, id == i)
      Yi <- matrix(i.dat$Y, nc = 1); Yi2 <- Yi^2
      t <- i.dat$time; ones <- matrix(1, nr = length(t))
      mi[i] <- nrow(Zi)
      # D
      D.cont[[i]] <- tcrossprod(Ebi.mat[i,])
      # e
      Eb0 <- Ebi.mat[i,1]; Eb1 <- Ebi.mat[i,2]
      Eb02 <- Ebi2.mat[i,1]; Eb12 <- Ebi2.mat[i,2]
      Eb0b1 <- Ebi2.mat[i,3]
      r <- Yi2 - Yi %*% Eb0 - Yi %*% Eb1 * t
      e[i] <- sum(r + Eb02 * ones + Eb12 * t + 2 * Eb0b1 * t)
    }
    D.new <- Reduce('+', D.cont)/n
    var.e.new <- (sum(e)/sum(mi))
    
    bh.new <- lambdaUpdate(data, bh, sf, Eexpbu)
    bh.new <- cbind(cumsum(bh.new[,1]), bh.new[,2])
    
    gamma.new <- gammaUpdate(dat.ord, sf, gamma, Ebi, Eexpbu, Ebexpbu, Ebbexpbu)
    
    params.new <- c(gamma.new, vech(D.new), var.e.new)
    diffs <- abs(params.new - params)
    diff <- max(diffs)
    print(sapply(params, round, 4))
    print(sapply(params.new, round, 4))
    message("Iteration ", iter, " max difference ", diff)
    
    # Set new parameters
    params <- params.new
    D <- D.new
    var.e <- var.e.new
    b <- t(Ebi.mat)
    bh <- bh.new; l0 <- bh.new[,1]
    gamma <- gamma.new
    iter <- iter+1
  }
  
}

surv.fit <- coxph(Surv(survtime, status) ~ Y, dat.ord[dat.ord$time == floor(dat.ord$survtime),])
bh <- basehaz(surv.fit)
sf <- survfit(surv.fit)

jm(data = dat.ord,
   D.init = SigmaGen(1, 0.5), bh = bh, sf = sf)
