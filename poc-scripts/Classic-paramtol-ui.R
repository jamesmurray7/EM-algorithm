#' ###############
#' Single run to see if Classic EM produces Uis under two approaches.
#' ###############

rm(list=ls())
setwd("~/Documents/PhD/EM-Algorithm/fns")
library(lme4)
source("simlong.R")
source("Classic-paramtol.R")
# Post-processing for U-i.
# Gamma_i = (N_i/(sigma.epsilon^2) + 1/(sigma.u^2))^-1
# u_i = Gamma_i * Z_i * (Y - Xb)/(sigma.epsilon^2)
x <- simlong(100,5)
X <- x$X; Y <- x$Y;Z <- x$Z
beta0 <- matrix(lm(Y~x1+x2, data = x$long.data)$coef,nc=1)


em.fit <- em(X,Y,Z,beta0,1,1)
sigma.e <- em.fit$sigma.0; sigma.u <- em.fit$sigma.1
var.e <- sigma.e^2; var.u <- sigma.u^2
betas <- em.fit$beta
dat <- x$long.data

# bi = ZiSigmae^-1Sigmau * (Yi- Xi73)

ui <- c() # This way from Laird
for(i in unique(dat$id)){
  i.data <- subset(dat, id == i)
  Yi <- i.data$Y; Xi <- cbind(1, i.data$x1, i.data$x2)
  Zi <- matrix(rep(1,5),ncol = 1)
  residi <- Yi - Xi %*% betas
  Sigma <- Zi %*% var.u %*% t(Zi) + c(var.e) * diag(length(Yi))
  Sigma.inv <- solve(Sigma)
  ui[i] <- t(Zi %*% var.u) %*% Sigma.inv %*% residi
}

ui2 <- c() # This from the website and given by Pete in meeting.
for(i in unique(dat$id)){
  i.data <- subset(dat, id == i)
  Yi <- i.data$Y; Xi <- cbind(1, i.data$x1, i.data$x2)
  Zi <- matrix(rep(1,5),ncol = 1)
  residi <- Yi - Xi %*% betas
  Gamma_i <- solve(length(Yi)/var.e + 1/var.u)
  test <- (Gamma_i %*% t(Zi) %*% residi)/var.e
  ui2[i] <- test
}

# Either method produces good estimates, it seems.
plot(ui, ranef(x$lmer.fit)$id$`(Intercept)`, pch = 19, cex = .5)
lines(-10:10,-10:10)
plot(ui2, ranef(x$lmer.fit)$id$`(Intercept)`, pch = 19, cex = .5)
lines(-10:10,-10:10)
