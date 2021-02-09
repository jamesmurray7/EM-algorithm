#' ###############
#' Simulated code from Van Niekerk et al. (2019)
#' Longitudinal linear predictor L = t^2 + v_i
#' Survival linear predictor S = \beta_s * v_i
#' Toy example...
#' ###############

rm(list = ls())
library(INLA)
num_subj <- 100 # Number of patients
num_times <- 6 # Number of longitudinal observations per person

#Generate time-dependent exponential event times with iid shared random effects term

#Survival part // (some unchanged from van Niekerk)
age <- floor(runif(num_subj,min=18,max=75))
u_i <- rnorm(num_subj, mean=0, sd=1.5) # Subject-specific random effects
coeff <- 0.5 # share parameter (nu?)

lambda <- exp(coeff*u_i+0.01*age)
s_i <- rexp(num_subj, rate=lambda) # Exp. survival times (Do Gompertz etc. later)
c_i <- rep(1, num_subj) #no censoring

# Longitudinal part // (departure from van Niekerk)
t <- 0:(num_times-1)
ID <- 1:num_subj
Y <- rep(t^2, num_subj) + rep(u_i, each = num_times) + rnorm(num_subj * num_times, 0, 0.5)
plot(rep(t, num_subj),Y) #view longitudional data

# Setting up for joint model
# Data comprised as a list in which each variable consists of N_l + N_s elements
N_l <- length(Y)
N_s <- length(s_i)
# > Zeroes for FEs if covariate NOT in that sub-model
# > NAs for REs if covariate NOT in that sub-model

linear.covariate <- list(muJ = as.factor(c(rep(1, N_l), rep(2, N_s))),
                               ltime = c(rep(t, num_subj), rep(0, N_s)),
                               stime = c(rep(0, N_l), s_i),
                               age = c(rep(0, N_l), age)) # Age NOT in longitudinal sub-model!

random.covariate<-list(ID_s=c(rep(NA, N_l), 1:N_s),
                       ID_l=c(rep(ID, each = num_times), rep(NA, N_s)),
                       V1 = c(rep(t, num_subj), rep(NA, N_s))) # Unsure what this V1 is doing (straight from vN)


# Create the outcome - longitudinal and survival.
joint.data <- c(linear.covariate, random.covariate)
Y_joint <- list(c(Y, rep(NA, N_s)),
                inla.surv(time = c(rep(NA, N_l), s_i), event = c(rep(NA,N_l), c_i)))
joint.data$Y <- Y_joint

# Model with INLA
formulaJ <-  Y ~ muJ + age + f(inla.group(V1, n = 5), model="rw2", scale.model = TRUE,
                              hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
           f(ID_l, model = "iid", hyper = list(prec = list(prior="pc.prec", param = c(1, 0.01)))) +
           f(ID_s, copy = "ID_l", fixed = FALSE)

inla.fit <- inla(formulaJ, family = c("gaussian","exponentialsurv"),
               data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE,config=TRUE))
names(inla.fit)

# Getting the longitudinal predictors ----
# Random effects u_i - INLA estimates vs actual
plot(u_i, inla.fit$summary.random$ID_l$mean, pch = 19)
lines(-10:10,-10:10)
# Fitted values - INLA values
plot(rep(t, num_subj), Y, type = "l")
lines(rep(t,num_subj), inla.fit$summary.fitted.values[1:N_l,1], type = "l",
      col = scales::alpha("blue", .5))
legend("topleft", legend = c("INLA", "Actual"), lty = 1, col = c("blue", "black"))

# Estimate of \gamma - Mean of REs in survival/longitudinal
mean(inla.fit$summary.random$ID_s[,2]/inla.fit$summary.random$ID_l[,2])

# Survival part ----
# Estimated survival times
inla.fit$summary.fixed[3,] # True value: 0.1
plot(s_i, log(inla.fit$summary.fitted.values[(N_l+1):(N_l+N_s),1])) # Not sure what INLA actually estimates here?
lines(-10:10,-10:10)

# Longitudinal fit with lme4 ----
library(lme4)
t_long <- rep(t, num_subj)
ID_l <- joint.data$ID_l[!is.na(joint.data$ID_l)]
lmer.fit <- lmer(Y ~ I(t_long^2) + (1|ID_l))
# Compare random effects
plot(x = 1:100, y = ranef(lmer.fit)$ID_l$`(Intercept)`, pch = 19)
points(x = 1:100, y = inla.fit$summary.random$ID_l$mean, pch = 19, col = "blue", cex = .5)

# Joint fit with joineRML ----
library(joineRML)
library(dplyr)
data.long <- data.frame(ID_l = ID_l, time = t_long, Y= Y)
data.surv <- data.frame(ID_s = joint.data$ID_s[!is.na(joint.data$ID_s)],
                        age = joint.data$age[joint.data$age > 0],
                        survtime = s_i,
                        status = 1)

jointdata <- left_join(data.long, data.surv, by = c("ID_l"="ID_s")) %>% 
  filter(time <= survtime)
  
jml.fit <- mjoint(
  formLongFixed = Y ~ I(time^2),
  formLongRandom = ~ 1 | ID_l,
  formSurv = Surv(survtime, status) ~ age,
  data = jointdata, timeVar = "time",
  verbose = T)
)

jml.fit$coefficients
# Estimate of \gamma matches, as does the coefficient for \beta_s (age)

