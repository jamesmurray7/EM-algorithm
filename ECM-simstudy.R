#' #########
#' Doing N simulations of ECM under different parameters
#' and seeing if it's (broadly) faster than EM.
#' Dependencies:
#'  fns/Classic-paramtol.R 
#'  fns/simlong.R
#'  fns/ECM-paramtol.R
#' #########

rm(list=ls())
library(lme4)
library(tidyverse)
library(rbenchmark)
# Load scripts from function folder ----
setwd("~/Documents/PhD/EM-Algorithm/fns")
# To simulate longitudinal data
source("simlong.R") 
# Param tolerance first, to rename!
source("./ECM-paramtol.R")
# And one on ll
source("./Classic-paramtol.R")


# Single run - check both working and give expected results.
set.seed(1995)
long <- simlong(num_subj = 100, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)
# Run classic.em
em(X, Y, Z, init.beta, init.var.0 = 1, init.var.1 = 1)
ecm(X, Y, Z, init.beta, init.var.0 = 1, init.var.1 = 1)


# Small sample ------------------------------------------------------------
rm(long)
long <- simlong(num_subj = 100, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)

benchmark(
  "EM" = {
    em.res <- em(X,Y,Z,init.beta,1,1)
  },
  "ECM" = {
    ecm.res <- ecm(X,Y,Z,init.beta,1,1)

  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# test replications elapsed relative user.self sys.self
# 2  ECM           10  28.399    1.000    25.095    3.102
# 1   EM           10  29.416    1.036    25.493    3.544

# Medium sample
rm(long, em.res, ecm.res)
long <- simlong(num_subj = 300, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)
benchmark(
  "EM" = {
    em.res <- em(X,Y,Z,init.beta,1,1)
  },
  "ECM" = {
    ecm.res <- ecm(X,Y,Z,init.beta,1,1)
    
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)

# test replications elapsed relative user.self sys.self
# 2  ECM           10 579.028    1.000   460.115   12.795
# 1   EM           10 597.859    1.033   442.219   14.999

# Simulation study --------------------------------------------------------

# Create function to simulate data and then fit EM to
sim.set <- crossing(
  num_subj = c(100, 250, 500),
  num_times = c(5,10)
) %>% 
  mutate(id = glue::glue("n = {num_subj}, m_i = {num_times}")) %>% 
  group_by(id) %>% 
  nest %>% ungroup

sim.data <- function(df){
  replicate(10, simlong(num_subj = df$num_subj, num_times = df$num_times), simplify = F)
}

sim.set2 <- sim.set %>% 
  mutate(
    simdata = map(data, sim.data)
  )

simulated.data <- flatten(sim.set2$simdata)
source("./ECMp2.R")

ecm.fit <- function(x){
  print(dim(x$long.data))
  X <- x$X; Y <- x$Y; Z <- x$Z
  beta0 <- matrix(lm(Y ~ x1 + x2, data = x$long.data)$coefficients,ncol=1)
  fit <- ecm(X, Y, Z, beta0, 1, 1, history = T)
  fit
}

simulated.fits <- lapply(simulated.data, ecm.fit)
save(simulated.fits, file = "../tempfitsECM.RData")

# Plot iteration histories ------------------------------------------------

# Populate a data list
id <- rep(sim.set$id, each = 10)
ids <- data.frame(id, grp.c=1:60)
iter.hists <- list()
for(i in 1:length(simulated.fits)){
  iter.hists[[i]] <- simulated.fits[[i]]$iter.hist
}
iter.hists <- bind_rows(iter.hists)

# Find where new id begins
iter.hists <- iter.hists %>% 
  mutate(grp = ifelse(iter == 0, 1, 0),
         grp.c=cumsum(grp)) %>% 
  left_join(., ids, "grp.c")

# Plot iteration histories
iter.hists %>% 
  pivot_longer(cols = c(X1:X3, var.0, var.1), names_to = "param", values_to = "estimate") %>% 
  ggplot(aes(x = iter, y = estimate, group = grp.c)) + 
  geom_line(alpha = .5) + 
  facet_grid(param~id, scales = "free_y") + 
  theme_bw()
ggsave("../GridIterHistECM.png")

# Plot time taken ---------------------------------------------------------
times <- data.frame(id = id, time = NA)
for(i in 1:length(simulated.fits)){
  times[i,2] <- simulated.fits[[i]]$time
}

times %>% 
  ggplot(aes(x = id, y = time)) + 
  geom_boxplot() + 
  facet_wrap(~id, scales = "free")+
  theme_bw() + 
  labs(y = "Elapsed time (s)", x = "") + 
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
ggsave("../ECMTimes.png")

# Compare time taken in EM and ECM ----------------------------------------
rm(list = setdiff(ls(), "id"))

setwd("../")
load("EMfits.RData")
em <- simulated.fits
load("tempfitsECM.RData")
ecm <- simulated.fits

times.em <- data.frame(id = id, time = NA)
times.ecm <- times.em

for(i in 1:nrow(times.em)){
  times.em[i, 2] <- em[[i]]$time
  times.ecm[i, 2] <- ecm[[i]]$time
}

times.em$alg <- "EM"; times.ecm$alg <- "ECM"
times <- rbind(times.em, times.ecm)

times %>% 
  ggplot(aes(x = alg, y = time)) + 
  geom_boxplot() + 
  facet_wrap(~id, scales = "free")
