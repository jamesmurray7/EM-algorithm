#' #########
#' Doing N simulations of EM under different parameters
#' and using the two different functions written.
#' Dependencies:
#'  fns/classic.R 
#'  fns/simlong.R
#'  fns/Classic-paramtol.R
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
source("./Classic-paramtol.R")
em.p <- em 
# And one on ll
source("./Classic.R")


# Single run - check both working and give expected results.
set.seed(1995)
long <- simlong(num_subj = 100, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)
# Run classic.em
em(X, Y, Z, init.beta, init.var.0 = 1, init.var.1 = 1)
em.p(X,Y,Z, init.beta, init.var.0 = 1, init.var.1 = 1)


# Small sample ------------------------------------------------------------
rm(long)
long <- simlong(num_subj = 100, num_times = 5)
X <- long$X; Y <- long$Y; Z <- long$Z
# Initial estimates of regression coefficients
init.beta <- matrix(lm(Y~x1+x2,data = long$long.data)$coef)

# Don't run this again unless desired (!)
benchmark(
  "ll.tol" = {
    ll.em <- em(X,Y,Z,init.beta,1,1)
  },
  "param.tol" = {
    param.em <- em.p(X,Y,Z,init.beta,1,1)

  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self")
)
# test replications elapsed relative user.self sys.self
# 1    ll.tol           10 105.022    3.953    95.531    8.881
# 2 param.tol           10  26.565    1.000    24.195    2.420

# Check if the estimates are consistent
run1 <- replicate(100, em.p(X,Y,Z,init.beta,1,1, history = T), simplify = F)

# Check items are always the same - no variation run-to-run...
lapply(run1, function(x) print(x$sigma.1))
lapply(run1, function(x) print(x$ll))

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

em.fit <- function(x){
  X <- x$X; Y <- x$Y; Z <- x$Z
  beta0 <- matrix(lm(Y ~ x1 + x2, data = x$long.data)$coefficients,ncol=1)
  fit <- em.p(X, Y, Z, beta0, 1, 1, history = T)
  fit
}

simulated.fits <- lapply(simulated.data, em.fit)


# Plot iteration histories ------------------------------------------------

# Populate a data list
id <- rep(sim.set$id, each = 10)
ids <- data.frame(id, grp.c=1:60)
iter.hists <- list()
for(i in 1:length(simulated.fits)){
  iter.hists[[i]] <- simulated.fits[[i]]$iter.hist
}
iter.hists <- bind_rows(iter.hists)
iter.hists$X4 <- sqrt(iter.hists$X4)
iter.hists$X5 <- sqrt(iter.hists$X5)

# Find where new id begins
iter.hists <- iter.hists %>% 
  mutate(grp = ifelse(iter == 0, 1, 0),
         grp.c=cumsum(grp)) %>% 
  left_join(., ids, "grp.c")

# Plot iteration histories
iter.hists %>% 
  pivot_longer(cols = X1:X5, names_to = "param", values_to = "estimate") %>% 
  ggplot(aes(x = iter, y = estimate, group = grp.c)) + 
  geom_line(alpha = .5) + 
  facet_grid(param~id, scales = "free_y") + 
  theme_bw()
ggsave("../GridIterHist.png")

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
ggsave("../ClassicTimes.png")
