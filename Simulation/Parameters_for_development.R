# Libraries
library(tidyverse)

# parameters to develop sampling function

n.sample = 100
times = seq(1, 100)

# baseline hazard
alpha.nt = 0.005 * times + 1
alpha.t = 0.006 * times + 1
lambda.or = 0.005 * times
# baseline mean
zeta.long = 0.005 * times + 2
# coefficients for predictor
# can depend on time
beta.nt = c(1, 1, 2)
beta.ntr = 2
beta.t = c(1, 1, 2)
nu.or = c(1, 1, 2)
eta.long = c(1, 2, 2)

# individual frailty --> later
frailty = FALSE

# settings
biv.binary = FALSE
long.trajectory = FALSE

force.mortality = 1

# cens.poten.rate is not really the censrate
cens.poten.rate = 0




