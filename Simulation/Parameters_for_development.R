# Libraries
library(tidyverse)

# parameters to develop sampling function

n.sample = 100
times = seq(0, 10, by = 0.1)

p = 3

# death probability
lambda.t = (0.0005 * times - 1)[-1]
theta.t = 0.05
varphi.t = runif(n = length(times) - 1, min = -0.3, max = 0.3)
xi.t = c(0.01, -0.2, -0.4)


# longitudinal trajectory
zeta.long = (0.0005 * times - 3)
eta.long = runif(n = length(times) - 1, min = -2, max = 1)
beta.long = c(-2, -0.2, 0.2)

# individual frailty --> later
frailty = FALSE

# settings
biv.binary = FALSE
long.trajectory = TRUE

force.mortality = 1
sd.long.trajectory = 1

# cens.poten.rate is not really the censrate
cens.poten.rate = 0




