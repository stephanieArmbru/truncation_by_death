##### EXPLORATION OF PACKAGE LONGITSEMICOMP #####

library(tidyverse)
set.seed(1234)

# Load functions ----------------------------------------------------------
devtools::load_all("Simulation/LongitSemiComp/R")

# for time-independent covariates
?LongitSCparam()

## ARGUMENTS
# longit.data
# A list with entries named risk.NT, risk.T, YNT, YT. See details below.
#
# times
# A vector of increasing times (for example, the interval partition points). This vector is used to construct the B-splines
#
# formula.NT
# A formula of the form ~ x1 + x2 where x1 and x2 are covariates to be used for the non-terminal event probability sub-model.
#
# formula.T
# A formula of the form ~ x1 + x3 where x1 and x3 are covariates to be used for for the non-terminal event probability sub-model.
#
# formula.OR
# A formula of the form ~ x1 + x4 where x1 and x4 are covariates to be used for for the odds ratio sub-model.
#
# formula.inter
# A formula of the form ~ x1 + x4 where x1 and x4 are covariates to be used for the interaction between non-terminal event and covariates sub-model.
#
# data
# A data.frame with the covariates specified in formula.NT, formula.T and formula.OR.
#
# epsOR
# How close it the OR allowed to be to one before assuming it equals to one. Default is 10^(-10)
#
# init
# Initial values for the parameters.
#
# maxit.optim
# For internal use of optim. Default is 50000.

## DATA STRUCTURE
# Each of the matrices in longit.data, risk.NT, risk.T, YNT, YT have a row for
# each unit and a column for each interval.
# Then, risk.NT and risk.T indicate whether the unit is at risk in each interval
# for each of the events, respectively.
# The matrices YNT and YT indicate whether the non-terminal and terminal event,
# respectively, were observed by the end of each interval.
# The function TimesToLongit can be used to obtain this representation of
# semicompeting risks time-to-event data.


# Example -----------------------------------------------------------------
times <- seq(1, 15, 1)
# partition-specific baseline risk
alpha.nt <- LongitSemiComp:::logit(dchisq(times,3, ncp =5)/2 + 0.025)
alpha.t <- LongitSemiComp:::logit(times*(0.075/10)  - 0.0005*(times/20)^2  + 0.05)
alpha.or <- 0.15 - times/10 + 0.75*(times/10)^2 + 0.3*(times/20)^3

plot(x = times, y= exp(alpha.or))
plot(x = times, y= LongitSemiComp:::expit(alpha.nt))
plot(x = times, y= LongitSemiComp:::expit(alpha.t))

# baseline covariate effects
beta.nt <- log(c(0.7, 3))
beta.t <- log(c(0.5, 1))
beta.or <- log(c(1.4, 1))
# effect of Y_NT on Y_T
beta.y <- log(1.4)

my.data <- SimLongitData(n.sample = 2000,
                         times = times,  beta.y,
                         alpha.nt, alpha.t, alpha.or,
                         beta.nt, beta.t, beta.or)

my.data %>% str()
# contains design matrix, YNT, YT, risk.NT and risk.T

# zero / one matrix for patients in row and time partition in column
my.data$YNT %>% head()
my.data$risk.NT %>% head()



longit.data <- my.data[-1]

X <- as.data.frame(my.data[1])


res <- LongitSC(longit.data = longit.data,
                times = times,
                formula.NT =  ~ X.1 + X.2,
                formula.T =  ~ X.1 + X.2,
                formula.OR = ~ X.1 + X.2,
                data = X, epsOR = 10^(-10),
                knots = 5, lambda = 1)

res

res_param <- LongitSCparam(longit.data = longit.data,
                           times = times,
                           formula.NT =  ~ X.1 + X.2,
                           formula.T =  ~ X.1 + X.2,
                           formula.OR = ~ X.1 + X.2,
                           data = X, epsOR = 10^(-10))
res_param

res_param %>% coef()




























