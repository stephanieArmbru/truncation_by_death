###### SIMULTATIONS ######


# Library -----------------------------------------------------------------
library(tidyverse)

# Functions
source("Simulation/Data_generating_mechanism.R")
set.seed(123)

# Parameter ---------------------------------------------------------------
# according to example in Nevo
times <- seq(0, 14, 0.5)
# size of sample
n <- 1000


## MORTALITY
lambda.t <- (logit(times*0.0002  + 0.0001*(times^2) + 0.005))[-1]
# positive effect of marker on death probability
theta.t <- 3
varphi.t <- runif(n = length(times) - 1, min = 0, max = 0.5)
xi.t <- c(0.02, -0.04, -0.03)

## LONGITUDINAL MARKER
zeta.long <- 0.005*(times)^2
eta.long <- runif(n = length(times) - 1, min = 0.5, max = 1)
beta.long <- c(-0.02, 0.03, -0.04)


# plot baseline death probability
ggplot(data = data.frame(times = times[-1],
                         lambda.t = lambda.t),
       aes(x = times,
           y = expit(lambda.t))) +
  geom_point() +
  theme_bw() +
  labs(y = "",
       title = "Baseline mortality (partition-specific)")

# plot longitudinal mean
ggplot(data = data.frame(times = times,
                         mean = zeta.long),
       aes(x = times,
           y = mean)) +
  geom_point() +
  theme_bw() +
  labs(y = "",
       title = "Baseline mean longitudinal marker trajectory")

# Simulation --------------------------------------------------------------
sim_dat <- SimulateJointLongData(n.sample = n,
                                 times = times,
                                 p = 3,
                                 lambda.t = lambda.t,
                                 theta.t = theta.t,
                                 varphi.t = varphi.t,
                                 xi.t = xi.t,
                                 zeta.long = zeta.long,
                                 eta.long = eta.long,
                                 beta.long = beta.long,

                                 sd.long.trajectory = 0.001)
# save data set
sim_data <- sim_dat$df.return %>%
  ungroup()

# Visualization -----------------------------------------------------------
# Death
sim_data$Death %>% table()

sim_data %>%
  group_by(ID) %>%
  reframe(dies = any(Death == 1)) %>%
  pull(dies) %>%
  table()

# time of death
sim_data %>%
  group_by(ID) %>%
  filter(Death == 1) %>%
  reframe(death_time = min(Time)) %>%
  pull(death_time) %>%
  table()

# patients at risk
n.risk <- sim_data %>%
  group_by(Time) %>%
  reframe(n_risk = sum(Death == 0)) %>%
  select(Time, n_risk)

ggplot(data = n.risk,
       aes(x = Time,
           y = n_risk)) +
  geom_line() +
  theme_bw() +
  labs(y = "Patients at risk")


# Longitudinal marker
# mean for all (beyond death)
mean.Y.traj.all <- sim_data %>%
  group_by(Time) %>%
  reframe(Time = unique(Time),
          mean_all = mean(na.omit(Marker)))
# mean for survivors only
mean.Y.traj.surv <- sim_data %>%
  filter(Death == 0) %>%
  group_by(Time) %>%
  reframe(Time = unique(Time),
          mean_surv = mean(na.omit(Marker)))

ggplot() +
  geom_line(data = mean.Y.traj.all,
            aes(x = Time,
                y = mean_all,
                color = "Immortal cohort")) +
  geom_line(data = mean.Y.traj.surv,
            aes(x = Time,
                y = mean_surv,
                color = "Mortal cohort")) +
  # geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(y = "Longitudinal Marker",
       color = "")

# Heterogeneity
# new ID
ID_rank_t0 <- sim_data %>%
  filter(Time == 0) %>%
  pull(Marker) %>%
  rank()

sim_data_hetero <- sim_data %>%
  group_by(Time) %>%
  mutate(ID_plot = ID_rank_t0) %>%
  ungroup()


ggplot(data = sim_data_hetero %>%
         filter(ID %in% seq(1, 10)) %>%
         filter(Time %in% c(0, 5, 10, 14)),
       aes(x = ID_plot,
           y = Marker,
           group = Time,
           color = ID %>% as.factor())) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ Time)

ggplot(data = sim_data_hetero %>%
         filter(ID %in% seq(1, 10)),
       aes(x = Time,
           y = Marker,
           color = ID %>% as.factor())) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "")





# OLD ---------------------------------------------------------------------
# alpha.nt <- logit(times*0.005 + 0.005*(times-2)^2 - (0.0002*(times + 1)^3) + 0.005)
# alpha.t <- logit(times*0.0075  + 0.001*(times^2) + 0.4)
#
# lambda.or <- 0.9 + 0.175*times - 0.02*times^2 #+ 0.3*(times/20)^3
# lambda.or[times >= 13] <- 0
#
# zeta.long <- 2 + 0.02*(times-2)^2


