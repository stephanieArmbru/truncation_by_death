###### SIMULTATIONS ######


# Library -----------------------------------------------------------------
library(tidyverse)

source("Data_generating_mechanism.R")

# Parameter ---------------------------------------------------------------
# according to example in Nevo
times <- seq(1,14,1)
alpha.nt <- logit(times*0.005 + 0.005*(times-2)^2 - (0.0002*(times + 1)^3) + 0.005)
alpha.t <- logit(times*0.0075  + 0.001*(times^2) + 0.4)

lambda.or <- 0.9 + 0.175*times - 0.02*times^2 #+ 0.3*(times/20)^3
lambda.or[times >= 13] <- 0

zeta.long <- 2 + 0.02*(times-2)^2

# plot baseline odds ratio
ggplot(data = data.frame(times = times,
                         lambda.or = lambda.or),
       aes(x = times,
           y = exp(lambda.or))) +
  geom_point() +
  theme_bw() +
  labs(y = "",
       title = "Baseline odds ratio")
# plot baseline hazard ratio
ggplot(data = data.frame(times = times,
                         alpha.nt = alpha.nt),
       aes(x = times,
           y = expit(alpha.nt))) +
  geom_point() +
  theme_bw() +
  labs(y = "",
       title = "Baseline hazard for non-terminal event")
ggplot(data = data.frame(times = times,
                         alpha.t = alpha.t),
       aes(x = times,
           y = expit(alpha.t))) +
  geom_point() +
  theme_bw() +
  labs(y = "",
       title = "Baseline hazard for terminal event")
# plot longitudinal mean
ggplot(data = data.frame(times = times,
                         mean = zeta.long),
       aes(x = times,
           y = mean)) +
  geom_point() +
  theme_bw() +
  labs(y = "",
       title = "Baseline mean longitudinal trajectory for survivors")

# coefficients: longitudinal covariate + covariates
beta.nt <- log(c(1, 0.7, 3))

beta.ntr <- log(1.4)
beta.t <- log(c(0.5, 0.5, 0.1))
nu.or <- log(c(0.2, 1, 1))

eta.long <- c(0.7, 0.2, 1)

sd.long.trajectory <- 1


# Simulation --------------------------------------------------------------
sim_data <- SimulateBivariateLongData(n.sample = 50,
                                      times = times,
                                      alpha.nt = alpha.nt,
                                      alpha.t = alpha.t,
                                      lambda.or = lambda.or,
                                      zeta.long = zeta.long,
                                      beta.nt = beta.nt,
                                      beta.ntr = beta.ntr,
                                      beta.t = beta.t,
                                      nu.or = nu.or,
                                      eta.long = eta.long,
                                      frailty = FALSE,
                                      biv.binary = FALSE,
                                      long.trajectory = FALSE,
                                      force.mortality = 1,
                                      sd.long.trajectory = sd.long.trajectory,
                                      cens.poten.rate = 0.1)
sim_data$df.return %>% head()
sim_data$cens %>% head()
sim_data$Y.traj.all %>% head()

sim_data$Y.traj.all$YT %>% is.na() %>% which()

sim_data$df.return$X.1 %>% range()
sim_data$df.return$X.2 %>% range()
sim_data$Y.traj.all$Ytraj %>% na.omit() %>% range()


# Visualization -----------------------------------------------------------
# mean for all (beyond death)
mean.Y.traj.all <- sim_data$Y.traj.all %>%
  group_by(times) %>%
  summarize(mean_all = mean(na.omit(Ytraj)))
# mean for survivors only
mean.Y.traj.surv <- sim_data$Y.traj.all %>%
  filter(YT == 1) %>%
  group_by(times) %>%
  summarize(mean_surv = mean(na.omit(Ytraj)))

ggplot() +
  geom_point(data = sim_data$Y.traj.all,
             aes(x = times,
                 y = Ytraj,
                 group = ID,
                 color = YT %>% as.factor())
             ) +
  geom_line(data = sim_data$Y.traj.all,
            aes(x = times,
                y = Ytraj,
                group = ID,
                color = YT %>% as.factor())) +
  geom_line(data = mean.Y.traj.all,
            aes(x = times,
                y = mean_all),
            color = "black") +
  geom_line(data = mean.Y.traj.surv,
            aes(x = times,
                y = mean_surv),
            color = "blue") +
  # geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(y = "longitudinal trajectory",
       color = "Death",
       caption = "Black: mean trajectory for all; Blue: mean trajectory for survivor")


ggplot() +
  geom_point(data = sim_data$df.return,
             aes(x = TIME,
                 y = Ytraj,
                 group = ID,
                 shape = YT %>% as.factor())
  ) +
  geom_line(data = sim_data$df.return,
            aes(x = TIME,
                y = Ytraj,
                group = ID,
                color = YT %>% as.factor())) +
  # geom_line(data = mean.Y.traj.all,
  #           aes(x = times,
  #               y = mean_all),
  #           color = "black") +
  # geom_line(data = mean.Y.traj.surv,
  #           aes(x = times,
  #               y = mean_surv),
  #           color = "blue") +
  # geom_smooth(method='lm', formula= y~x) +
  theme_bw() +
  labs(y = "longitudinal trajectory",
       color = "Death",
       caption = "Black: mean trajectory for all; Blue: mean trajectory for survivor")





