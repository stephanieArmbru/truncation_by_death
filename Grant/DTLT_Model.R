####### DISCRETE TIME LONGITUDINAL TRAJECTORY MODEL #######

# Library -----------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(rlist)
library(ggpubr)

# discrete time partition
truncation_time <- 20
granularity <- 10

# sample size
n <- 100

# simulation time period
time_sequence <- seq(0, truncation_time, by = 0.1)


set.seed(1234)


# Functions ---------------------------------------------------------------
mean_trend <- function(t) {
  t * 0.2
}

expit <- function(mu) {
  exp(mu) / (1 + exp(mu))
}

long_traj <- function(x, t, beta, n, has_trend = 0) {
  rnorm(n = n,
        mean = has_trend * mean_trend(t) + x %*% beta,
        sd = 0.1)
}

# Data Simulation ---------------------------------------------------------
time_partition <- seq(1, truncation_time,
                      by = granularity)

# baseline covariates
X <- matrix(c(rep(1, n),
              rnorm(n, mean = 0, sd = 1),
              rnorm(n, mean = 0, sd = 1)),

            nrow = n, ncol = 3, byrow = FALSE)

# Simulation longitudinal trajectory
beta_traj <- c(30, -25, 25)

traj <- list()
j <- 1
for (t in time_sequence) {
  traj[[j]] <- data.frame(id = seq(1, n),
                          time = time_sequence[j],
                          traj = long_traj(X, t, beta_traj, n = n))

  j <- j + 1
}

traj_df <- list.rbind(traj)

ggplot(data = traj_df %>% filter(id %in% sample(seq(1, n), 10)),
       aes(x = time,
           y = traj,
           group = id,
           color = id %>% as.factor())) +
  geom_point() +
  geom_line()



# Simulation event time for death
# Exponential distribution
# global parameters
# beta_death <- c(25, 15, -11, 25)

beta_death <- c(-5.5, 0.2, -0.2, 2)

# local parameters for each time partition subset
# loc_beta <- runif(n = length(time_sequence),
#                   min = -2, max = 2)
loc_beta <- runif(n = length(time_sequence),
                      min = -0.5, max = 0.5)

# data_death <- data.frame(id = seq(1, n),
#                          time = rep(NA, n))
data_death <- data.frame(id = seq(1, n),
                         dead = rep(0, n),
                         time = rep(NA, n))

for (j in seq(1, length(time_sequence) - 1)) {
  long_X <- (traj[[j]]$traj - mean(traj[[j]]$traj)) / sd(traj[[j]]$traj)
  x <- cbind(X, long_X)

  delta_death <- rbinom(n = n, size = 1,
                        prob = expit(x %*% beta_death + loc_beta[j] * long_X))

  # time of death uniformly over time partition interval
  time_death_delta <- runif(n = n,
                            min = time_sequence[j],
                            max = time_sequence[j+1])

  # time_death_delta <-  time_sequence[j] + rexp(n = n,
  #                                              rate = 1 / abs(x %*% beta_death + loc_beta[j] * long_X))

  data_death <- data_death %>%
    mutate(dead = delta_death) %>%
    mutate(time = case_when(is.na(time) & (delta_death == 1) ~ time_death_delta,
                            !is.na(time) ~ time,
                            .default = NA))

    # mutate(prelim_time = time_death_delta) %>%
    # mutate(time = case_when(is.na(time) & (prelim_time <= time_sequence[j + 1]) ~ prelim_time,
    #                         !is.na(time) ~ time,
    #                      .default = NA))
}

data_death <- data_death %>%
  dplyr::select(id, time) %>%
  mutate(time = ifelse(is.na(time), truncation_time + 1, time)) %>%
  group_by(id) %>%
  mutate(time_death = min(time, truncation_time),
         delta_death = time <= truncation_time)


data_death$delta_death %>% table()
data_death$time_death %>% range()

ggplot(data = data_death %>% filter(delta_death == 1),
       aes(x = time_death)) +
  geom_density()


# cumulative hazard
cum_haz <- 0
cum_prob_survival <- 1

# trajectory up to death
traj_death <- list()

# probability
P <- list()
P[[1]] <- data.frame(id = seq(1, n),
                     time = time_sequence[1],
                     delta_death = 0,
                     surv_prob = 1)

# landmarking style
for (j in seq(1, length(time_sequence) - 1)) {
#
#   # longitudinal value at landmarking time point
#   traj_death[[j]] <-  traj[[j]] %>%
#     mutate(traj_trunc = ifelse(data_death$time_death >= time_sequence[j], traj, NA),
#            delta_death = ifelse(data_death$time_death > time_sequence[j], 0, 1))

  # Design matrix
  long_X <- (traj[[j]]$traj - mean(traj[[j]]$traj)) / sd(traj[[j]]$traj)
  x <- cbind(X, long_X)

  # probability for death
  prob_death <- 1 - exp(-(1 / abs(x %*% beta_death + loc_beta[j] * long_X)) %*% (time_sequence[t + 1] - time_sequence[t]))

  cum_prob_survival <- cum_prob_survival * (1 - prob_death)


  P[[j + 1]] <- data.frame(id = seq(1, n),
                       time = time_sequence[j + 1],
                       delta_death = ifelse(data_death$time_death > time_sequence[j], 0, 1),
                       surv_prob = cum_prob_survival)
}

# format to data frame
P_df <- list.rbind(P)

res_df <- inner_join(traj_df,
                     P_df,
                     by = c("id", "time"))


# Occurrence --------------------------------------------------------------
# Create a plot with vertical lines and points
ggplot(data_death[sample(seq(1, n), 10), ],
       aes(x = time,
           y = id)) +
  geom_segment(aes(x = 0,
                   xend = pmin(time, truncation_time),
                   y = id,
                   yend = id),
               color = "#696969", size = 1) +
  geom_point(size = 3) +
  xlim(c(0, truncation_time)) +
  labs(x = "Time",
       y = "Patient",
       color = "") +
  theme_bw() +
  theme(legend.position = "bottom")


# Visualize trajectory and mortality --------------------------------------
# identify good example patients
P_df %>%
  group_by(id) %>%
  summarize(min(surv_prob %>% na.omit()))
P_df %>% filter(delta_death == 1)
P_df %>% filter(delta_death == 0 & time == truncation_time)


example_ind <- c(13, 74, 44, 99)

# scaling for visualization
scale_factor <- traj_df$traj %>%
  na.omit() %>%
  max() %>%
  ceiling()

g1 <- ggplot() +
  geom_smooth(data = traj_df %>%
                filter(id %in% example_ind) %>%
                na.omit(),
              aes(x = time,
                  y = traj / scale_factor,
                  color = "Trajectory"),
              se = F) +
  geom_smooth(data = P_df %>%
                filter(id %in% example_ind),
              aes(x = time,
                  y = 1 - surv_prob,
                  color = "Mortality Probability"),
              se = F) +
  # geom_vline(data = data_death[example_ind, ],
  #            mapping = aes(xintercept = ifelse(time <= truncation_time, time, NA)),
  #            color = "black",
  #            size = 1, linetype = "dashed") +
  scale_y_continuous(
    name = "Mortality Probability",
    sec.axis = sec_axis(~ . * scale_factor, name = "Longitudinal Trajectory"),
    # limits = c(0, 1)
  ) +
  scale_color_manual(values = c("Trajectory" = "forestgreen",
                                "Mortality Probability" = "#DC267F")) +
  labs(x = "Time",
       y = "Patient",
       color = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_grid(~ id)


g2 <- ggplot(data_death[example_ind, ],
             aes(x = time,
                 y = 1)) +
  geom_segment(aes(x = 0,
                   xend = pmin(time, truncation_time),
                   y = 1,
                   yend = 1),
               color = "#696969", size = 1) +
  geom_point(size = 3) +
  xlim(c(0, truncation_time)) +
  labs(y = "",
       x = "Time") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        axis.title.y = element_blank(),   # Remove y-axis title
        axis.title.y.right = element_blank() # Remove secondary y-axis title
  ) +
  facet_grid(~ id)


ggarrange(g1, g2,
          ncol = 1, nrow = 2,
          common.legend = TRUE,
          legend = "bottom",
          heights = c(3, 1),
          align = "v")

ggsave("Grant/joint_traj.pdf",
       height = 10, width = 20, unit = "cm")

g1
ggsave("Grant/joint_traj_only.pdf",
       height = 10, width = 20, unit = "cm")


# Population average for mortal and immortal cohort -----------------------

res_df %>%
  filter(time == 0) %>%
  summarize(traj = mean(traj),
            traj_delta = mean(traj * (1 - delta_death)))

res_df %>%
  filter(time == 0 & traj > 20 & traj < 30)


traj_mortal_immortal <- res_df %>%
  group_by(time) %>%
  summarize(mean_mortal = mean(traj * (1 - delta_death)),
            mean_immortal = mean(traj))

g_mean <- ggplot(data = traj_mortal_immortal) +
  geom_line(aes(x = time,
                y = mean_mortal,
                color = "Mortal cohort"),
            linewidth = 1) +
  geom_line(aes(x = time,
                y = mean_immortal,
                color = "Immortal cohort"),
            linewidth = 1) +
  geom_line(data = res_df %>%
              filter(id == 6),
            aes(x = time,
            y = traj),
            color = "gray",
            linewidth = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(y = "Longitudinal Trajectory",
       x = "Time",
       color = "") +
  scale_color_manual(values = c("Immortal cohort" = "#66C266",
                                "Mortal cohort" = "#154F15"))

g_mean


# Heterogeneity -----------------------------------------------------------
# according to PIPER-ICD grant

# order according to longitudinal marker at time 0
ordering_t0 <- res_df %>%
  filter(time == 0) %>%
  pull(traj) %>%
  rank()

res_df <- res_df %>%
  group_by(time) %>%
  mutate(id_ordered = ordering_t0)

# at time point 0
res_t0 <- res_df %>%
  filter(time == 0) %>%
  mutate(plot_color = as.character(delta_death))

# at time point 10
res_t1 <- res_df %>%
  filter(time == 10) %>%
  mutate(plot_color = as.character(delta_death))

# at time point 20
res_t2 <- res_df %>%
  filter(time == 19) %>%
  mutate(plot_color = as.character(delta_death))

# time 0
g_t0_mortality <- ggplot() +
  geom_line(data = res_t0,
           aes(x = id_ordered,
               y = 1 - surv_prob,
               color = plot_color)) +
  scale_color_manual(values = c("0" = "#DC267F",
                                "1" = "#CACACA"),
                     labels = c("0" = "Alive",
                                "1" = "Dead")) +
  labs(x = "Patient ID", y = "Death Probability",
       color = "", fill = "",
       title = expression(tau[k] *"=0")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylim(c(0, 1))

g_t0_marker <- ggplot() +
  geom_point(data = res_t0,
            aes(x = id_ordered,
                y = traj,
                color = plot_color),
            alpha = 0.7) +
  scale_color_manual(values = c("0" = "forestgreen",
                                "1" = "#CACACA"),
                     labels = c("0" = "Alive",
                                "1" = "Dead")) +
  labs(x = "Patient ID", y = "Longitudinal marker",
       color = "", fill = ""
       ) +
  theme_bw() +
  theme(legend.position = "bottom")

# time 10
g_t1_mortality <- ggplot() +
  geom_col(data = res_t1,
            aes(x = id_ordered,
                y = 1 - surv_prob,
                fill = plot_color)) +
  scale_fill_manual(values = c("0" = "#DC267F",
                                "1" = "#CACACA"),
                     labels = c("0" = "Alive",
                                "1" = "Dead")) +
  labs(x = "Patient ID", y = "Death Probability",
       color = "", fill = "",
       title = expression(tau[k] *"=10")) +
  theme_bw() +
  theme(legend.position = "bottom")

g_t1_marker <- ggplot() +
  geom_point(data = res_t1,
             aes(x = id_ordered,
                 y = traj,
                 color = plot_color),
             alpha = 0.7) +
  scale_color_manual(values = c("0" = "forestgreen",
                                "1" = "#CACACA"),
                     labels = c("0" = "Alive",
                                "1" = "Dead")) +
  labs(x = "Patient ID", y = "Longitudinal marker",
       color = "", fill = ""
       ) +
  theme_bw() +
  theme(legend.position = "bottom")


# time 19
g_t2_mortality <- ggplot() +
  geom_col(data = res_t2,
            aes(x = id_ordered,
                y = 1 - surv_prob,
                fill = plot_color)) +
  scale_fill_manual(values = c("0" = "#DC267F",
                                "1" = "#CACACA"),
                     labels = c("0" = "Alive",
                                "1" = "Dead")) +
  labs(x = "Patient ID", y = "Death Probability",
       color = "", fill = "",
       title = expression(tau[k] *"=19")) +
  theme_bw() +
  theme(legend.position = "bottom")

g_t2_marker <- ggplot() +
  geom_point(data = res_t2,
             aes(x = id_ordered,
                 y = traj,
                 color = plot_color),
             alpha = 0.7) +
  scale_color_manual(values = c("0" = "forestgreen",
                                "1" = "#CACACA"),
                     labels = c("0" = "Alive",
                                "1" = "Dead")) +
  labs(x = "Patient ID", y = "Longitudinal marker",
       color = "", fill = ""
       ) +
  theme_bw() +
  theme(legend.position = "bottom")



g <- ggarrange(g_t0_mortality,
          g_t1_mortality,
          g_t2_mortality,

          g_t0_marker,
          g_t1_marker,
          g_t2_marker,
          ncol = 3, nrow = 2,
          heights = c(1, 0.9))

g



# Save --------------------------------------------------------------------
ggsave(plot = g,
       filename = "Grant/joint_traj_at_landmark_times_no_mean.pdf",
       height = 13,
       width = 18,
       unit = "cm")

ggsave(filename = "Grant/mean_traj.pdf",
       plot = g_mean,
       height = 7,
       width = 10, unit = "cm")

ggarrange(g_time,
          g_mean + labs(title = "Population averaged"),
          # align = "h",
          ncol = 2,
          widths = c(3, 1))
ggsave("Grant/joint_traj_at_landmark_times.pdf",
       height = 10,
       width = 30,
       unit = "cm")

