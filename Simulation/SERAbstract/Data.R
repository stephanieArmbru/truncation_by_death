### ESTIMATION OF STANDARD JOINT MODEL AND MARGINAL IEE MODEL ###

rm(list = ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(nlme)
library(ggpubr)

library(JM)
library(geepack)


source("Simulation/ShinySimulation/R/Data_generating_mechanism.R")
source("Simulation/ShinySimulation/R/Simulation_analysis_helpers.R")


# Define your custom theme
custom_theme <- theme(
  legend.position = "bottom"
)

# Set the custom theme as the default for all plots
theme_set(theme_bw() + custom_theme)


# Parameters --------------------------------------------------------------

## Global parameters
n <- 10000  # Sample size
times <- seq(0, 10, by = 0.25)  # Discrete Time Partition
N_times <- length(times)

# Mortality Parameters
slope_t_threshold <- -3 # Threshold for slope for add. marker effect
long_t_threshold <- 21 # Threshold for marker value for add. marker effect (scaled by distance to threshold)


lambda_t <- rep(logit(0.00000001), N_times - 1)  # Logit baseline mortality probability
theta_t <- 0  # Global linear marker effect on mortality probability
varphi_t <- c(0, 0, 0, 0,
              rep(0.0, N_times - 5))  # Local linear marker effect on mortality probability
varphi_ReLU_t <- rep(2, N_times - 1) # Local marker effect acc. to threshold
varphi_slope_t <- rep(15, N_times - 1) # Local marker slope effect acc. to threshold
xi_t <- c(0, 0, 0)  # Time-independent covariate effect on mortality probability

# Longitudinal Marker Parameters
long_threshold <- 24  # Threshold for marker change
zeta_long <- sort(pweibull(seq(0, 2.5, length.out = N_times),
                      shape = 2.5, scale = 1) * 2 + 23,
                  decreasing = T)
  # c(-(0.15 * times[1:17])^2 + 25,
  #              rep(min(-(0.2 * times[1:17])^2 + 25), N_times - 17)) # Baseline marker trend
zeta_ReLU_long <- c(rep(-1.5, N_times))  # Threshold trend effect

# stationary
eta_long <- sort(runif(n = N_times - 1, min = 0.7, max = 0.95),
                 decreasing = T)  # Local autoregressive component
slope_threshold <- -1
eta_slope_long <- rep(-1.5, N_times - 1)  # Local autoregressive slope effect


beta_long <- c(0, 0, 0)  # Time-independent covariate effect on marker
sd_long_trajectory <- 0.25  # Standard deviation for longitudinal marker


# Expit -------------------------------------------------------------------
N <- 2000
Time <- 5
expit_df <- data.frame(x = seq(-10, 10, length.out = N),
                       marker = rnorm(N, mean = zeta_long[Time], sd = sd_long_trajectory)) %>%
  mutate(y = expit(x),
         pred = lambda_t[Time] + (max(zeta_long) - marker) * (varphi_t[Time] + theta_t),
         y_example = expit(pred))

expit_df %>% head()

ggplot(data = expit_df,
       aes(x = x,
           y = y)) +
  geom_line() +
  theme_bw()

ggplot(data = expit_df,
       aes(x = marker,
           y = y_example)) +
  geom_line() +
  theme_bw()


# Simulation --------------------------------------------------------------
# Simulate the data
data <- SimulateJointLongData(
  n.sample = n,
  times = times,
  p = length(beta_long),

  slope.t.threshold = slope_t_threshold,
  long.t.threshold = long_t_threshold,
  lambda.t = lambda_t,
  theta.t = theta_t,
  varphi.t = varphi_t,
  varphi.slope.t = varphi_slope_t,
  varphi.ReLU.t = varphi_ReLU_t,
  xi.t = xi_t,


  long.threshold = long_threshold,
  zeta.long = zeta_long,
  zeta.ReLU.long = zeta_ReLU_long,
  eta.long = eta_long,
  slope.threshold = slope_threshold,
  eta.slope.long = eta_slope_long,
  beta.long = beta_long,

  sd.long.trajectory = sd_long_trajectory
)$df.return

plot_baseline_mortality_risk(times, lambda_t)
plot_longitudinal_trend(times, zeta_long)

# patients at risk
plot_patients_at_risk(data)

# calibration check
# plot_calibration(data)
estimate_d_calibration(data)

# mortality probability
plot_mortality_prob(data)
plot_cum_mortality_prob(data)

plot_heterogeneity(data)


plot_immortal_mortal_longitudinal_trajectory(data)
plot_diff_moral_immortal(data)



# Zero scenario -----------------------------------------------------------
# Mortality Parameters
slope_t_threshold <- -Inf # Threshold for slope for add. marker effect
long_t_threshold <- 1 # Threshold for marker value for add. marker effect (scaled by distance to threshold)


lambda_t <- rep(logit(0.00000001), N_times - 1)  # Logit baseline mortality probability
theta_t <- 0  # Global linear marker effect on mortality probability
varphi_t <- c(0, 0, 0, 0,
              rep(0.0, N_times - 5))  # Local linear marker effect on mortality probability
varphi_ReLU_t <- rep(0, N_times - 1) # Local marker effect acc. to threshold
varphi_slope_t <- rep(0, N_times - 1) # Local marker slope effect acc. to threshold
xi_t <- c(0, 0, 0)  # Time-independent covariate effect on mortality probability


data_zero <- SimulateJointLongData(
  n.sample = n,
  times = times,
  p = length(beta_long),

  slope.t.threshold = slope_t_threshold,
  long.t.threshold = long_t_threshold,
  lambda.t = lambda_t,
  theta.t = theta_t,
  varphi.t = varphi_t,
  varphi.slope.t = varphi_slope_t,
  varphi.ReLU.t = varphi_ReLU_t,
  xi.t = xi_t,


  long.threshold = long_threshold,
  zeta.long = zeta_long,
  zeta.ReLU.long = zeta_ReLU_long,
  eta.long = eta_long,
  slope.threshold = slope_threshold,
  eta.slope.long = eta_slope_long,
  beta.long = beta_long,

  sd.long.trajectory = sd_long_trajectory
)$df.return

plot_baseline_mortality_risk(times, lambda_t)
plot_longitudinal_trend(times, zeta_long)

plot_patients_at_risk(data_zero)

plot_heterogeneity(data_zero)

plot_calibration(data_zero)
estimate_d_calibration(data_zero)

plot_cum_mortality_prob(data_zero)

plot_immortal_mortal_longitudinal_trajectory(data_zero)

# prevalence of death correct
(data_zero %>%
    group_by(ID) %>%
    filter(Time == max(Time)) %>%
    pull(Death) %>%
    sum()) / (data_zero$ID %>% unique() %>% length())


# Formatting data ---------------------------------------------------------
# IMMORTAL DATA
# format data for analysis
# due to the autoregressive structure, we have to also regress on the lagged marker value
# lagged marker
data_analysis <- data %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(lagged_Marker = c(0, na.omit(lag(Marker)))) %>%
  ungroup() %>%
  mutate(Time_Dummy = as.factor(Time))


# SURVIVOR DATA
data_analysis_surv <- data_analysis %>%
  filter(Death != 1)


# Regression --------------------------------------------------------------
## ON ALL DATA
# piecewise time trend only
LME_long_dummy <- lme(Marker ~ Time_Dummy,
                      random = ~ 1 | ID,
                      data = data_analysis)

# autoregressive additive component
LME_long <- lme(Marker ~ lagged_Marker + Time_Dummy,
                random = ~ 1 | ID,
                data = data_analysis)

# interaction autoregressive component and time
LME_long_ict <- lm(Marker ~ lagged_Marker * Time_Dummy,
                  data = data_analysis)


# predictions according to LME models
data_analysis$LME_long_pred <- predict(LME_long)
data_analysis$LME_long_dummy_pred <- predict(LME_long_dummy)
data_analysis$LME_long_ict_pred <- predict(LME_long_ict)

# sample IDs for visualization
sampled_ID <- sample(seq(1, length(unique(data_analysis$ID))),
                     5)

# visualize individual trajectories
g_full <- ggplot(data = data_analysis %>%
         filter(ID %in% sampled_ID),
       aes(x = Time,
           y = Marker,
           group = ID,
           color = ID %>% as.factor())) +
  geom_line(aes(linetype = "True trajectory")) +
  geom_line(aes(y = LME_long_pred,
                linetype = "LME prediction")) +
  geom_line(aes(y = LME_long_dummy_pred,
                linetype = "LME Dummy prediction")) +
  geom_line(aes(y = LME_long_ict_pred,
                linetype = "LME Interaction prediction")) +
  labs(color = "", linetype = "",
       title = "Immortal cohort")
g_full


# visualize mean trajectories
data_analysis_LME_mean <- data_analysis %>%
group_by(Time) %>%
  reframe(true_mean = mean(Marker),
          LME_mean = mean(LME_long_pred),
          LME_dummy_mean = mean(LME_long_dummy_pred),
          LME_ict_mean = mean(LME_long_ict_pred))

g_mean_full <- ggplot(data = data_analysis_LME_mean,
       aes(x = Time,
           y = true_mean,
           color = "True mean trajectory")) +
  geom_line() +
  geom_line(aes(y = LME_mean,
                color = "LME mean trajectory")) +
  geom_line(aes(y = LME_dummy_mean,
                color = "LME dummy only mean trajectory")) +
  geom_line(aes(y = LME_ict_mean,
                color = "LME interaction mean trajectory")) +
  labs(color = "")

## ON SURVIVAL DATA ONLY
# piecewise time trend only
LME_long_surv_dummy <- lme(Marker ~ Time_Dummy,
                      random = ~ 1 | ID,
                      data = data_analysis_surv)

# autoregressive additive component
LME_long_surv <- lme(Marker ~ lagged_Marker + Time_Dummy,
                random = ~ 1 | ID,
                data = data_analysis_surv)

# interaction autoregressive component and time
LME_long_surv_ict <- lm(Marker ~ lagged_Marker * Time_Dummy,
                   data = data_analysis_surv)


# predictions according to LME models; extrapolate beyond death
data_analysis$LME_long_surv_dummy_pred <- predict(LME_long_surv_dummy,
                                            newdata = data_analysis)
data_analysis$LME_long_surv_pred <- predict(LME_long_surv,
                                                  newdata = data_analysis)
data_analysis$LME_long_surv_ict_pred <- predict(LME_long_surv_ict,
                                                newdata = data_analysis)


# visualize individual trajectories
g_surv <- ggplot(data = data_analysis %>%
                   filter(ID %in% sampled_ID),
                 aes(x = Time,
                     y = Marker,
                     group = ID,
                     color = ID %>% as.factor())) +
  geom_line(aes(linetype = "True trajectory")) +
  geom_line(aes(y = LME_long_surv_pred,
                linetype = "LME prediction")) +
  geom_line(aes(y = LME_long_surv_dummy_pred,
                linetype = "LME Dummy prediction")) +
  geom_line(aes(y = LME_long_surv_ict_pred,
                linetype = "LME Interaction prediction")) +
  labs(color = "", linetype = "",
       title = "Extrapolation beyond death")
g_surv

# compare subject-specific predicted trajectories
# based on immortal data and survivor data only
ggarrange(g_full, g_surv,
          common.legend = T)


# visualize mean trajectories
data_analysis_LME_mean_surv <- data_analysis %>%
  group_by(Time) %>%
  reframe(true_mean = mean(Marker),
          LME_mean = mean(LME_long_surv_pred),
          LME_dummy_mean = mean(LME_long_surv_dummy_pred),
          LME_ict_mean = mean(LME_long_surv_ict_pred))

g_mean_surv <- ggplot(data = data_analysis_LME_mean_surv,
                      aes(x = Time,
                          y = true_mean,
                          color = "True mean trajectory")) +
  geom_line() +
  geom_line(aes(y = LME_mean,
                color = "LME mean trajectory")) +
  geom_line(aes(y = LME_dummy_mean,
                color = "LME dummy only mean trajectory")) +
  geom_line(aes(y = LME_ict_mean,
                color = "LME interaction mean trajectory")) +
  labs(color = "")

# compare mean trajectory
# based on immortal data and survivor data only
ggarrange(g_mean_full, g_mean_surv,
          common.legend = T)

# Logistic regression -----------------------------------------------------
# with baseline marker value



# landmarking style; time-varying marker as input



LM_death <- glm(Death ~ Marker + Time,
                data = data,
                family = binomial)
LM_death %>% summary()


# Estimating equation ----------------------------------------------------
IEE_mod_long <- geeglm(Marker ~ lagged_Marker + Time,
                            data = data,
                            family = gaussian,
                            id = ID,
                            weights = lagged_Marker + Time,
                            corstr = "independence")

IEE_mod_long_surv <- geeglm(Marker ~ lagged_Marker + Time,
                       data = data,
                       family = gaussian,
                       id = ID,
                       weights = (1 - Death) * (lagged_Marker + Time),
                       corstr = "independence")
IEE_mod_long %>% summary()

IEE_mod_death <- geeglm(Death ~ Marker + Time,
                       data = data,
                       family = binomial,
                       id = ID,
                       weights = Marker + Time,
                       corstr = "independence")
IEE_mod_death %>% summary()


# CART --------------------------------------------------------------------





# Joint model -------------------------------------------------------------
# according to Rizopoulos
# The joint model models (1) the longitudinal marker as a subject specific
# linear mixed model where the random effect follows from a Gaussian distribution
# and (2) the survival probability to depend on the random effect.
# Conditional on the random effect, the marker and the survival probability are
# independent.
# The marginal likelihood, marginalized over the random effect, can be obtained by
# Gauss-Hermite quadrature (as no closed form solution exists).

# We consider a piecewise constant baseline hazard (which aligns with the simulation
# truth) which is estimated based on a B-splines.


# format data to have one row per patient with death time
data.id.dead <- data %>%
  group_by(ID) %>%
  mutate(status = ifelse(any(Death == 1), 1, 0)) %>%
  filter(status == 1) %>%
  filter(Time == min(Time))

data.id.alive <- data %>%
  group_by(ID) %>%
  mutate(status = ifelse(any(Death == 1), 1, 0)) %>%
  filter(status == 0) %>%
  filter(Time == max(Time))

data.id <- rbind(data.id.dead,
                 data.id.alive)

# due to the non-stationary autoregressive structure, we cannot fit a converging linear model
# we render the marker trajectory stationary by lag-1 differencing
data <- data %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(Marker_lagged = c(NA, diff(Marker, 1))) %>%
  na.omit()


# visualize individual death time
ggplot(data = data.id,
       aes(x = Mort_Time)) +
  geom_density() +
  theme_bw()


# jointModel(lmeObject = lme(Marker_lagged ~ Time + X1 + X2 + X3,
#                            random = ~ 1 | ID, data = data),
#            survObject = coxph(Surv(Mort_Time, status) ~ 1 + X1 + X2 + X3,
#                               data = data.id, x = TRUE),
#            timeVar = "Time",
#            parameterization = "value",
#            method = "piecewise-PH-GH",
#            control = list(knots = 15,
#                           verbose = TRUE))


JM_mod <- jointModel(lmeObject = lme(Marker ~ Mort_Time,
                           random = ~ 1 | ID, data = data),
           survObject = coxph(Surv(Mort_Time, status) ~ 1,
                              data = data.id, x = TRUE),
           timeVar = "Mort_Time",
           # parameterization = "value",
           method = "Cox-PH-GH"

           # control = list(knots = 15,
           #                verbose = TRUE)
           )
summary(JM_mod)


# Marginal model ----------------------------------------------------------
# according to Kurland 2005
library(geepack)

data %>% colnames()

# missing weights, but they are multidimensional so I need to use a
# different GEE function / package
IEE_mod <- geeglm(Marker ~ Time,
                  data = data, family = gaussian,
                  id = ID,
                  weights = 1 - Death,
                  corstr = "independence")

summary(IEE_mod)

IEE_coef <- IEE_mod$coefficients

IEE_mortal_trajectory <- data.frame(Time = times) %>%
  mutate(Y_mort = IEE_coef[1] + IEE_coef[2] * Time +
           IEE_coef[3] * 0.5 + IEE_coef[4] * 0.5)



# Visualize ---------------------------------------------------------------

plot_immortal_mortal_longitudinal_trajectory(data) +
  geom_line(data = IEE_mortal_trajectory,
            aes(x = Time,
              y = Y_mort,
              color = "marginal model"))




example_patient <- data %>%
  filter(ID == 15) %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(lag_Mort_Prob = c(0, na.omit(lag(Mort_Prob))),
         Surv = 1 - lag_Mort_Prob,
         cum_Surv = cumprod(Surv),
         marg_Mort_Prob = cum_Surv * Mort_Prob,
         cum_Mort_Prob = cumsum(marg_Mort_Prob)) %>%
  filter(Time != 0)

example_patient_death <- example_patient %>%
  filter(Death == 1) %>%
  filter(Time == min(Time)) %>%
  pull(Time)


g_marker <- plot_immortal_mortal_longitudinal_trajectory(data) +
  geom_line(data = example_patient %>%
              filter(Time <= example_patient_death),
            aes(x = Time,
                y = Marker,
                color = "patient-specific trajectory")) +
  geom_line(data = example_patient %>%
              filter(Time >= example_patient_death),
            aes(x = Time,
                y = Marker,
                color = "patient-specific trajectory"),
            linetype = "dashed") +
  geom_point(data = example_patient %>%
               filter(Death == 1) %>% filter(Time == min(Time)),
             aes(x = Time,
                 y = Marker,
                 color = "Death"),
             shape = 4, size = 3) +
  scale_color_manual(values = c("Mortal cohort" = "#00008B",
                                "Immortal cohort" = "#87CEEB",
                                "Death" = "red",
                                "patient-specific trajectory" = "#A9A9A9"))

g_mortality <- ggplot() +
  geom_line(data = example_patient %>%
              filter(Time <= example_patient_death),
            aes(x = Time,
                y = cum_Mort_Prob,
                color = "patient-specific probability")) +
  geom_line(data = example_patient %>%
              filter(Time >= example_patient_death),
            aes(x = Time,
                y = cum_Mort_Prob,
                color = "patient-specific probability"),
            linetype = "dashed") +
  geom_point(data = example_patient %>%
               filter(Death == 1) %>% filter(Time == min(Time)),
             aes(x = Time,
                 y = cum_Mort_Prob,
                 color = "Death"),
             shape = 4, size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Time",
       y = "Cumulative Death Probability",
       color = "") +
  scale_color_manual(values = c("Death" = "red",
                                "patient-specific probability" = "#A9A9A9"))

ggarrange(g_marker,
          g_mortality,
          nrow = 2, ncol = 1,
          heights = c(2, 1))


g_marker <- plot_immortal_mortal_longitudinal_trajectory(data) +
  geom_line(data = example_patient %>%
              filter(Time <= example_patient_death),
            aes(x = Time,
                y = Marker,
                color = "patient-specific trajectory")) +
  geom_line(data = example_patient %>%
              filter(Time >= example_patient_death),
            aes(x = Time,
                y = Marker,
                color = "patient-specific trajectory"),
            linetype = "dashed") +
  geom_point(data = example_patient %>%
               filter(Death == 1) %>% filter(Time == min(Time)),
             aes(x = Time,
                 y = Marker,
                 color = "Death"),
             shape = 4, size = 3) +
  geom_line(data = example_patient %>%
              filter(Time <= example_patient_death),
            aes(x = Time,
                y = cum_Mort_Prob,
                color = "patient-specific probability")) +
  geom_line(data = example_patient %>%
              filter(Time >= example_patient_death),
            aes(x = Time,
                y = cum_Mort_Prob * 4,
                color = "patient-specific probability"),
            linetype = "dashed") +
  geom_point(data = example_patient %>%
               filter(Death == 1) %>% filter(Time == min(Time)),
             aes(x = Time,
                 y = cum_Mort_Prob * 4,
                 color = "Death"),
             shape = 4, size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Time",
       y = "Cumulative Death Probability",
       color = "") +
  scale_color_manual(values = c("Mortal cohort" = "#00008B",
                                "Immortal cohort" = "#87CEEB",
                                "Death" = "red",
                                "patient-specific trajectory" = "#A9A9A9",
                                "Death" = "red",
                                "patient-specific probability" = "#A9A9A9"))


