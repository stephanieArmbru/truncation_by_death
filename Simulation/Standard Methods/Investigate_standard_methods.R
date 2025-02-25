### ESTIMATION OF STANDARD JOINT MODEL AND MARGINAL IEE MODEL ###

rm(list = ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(ggpubr)

# packages for generalized linear mixed models
library(nlme)
library(lme4)

# packages for GEE
library(gee)
library(geepack)

# package for joint models
library(JM)



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
slope_t_threshold <- -5 # Threshold for slope for add. marker effect
long_t_threshold <- 15 # Threshold for marker value for add. marker effect (scaled by distance to threshold)


lambda_t <- rep(logit(0.02), N_times - 1)  # Logit baseline mortality probability
theta_t <- 0.02  # Global linear marker effect on mortality probability
varphi_t <- rep(0.05, N_times - 1)  # Local linear marker effect on mortality probability
varphi_ReLU_t <- rep(1.5, N_times - 1) # Local marker effect acc. to threshold
varphi_slope_t <- rep(2, N_times - 1) # Local marker slope effect acc. to threshold
xi_t <- c(-0.1, 0.2, -0.3)  # Time-independent covariate effect on mortality probability

# Longitudinal Marker Parameters
long_threshold <- 24  # Threshold for marker change
zeta_long <- sort(25 - 2 / (1 + exp(-2.5 * (seq(0, 2.5, length.out = 41) - 2))),
                  decreasing = T)
  # c(-(0.15 * times[1:17])^2 + 25,
  #              rep(min(-(0.2 * times[1:17])^2 + 25), N_times - 17)) # Baseline marker trend
zeta_ReLU_long <- c(rep(-1.5, N_times))  # Threshold trend effect

# stationary
eta_long <- rep(0.9, 40)
  # runif(n = N_times - 1, min = 0.85, max = 0.95) # Local autoregressive component
slope_threshold <- -1
eta_slope_long <- rep(-1.5, N_times - 1)  # Local autoregressive slope effect


beta_long <- c(-0.2, 0.3, 0.5)  # Time-independent covariate effect on marker
sd_long_trajectory <- 0.25  # Standard deviation for longitudinal marker


# Expit -------------------------------------------------------------------
N <- 2000
Time <- 5
expit_df <- data.frame(x = seq(-10, 10, length.out = N),
                       marker = rnorm(N, mean = zeta_long[Time], sd = sd_long_trajectory)) %>%
  mutate(y = expit(x - 2),
         pred = lambda_t[Time] - ((theta_t +
                      varphi_t[Time - 1]) * (zeta_long[Time] - marker )) +
           varphi_ReLU_t[Time] + varphi_slope_t[Time],
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
res <- SimulateJointLongData(
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
)

# simulated data for immortal cohort
data <- res$df.return.IM

plot_baseline_mortality_risk(times, lambda_t)

# baseline marker trajectory
lm <- lm(zeta_long ~ times)
plot_longitudinal_trend(times, zeta_long) +
  geom_abline(aes(intercept = lm$coefficients[1],
                  slope = lm$coefficients[2]))




# patients at risk
plot_patients_at_risk(data)

# calibration check
# plot_calibration(data)
D_calibration(data)

# mortality probability
plot_mortality_prob(data)
plot_cum_mortality_prob(data)

# longitudinal marker trajectory
plot_heterogeneity(data)


plot_immortal_mortal_longitudinal_trajectory(data)
plot_diff_moral_immortal(data)

# # Zero scenario -----------------------------------------------------------
# # Mortality Parameters
# slope_t_threshold <- -Inf # Threshold for slope for add. marker effect
# long_t_threshold <- 1 # Threshold for marker value for add. marker effect (scaled by distance to threshold)
#
#
# lambda_t <- rep(logit(0.00000001), N_times - 1)  # Logit baseline mortality probability
# theta_t <- 0  # Global linear marker effect on mortality probability
# varphi_t <- c(0, 0, 0, 0,
#               rep(0.0, N_times - 5))  # Local linear marker effect on mortality probability
# varphi_ReLU_t <- rep(0, N_times - 1) # Local marker effect acc. to threshold
# varphi_slope_t <- rep(0, N_times - 1) # Local marker slope effect acc. to threshold
# xi_t <- c(0, 0, 0)  # Time-independent covariate effect on mortality probability
#
#
# data_zero <- SimulateJointLongData(
#   n.sample = n,
#   times = times,
#   p = length(beta_long),
#
#   slope.t.threshold = slope_t_threshold,
#   long.t.threshold = long_t_threshold,
#   lambda.t = lambda_t,
#   theta.t = theta_t,
#   varphi.t = varphi_t,
#   varphi.slope.t = varphi_slope_t,
#   varphi.ReLU.t = varphi_ReLU_t,
#   xi.t = xi_t,
#
#
#   long.threshold = long_threshold,
#   zeta.long = zeta_long,
#   zeta.ReLU.long = zeta_ReLU_long,
#   eta.long = eta_long,
#   slope.threshold = slope_threshold,
#   eta.slope.long = eta_slope_long,
#   beta.long = beta_long,
#
#   sd.long.trajectory = sd_long_trajectory
# )$df.return
#
# plot_baseline_mortality_risk(times, lambda_t)
# plot_longitudinal_trend(times, zeta_long)
#
# plot_patients_at_risk(data_zero)
#
# plot_heterogeneity(data_zero)
#
# plot_calibration(data_zero)
# estimate_d_calibration(data_zero)
#
# plot_cum_mortality_prob(data_zero)
#
# plot_immortal_mortal_longitudinal_trajectory(data_zero)
#
# # prevalence of death correct
# (data_zero %>%
#     group_by(ID) %>%
#     filter(Time == max(Time)) %>%
#     pull(Death) %>%
#     sum()) / (data_zero$ID %>% unique() %>% length())
#
#


# Formatting data ---------------------------------------------------------
# IMMORTAL DATA
# format data for analysis
# due to the autoregressive structure, we have to also regress on the lagged marker value
# lagged marker
# arrange data for clustering in GEE methods
data_train <- data %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(lagged_Marker = c(0, na.omit(lag(Marker)))) %>%
  ungroup() %>%
  mutate(Time_Dummy = as.factor(Time)) %>%
  arrange(ID)

# SURVIVOR DATA
data_train_surv <- data_train %>%
  group_by(ID) %>%
  mutate(obs_post_death = cumsum(Death)) %>%
  filter(obs_post_death <= 1)


# sample IDs for visualization
sampled_ID <- sample(seq(1, length(unique(data_train$ID))),
                     5)

# Validation data set -----------------------------------------------------
# validation data set
n_val <- 10000
data_val <- SimulateJointLongData(
  n.sample = n_val,
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
)$df.return %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(lagged_Marker = c(0, na.omit(lag(Marker)))) %>%
  ungroup() %>%
  mutate(Time_Dummy = as.factor(Time))

plot_heterogeneity(data_val)
plot_immortal_mortal_longitudinal_trajectory(data_val)

plot_cum_mortality_prob(data_val)


# Mortal and immortal mean
data_val_mean <- MorImmorMean(data_val)

data_val_mean_ID <- do.call(rbind, replicate(sampled_ID %>% length(),
                         data_val_mean,
                         simplify = FALSE)) %>%
  arrange(Time) %>%
  mutate(ID = rep(sampled_ID, length(times) - 1))


# 1. Longitudinal Marker --------------------------------------------------

# 1.1 Immortal cohort -----------------------------------------------------
# piecewise time trend only
LMM_long_dummy <- lme(Marker ~ Time_Dummy + X1 + X2 + X3,
                      random = ~ 1 | ID,
                      data = data_train)

# autoregressive additive component (landmarking style)
LMM_long <- lme(Marker ~ lagged_Marker + Time_Dummy + X1 + X2 + X3,
                random = ~ 1 | ID,
                data = data_train)

# predictions according to LMM models
data_val$LMM_long_dummy_pred <- predict(LMM_long_dummy,
                                        newdata = data_val,
                                        level = 0)
data_val$LMM_long_pred <- predict(LMM_long,
                                  newdata = data_val,
                                  level = 0)

data_val$LMM_long_bl_pred <- predict_baseline(LMM_long$coefficients$fixed,
                                              newdata = data_val)

# landmarking style; independent error between
GEE_lagged_long <- geeglm(Marker ~ lagged_Marker + Time_Dummy + X1 + X2 + X3,
                          id = ID,
                          data = data_train,
                          family = gaussian,
                          corstr = "independence")

data_val$GEE_lagged_long_pred <- predict(GEE_lagged_long,
                                         newdata = data_val)
data_val$GEE_lagged_long_bl_pred <- predict_baseline(GEE_lagged_long$coefficients,
                                                     newdata = data_val)

# GEE based estimation with informative covariance matrix
GEE_ar1_long <- geeglm(Marker ~ Time_Dummy + X1 + X2 + X3,
                       id = ID,
                       data = data_train,
                       family = gaussian,
                       corstr = "ar1",
)
data_val$GEE_long_ar1_pred <- predict(GEE_ar1_long,
                                      newdata = data_val)


# visualize individual trajectories
g_long_immortal_cohort <- ggplot(data = data_val %>%
                                   filter(ID %in% sampled_ID),
                                 aes(x = Time,
                                     y = Marker,
                                     group = ID
                                     # color = ID %>% as.factor()
                                     )) +
  geom_line(aes(color = "True trajectory")) +
  geom_line(aes(y = LMM_long_pred,
                color = "LMM landmarking prediction")) +
  geom_line(aes(y = LMM_long_dummy_pred,
                color = "LMM prediction (trend only)")) +
  geom_line(aes(y = LMM_long_bl_pred,
                color = "LMM baseline prediction")) +
  geom_line(aes(y = GEE_lagged_long_pred,
                color = "GEE landmarking prediction (IWC)")) +
  geom_line(aes(y = GEE_lagged_long_bl_pred,
                color = "GEE baseline prediction (IWC)")) +
  geom_line(aes(y = GEE_long_ar1_pred,
                color = "GEE prediction (AR1)")) +
  geom_point(data = data_val %>%
               filter(ID %in% sampled_ID) %>%
               filter(Death == 1) %>%
               group_by(ID) %>%
               filter(Time == min(Time)),
             aes(x = Time,
                 y = Marker),
             size = 2) +
  labs(color = "", linetype = "",
       title = "Immortal cohort") +
    facet_grid( ~ ID) +
  geom_line(data = data_val_mean_ID,
            aes(x = Time,
                y = mean_all,
                linetype = "Immortal mean"),
            color = "#333333") +
  geom_line(data = data_val_mean_ID,
            aes(x = Time,
                y = mean_surv,
                linetype = "Mortal mean"),
            color = "#333333")

g_long_immortal_cohort

# 1.2 Extrapolation beyond Death ------------------------------------------
## ON SURVIVAL DATA ONLY
# piecewise time trend only
LMM_long_surv_dummy <- lme(Marker ~ Time_Dummy + X1 + X2 + X3,
                           random = ~ 1 | ID,
                           data = data_train_surv)

# autoregressive additive component
LMM_long_surv <- lme(Marker ~ lagged_Marker + Time_Dummy + X1 + X2 + X3,
                     random = ~ 1 | ID,
                     data = data_train_surv)


# predictions according to LMM models; extrapolate beyond death
data_val$LMM_long_surv_dummy_pred <- predict(LMM_long_surv_dummy,
                                             newdata = data_val)
data_val$LMM_long_surv_pred <- predict(LMM_long_surv,
                                       newdata = data_val)
data_val$LMM_long_surv_bl_pred <- predict_baseline(LMM_long_surv$coefficients$fixed,
                                                   newdata = data_val)



# GEE based estimation with informative covariance matrix
GEE_ar1_surv_long <- geeglm(Marker ~ Time_Dummy + X1 + X2 + X3,
                            id = ID,
                            data = data_train_surv,
                            family = gaussian,
                            corstr = "ar1",
)
data_val$GEE_long_ar1_surv_pred <- predict(GEE_ar1_surv_long,
                                           newdata = data_val)


# visualize individual trajectories
g_long_extra_death <- ggplot(data = data_val %>%
                       filter(ID %in% sampled_ID),
                     aes(x = Time,
                         y = Marker,
                         group = ID)) +
  geom_line(aes(color = "True trajectory")) +
  geom_line(aes(y = LMM_long_surv_pred,
                color = "LMM landmarking prediction")) +
  geom_line(aes(y = LMM_long_surv_dummy_pred,
                color = "LMM prediction (trend only)")) +
  geom_line(aes(y = LMM_long_surv_bl_pred,
                color = "LMM baseline prediction")) +
  geom_line(aes(y = GEE_long_ar1_surv_pred,
                color = "GEE prediction (AR1)")) +
  geom_point(data = data_val %>%
               filter(ID %in% sampled_ID) %>%
               filter(Death == 1) %>%
               group_by(ID) %>%
               filter(Time == min(Time)),
             aes(x = Time,
                 y = Marker),
             size = 2) +
  labs(color = "", linetype = "",
       title = "Extrapolation beyond death") +
  facet_grid(~ ID) +
  geom_line(data = data_val_mean_ID,
            aes(x = Time,
                y = mean_all,
                linetype = "Immortal mean"),
            color = "#333333") +
  geom_line(data = data_val_mean_ID,
            aes(x = Time,
                y = mean_surv,
                linetype = "Mortal mean"),
            color = "#333333")

g_long_extra_death

# 1.3 Mortal cohort -------------------------------------------------------
# ON SURVIVOR DATA

# regression conditional on survival
RCA_long <- lm(Marker ~ Time_Dummy + X1 + X2 + X3,
               data = data_train_surv)
RCA_lagged_long <- lm(Marker ~ lagged_Marker + Time_Dummy + X1 + X2 + X3,
                      data = data_train_surv)

data_val$RCA_long_pred <- predict(RCA_long,
                                  newdata = data_val)
data_val$RCA_lagged_long_pred <- predict(RCA_lagged_long,
                                         newdata = data_val)
data_val$RCA_lagged_long_bl_pred <- predict_baseline(RCA_lagged_long$coefficients,
                                                     newdata = data_val)


# landmarking style; independent error between
GEE_surv_long <- geeglm(Marker ~ Time_Dummy + X1 + X2 + X3,
                               id = ID,
                               data = data_train_surv,
                               family = gaussian,
                               corstr = "independence")

GEE_lagged_surv_long <- geeglm(Marker ~ lagged_Marker + Time_Dummy + X1 + X2 + X3,
                               id = ID,
                               data = data_train_surv,
                               family = gaussian,
                               corstr = "independence")

data_val$GEE_long_surv_pred <-  predict(GEE_surv_long,
                                        newdata = data_val)
data_val$GEE_lagged_surv_long_pred <- predict(GEE_lagged_surv_long,
                                              newdata = data_val)
data_val$GEE_lagged_surv_long_bl_pred <- predict_baseline(GEE_lagged_surv_long$coefficients,
                                                          newdata = data_val)


# visualize individual trajectories
# visualize individual trajectories
g_long_mortal_cohort <- ggplot(data = data_val %>%
                                 filter(ID %in% sampled_ID),
                               aes(x = Time,
                                   y = Marker,
                                   group = ID)) +
  geom_line(aes(color = "True trajectory")) +
  geom_line(aes(y = RCA_long_pred,
                color = "RCA prediction (trend only)")) +
  geom_line(aes(y = RCA_lagged_long_pred,
                color = "RCA landmarking prediction")) +
  geom_line(aes(y = RCA_lagged_long_bl_pred,
                color = "RCA baseline prediction")) +
  geom_line(aes(y = GEE_long_surv_pred,
                color = "GEE prediction (trend only, IWC)")) +
  geom_line(aes(y = GEE_lagged_surv_long_pred,
                color = "GEE landmarking prediction (IWC)")) +
  geom_line(aes(y = LMM_long_surv_bl_pred,
                color = "GEE baseline prediction (IWC)")) +
  geom_point(data = data_val %>%
               filter(ID %in% sampled_ID) %>%
               filter(Death == 1) %>%
               group_by(ID) %>%
               filter(Time == min(Time)),
             aes(x = Time,
                 y = Marker),
             size = 2) +
  labs(color = "", linetype = "",
       title = "Partly conditional trajectory") +
  facet_grid(~ ID) +
  geom_line(data = data_val_mean_ID,
            aes(x = Time,
                y = mean_all,
                linetype = "Immortal mean"),
            color = "#333333") +
  geom_line(data = data_val_mean_ID,
            aes(x = Time,
                y = mean_surv,
                linetype = "Mortal mean"),
            color = "#333333")
g_long_mortal_cohort

# 1.4 Compare -------------------------------------------------------------
g_long_immortal_cohort +
  ylim(c(3, 30))

g_long_extra_death +
  ylim(c(3, 30))

g_long_mortal_cohort +
  ylim(c(3, 30))


ggarrange(g_long_immortal_cohort,
          g_long_extra_death,
          g_long_mortal_cohort,
          nrow = 3)



# 2. Mortality ------------------------------------------------------------

# Logistic regression -----------------------------------------------------
# # with baseline marker value
#
#
# ## SURVIVOR DATA ONLY
# # landmarking style with random effect
# LM_death_surv <- gLMMr(Death ~ Marker + Time_Dummy + X1 + X2 + X3 + (1 | ID),
#                        data = data_train_surv,
#                        family = binomial)
# LM_death_qrd_surv <- gLMMr(Death ~ Marker + I(Marker^2) + Time_Dummy + X1 + X2 + X3 + (1 | ID),
#                            data = data_train_surv,
#                            family = binomial)
# LM_death_lag_surv <- gLMMr(Death ~ Marker + lagged_Marker + Time_Dummy + X1 + X2 + X3 + (1 | ID),
#                            data = data_train_surv,
#                            family = binomial)
#
#
# data_val$LM_death_surv <- predict(LM_death_surv,
#                              newdata = data_val,
#                              re.form = NA,
#                              type = "response")
# data_val$LM_death_qrd_surv <- predict(LM_death_qrd_surv,
#                                  newdata = data_val,
#                                  re.form = NA,
#                                  type = "response")
# data_val$LM_death_lag_surv <- predict(LM_death_lag_surv,
#                                  newdata = data_val,
#                                  re.form = NA,
#                                  type = "response")
#
# ggplot(data = data_val %>%
#          filter(ID %in% sampled_ID),
#        aes(x = Time,
#            y = Death,
#            group = ID,
#            color = ID %>% as.factor())) +
#   geom_line(aes(linetype = "Death indicator")) +
#   geom_line(aes(y = LM_death_surv,
#                 linetype = "Landmarking model")) +
#   geom_line(aes(y = LM_death_qrd_surv,
#                 linetype = "Landmarking model (qrd. marker)")) +
#   geom_line(aes(y = LM_death_lag_surv,
#                 linetype = "Landmarking model (lagged marker)")) +
#   labs(color = "", linetype = "")
#
#
#
#

# Mortality estimating equations ------------------------------------------




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
data_id_dead <- data_train %>%
  group_by(ID) %>%
  mutate(status = ifelse(any(Death == 1), 1, 0)) %>%
  filter(status == 1) %>%
  filter(Time == min(Time))

data_id_alive <- data_train %>%
  group_by(ID) %>%
  mutate(status = ifelse(any(Death == 1), 1, 0)) %>%
  filter(status == 0) %>%
  filter(Time == max(Time))

data_train_id <- rbind(data_id_dead,
                       data_id_alive)


JM_pw_mod <- jointModel(lmeObject = lme(Marker ~ Mort_Time + X1 + X2 + X3,
                                        random = ~ 1 | ID,
                                        data = data_train_surv),
           survObject = coxph(Surv(Mort_Time, status) ~ X1 + X2 + X3,
                              data = data_train_id,
                              x = TRUE),
           timeVar = "Mort_Time",
           parameterization = "value",
           method = "piecewise-PH-GH",
           control = list(knots = 15,
                          verbose = TRUE))


JM_mod <- jointModel(lmeObject = lme(Marker ~ Mort_Time + X1 + X2 + X3,
                           random = ~ 1 | ID,
                           data = data_train_surv),
           survObject = coxph(Surv(Mort_Time, status) ~ 1 + X1 + X2 + X3,
                              data = data_train_id,
                              x = TRUE),
           timeVar = "Mort_Time",
           parameterization = "value",
           method = "Cox-PH-GH",

           control = list(verbose = TRUE)
           )
summary(JM_mod)


# Marginal model ----------------------------------------------------------
# according to Kurland 2005

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


