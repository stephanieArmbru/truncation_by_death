### WORKING ON THE ESTIION OF THE NEW METHOD ###

rm(list = ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(Rcpp)
set.seed(1234)

# Load functions ----------------------------------------------------------
# devtools::load_all("Simulation/NevoPackage/R")

source("Simulation/ShinySimulation/R/Data_generating_mechanism.R")
source("Simulation/ShinySimulation/R/Simulation_analysis_helpers.R")

source("Simulation/LongitMarkerTerminal/R/LongitMTparam.R")

source("Simulation/LongitMarkerTerminal/R/ParamLongitMTLogLik.R")
source("Simulation/LongitMarkerTerminal/R/GradParamLongitMTLogLik.R")

source("Simulation/LongitMarkerTerminal/R/Helper.R")


# RCPP
# Load C++ function
# sourceCpp("Simulation/LongitMarkerTerminal/Source/ParamLongitMTLogLik.cpp")
# sourceCpp("Simulation/LongitMarkerTerminal/Source/GradParamLongitMTLogLik.cpp")


# Define your custom theme
custom_theme <- theme(
  legend.position = "bottom"
)

# Set the custom theme as the default for all plots
theme_set(theme_bw() + custom_theme)



# Simulation parameters ---------------------------------------------------
## Global parameters
n <- 1000  # Sample size
times <- seq(0, 10, by = 0.25)  # Discrete Time Partition
N_times <- length(times)

# Mortality Parameters
slope_t_threshold <- -2 # Threshold for slope for add. marker effect
long_t_threshold <- 21 # Threshold for marker value for add. marker effect (scaled by distance to threshold)


lambda_t <- rep(logit(0.02), N_times - 1)  # Logit baseline mortality probability
theta_t <- 0.02  # Global linear marker effect on mortality probability
varphi_t <- 0.05 + 0.001 * times[-1]  # Local linear marker effect on mortality probability
varphi_ReLU_t <- rep(1.5, N_times - 1) # Local marker effect acc. to threshold
varphi_slope_t <- rep(2, N_times - 1) # Local marker slope effect acc. to threshold
xi_t <- c(1, -2, 3)  # Time-independent covariate effect on mortality probability

# Longitudinal Marker Parameters
long_threshold <- 24  # Threshold for marker change
zeta_long <- 25 - 0.015 * times^2
# c(-(0.15 * times[1:17])^2 + 25,
#              rep(min(-(0.2 * times[1:17])^2 + 25), N_times - 17)) # Baseline marker trend
zeta_ReLU_long <- c(rep(-1.5, N_times))  # Threshold trend effect

# stationary
eta_long <- runif(n = N_times - 1, min = 0.9, max = 0.9) # Local autoregressive component
slope_threshold <- -1
eta_slope_long <- rep(-1.5, N_times - 1)  # Local autoregressive slope effect


beta_long <- c(0.5, 0.2, -0.1)  # Time-independent covariate effect on marker
sd_long_trajectory <- 0.25  # Standard deviation for longitudinal marker


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

res_df <- res$df.return.M %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(Marker_lagged = c(0, na.omit(dplyr::lag(Marker))),
         Marker_lagged_2 = c(0, 0, na.omit(dplyr::lag(Marker,
                                               n = 2)))
  ) %>%
  ungroup() %>%
  arrange(Time, ID)

# add historical transforions for longitudinal and terminal event model
data.list <- list(YT = res_df$Death,
                  Ylong = res_df$Marker,
                  risk.T = 1 - res_df$Death,

                  # number of discrete time interval (not true time)
                  TM = res_df$TM,
                  Time = res_df$Time,
                  Time.sq = res_df$Time^2,
                  ID = res_df$ID,


                  Ylong.lagged = res_df$Marker_lagged,
                  Ylong.lagged.2 = res_df$Marker_lagged_2,
                  TimeYlong.lagged = res_df$Marker_lagged_2 * res_df$Time,
                  X1 = res_df$X1,
                  X2 = res_df$X2,
                  X3 = res_df$X3)


plot_longitudinal_trend(times, zeta_long)
plot_heterogeneity(data = res$df.return.IM)
plot_patients_at_risk(data = res$df.return.IM)
plot_mortality_prob(data = res$df.return.IM)


# Optimization ------------------------------------------------------------
# longitudinal trend
LMM.long <- lm(Marker ~ Time + I(Time^2) + Marker_lagged + Marker_lagged_2 + I(Time * Marker_lagged) + X1 + X2 + X3,
   data = res_df)

LR.T <- glm(Death ~ Time +  I(Time^2)  + Marker_lagged + Marker_lagged_2 + I(Time * Marker_lagged) + X1 + X2 + X3,
            data = res_df,
            family = binomial)

zeta <- c(LMM.long$coefficients["(Intercept)"],
          LMM.long$coefficients["Time"],
          LMM.long$coefficients["I(Time^2)"])
eta <- c(LMM.long$coefficients["Marker_lagged"],
         LMM.long$coefficients["Marker_lagged_2"],
         LMM.long$coefficients["I(Time * Marker_lagged)"])
beta <- c(LMM.long$coefficients["X1"],
          LMM.long$coefficients["X2"],
          LMM.long$coefficients["X3"])

sigma.sq <- var(res_df$Marker - predict(LMM.long, res_df))

lambda <- LR.T$coefficients[("(Intercept)")]
theta <- c(LR.T$coefficients[("Marker_lagged")],
           LR.T$coefficients[("Marker_lagged_2")])
varphi <- c(LR.T$coefficients[("Time")],
            LR.T$coefficients[("I(Time^2)")],
            LR.T$coefficients[("I(Time * Marker_lagged)")])
xi <- c(LR.T$coefficients[("X1")],
        LR.T$coefficients[("X2")],
        LR.T$coefficients[("X3")])


# parameters for optim()
init <- c(zeta, eta, beta, sigma.sq,
          lambda, theta, varphi, xi)

longit.data <- data.list

formula.long <- "~ X1 + X2 + X3" %>% as.formula()
formula.long.AR <- "~ Ylong.lagged + Ylong.lagged.2 + TimeYlong.lagged" %>% as.formula()
formula.long.trend <- "~ Time + Time.sq" %>% as.formula()

formula.T <- "~ X1 + X2 + X3" %>% as.formula()
formula.glob.dep <- "~ Ylong.lagged + Ylong.lagged.2" %>% as.formula()
formula.loc.dep <- "~ Time + Time.sq + TimeYlong.lagged" %>% as.formula()


optimResults <- LongitMTparam(longit.data,
                              times,
                              formula.long,
                              formula.long.AR,
                              formula.long.trend,

                              formula.T,
                              formula.glob.dep,
                              formula.loc.dep,

                              init = init,
                              maxit.optim = 50000)



# Prediction --------------------------------------------------------------
# simulate validation data set
res_validation <- SimulateJointLongData(
  n.sample = 10000,
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

res_val_df <- res_validation$df.return.M %>%
  group_by(ID) %>%
  arrange(Time) %>%
  mutate(Marker_lagged = c(0, na.omit(dplyr::lag(Marker))),
         Marker_lagged_2 = c(0, 0, na.omit(dplyr::lag(Marker,
                                                      n = 2)))
  ) %>%
  ungroup() %>%
  arrange(Time, ID)

# add historical transforions for longitudinal and terminal event model
data.list.val <- list(YT = res_val_df$Death,
                  Ylong = res_val_df$Marker,

                  # number of discrete time interval (not true time)
                  TM = res_val_df$TM,
                  Time = res_val_df$Time,
                  Time.sq = res_val_df$Time^2,
                  ID = res_val_df$ID,


                  Ylong.lagged = res_val_df$Marker_lagged,
                  Ylong.lagged.2 = res_val_df$Marker_lagged_2,
                  TimeYlong.lagged = res_val_df$Marker_lagged_2 * res_val_df$Time,
                  X1 = res_val_df$X1,
                  X2 = res_val_df$X2,
                  X3 = res_val_df$X3)






# Landmark prediction
LandmarkPred <- predictLongitMTparam.Landmark(par = optimResults$est,
                                     longit.data = data.list.val,
                                     formula.long,
                                     formula.long.AR,
                                     formula.long.trend,

                                     formula.T,
                                     formula.glob.dep,
                                     formula.loc.dep)

res_val_df$Pred_Mort_Prob <- LandmarkPred$PredProb
res_val_df$Pred_Marker <- LandmarkPred$MarkerPred


ggplot(data = res_val_df %>%
         filter(ID %in% seq(1, 10),
                TM >= 3),
       aes(x = Time,
           color = ID %>% as.factor())) +
  geom_line(aes(y = Mort_Prob,
                linetype = "True event probability")) +
  geom_line(aes(y = Pred_Mort_Prob,
                linetype = "Predicted event probability")) +
  labs(color = "",
       linetype = "",
       y = "Mortality Probability")

ggplot(data = res_val_df %>%
         filter(ID %in% seq(1, 10),
                TM >= 3),
       aes(x = Time,
           color = ID %>% as.factor())) +
  geom_line(aes(y = Marker,
                linetype = "True traj.")) +
  geom_line(aes(y = Pred_Marker,
                linetype = "Pred. traj.")) +
  labs(color = "",
       linetype = "")


# Baseline prediction
Xlong <- res_validation$X %>%
  select(-ID)
XT <- res_validation$X %>%
  select(-ID)


BaselinePred <- predictLongitMTparam.Baseline(par = optimResults$est,
                                              Xlong = Xlong,
                                              XT = XT,
                                              Ylong = res_val_df$Marker,
                                              YT = res_val_df$Death,

                                              p.long.AR = length(eta),
                                              p.long.Trend = length(zeta),
                                              p.T.globdep = length(theta),
                                              p.T.locdep = length(varphi),

                                              # time and ID information
                                              ID = res_val_df$ID,
                                              TM = res_val_df$TM,
                                              time = unique(res_val_df$Time))

ggplot(data = BaselinePred %>%
         filter(ID %in% seq(1, 10),
                Time >= 0.5),
       aes(x = Time,
           y = Marker,
           group = ID,
           color = ID %>% as.factor(),
           linetype ="Pred. traj.")) +
  geom_line() +
  geom_line(data = res_val_df %>%
              filter(ID %in% seq(1, 10),
                     TM >= 3),
            aes(x = Time,
                y = Marker,
                color = ID %>% as.factor(),
                linetype = "True traj.")) +
  labs(color = "",
       linetype = "",
       y = "Marker [trunc. y-axis]") +
  ylim(-40, 40)

ggplot(data = BaselinePred %>%
         filter(ID %in% seq(1, 10),
                Time >= 0.5),
       aes(x = Time,
           y = EventProb,
           group = ID,
           color = ID %>% as.factor(),
           linetype ="Pred. event prob.")) +
  geom_line() +
  geom_line(data = res_val_df %>%
              filter(ID %in% seq(1, 10),
                     TM >= 3),
            aes(x = Time,
                y = Mort_Prob,
                color = ID %>% as.factor(),
                linetype = "True event prob.")) +
  labs(color = "",
       linetype = "",
       y = "Mortality Probability")

# Manual check ------------------------------------------------------------
# generate parameters for function
init <- NULL
maxit.optim <- 5000


