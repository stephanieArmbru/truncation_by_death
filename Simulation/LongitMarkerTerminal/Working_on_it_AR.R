### AR MODELS ###

rm(list = ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(Rcpp)
library(ggpubr)
set.seed(1234)

# Load functions ----------------------------------------------------------
source("Simulation/ShinySimulation/R/Simulation_analysis_helpers.R")

source("Simulation/LongitMarkerTerminal/R/LongitMTparam.R")

source("Simulation/LongitMarkerTerminal/R/ParamLongitMTLogLik.R")
source("Simulation/LongitMarkerTerminal/R/GradParamLongitMTLogLik.R")

source("Simulation/LongitMarkerTerminal/R/Helper.R")
source("Simulation/LongitMarkerTerminal/R/SimDatARLongitTerminal.R")



# Define your custom theme
custom_theme <- theme(
  legend.position = "bottom"
)

# Set the custom theme as the default for all plots
theme_set(theme_bw() + custom_theme)



# Simulation parameters ---------------------------------------------------
## Global parameters
n <- 100000  # Sample size
times <- seq(0, 10, by = 0.25)  # Discrete Time Partition
N_times <- length(times)

# Mortality Parameters
slope_t_threshold <- -10 # Threshold for slope for add. marker effect
long_t_threshold <- 21 # Threshold for marker value for add. marker effect (scaled by distance to threshold)

lambda_t <- rep(logit(0.02), N_times - 1)  # Logit baseline mortality probability
theta_t <- 0.02  # Global linear marker effect on mortality probability
varphi_t <- 0.05 + 0.001 * times[-1]  # Local linear marker effect on mortality probability
varphi_ReLU_t <- rep(1.5, N_times - 1) # Local marker effect acc. to threshold
varphi_slope_t <- rep(2, N_times - 1) # Local marker slope effect acc. to threshold
xi_t <- c(1, -2, 3)  # Time-independent covariate effect on mortality probability

# Longitudinal Marker Parameters
zeta_long <- c(25, rep(3, N_times - 1))
eta_long <- rep(0.85, N_times - 1) # Local autoregressive component
beta_long <- c(0.5, 0.2, -0.1)  # Time-independent covariate effect on marker
sd_long_trajectory <- 0.25  # Standard deviation for longitudinal marker


# Simulation --------------------------------------------------------------
data <- SimDatARLongitTerminal(n.sample = n,
                               times = times,
                               p = 3,
                               lambda.t = lambda_t,
                               slope.t.threshold = slope_t_threshold,
                               long.t.threshold = long_t_threshold,
                               theta.t = theta_t,
                               varphi.t = varphi_t,
                               varphi.slope.t = varphi_slope_t,
                               varphi.ReLU.t = varphi_ReLU_t,
                               xi.t = xi_t,
                               zeta.long = zeta_long,
                               eta.long = eta_long,
                               beta.long = beta_long,
                               sd.long.trajectory = sd_long_trajectory)

plot_patients_at_risk(data$df.return.M)

idx <- sample(seq(1, n), size = 10)
g1 <- plot_heterogeneity(data = data$df.return.M,
                   sample_IDs = idx)
g2 <- plot_cum_mortality_prob(data$df.return.M,
                        sample_IDs = idx)
ggarrange(g1, g2,
          common.legend = TRUE)

plot_mortality_prob(data = data$df.return.M,
                    sample_IDs = idx)

plot_immortal_mortal_longitudinal_trajectory(data = data$df.return.IM)













