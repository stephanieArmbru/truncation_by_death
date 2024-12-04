############################### DATA SIMULATION ################################

## Stephanie Armbruster
## 10/20/2024

## This script simulates a data set on which the methods are applied.

# Library -----------------------------------------------------------------
set.seed(123)
library(tidyverse)
library(rlist)

# setwd("Log_Reg_Simulation")
source("Functions.R")

# Simulation --------------------------------------------------------------
# true OR coefficients
beta0 <- 1
beta1 <- 2.5
beta2 <- 4
beta12 <- 3

# sample size
N_train <- 100000
N_val <- 10000

train_dat <- simDat(N = N_train,
                    beta0 = beta0, beta1 = beta1, beta2 = beta2,
                    beta12 = beta12)

val_dat <- simDat(N = N_val,
                  beta0 = beta0, beta1 = beta1, beta2 = beta2,
                  beta12 = beta12)

# transform into data frames
train_dat_df <- data.frame(Y = train_dat$Y) %>%
  cbind(train_dat$X)

val_dat_df <- data.frame(Y = val_dat$Y) %>%
  cbind(val_dat$X)


# distribution of probabilities
ggplot(data = data.frame(x = expit(train_dat$lin_pred)),
       aes(x = x)) +
  geom_density() +
  theme_bw() +
  labs(x = "Probability")

ggplot(data = data.frame(x = expit(train_dat$lin_pred),
                         y = train_dat$Y),
       aes(x = x,
           y = y)) +
  geom_point() +
  theme_bw() +
  labs(x = "Probability",
       y = "Outcome")


# Grid search -------------------------------------------------------------
# grid over which the search is performed
grid_search <- expand.grid(beta0 = seq(-3, 3, by = 0.5),
                           beta1 = seq(-2, 10, by = 0.5),
                           beta2 = seq(-2, 10, by = 0.5),
                           beta12 = seq(-2, 10, by = 0.5)) %>%
  as.matrix()

# split grid search into subsearchs
split_size <- 150
row_indices <- rep(1:split_size,
                   each = ceiling(nrow(grid_search) / split_size))

grid_search_subsets_list <- split(grid_search,
                                  f = row_indices[seq(1, nrow(grid_search))]) %>%
  lapply(function(l) matrix(l, ncol = ncol(grid_search)))

# check if grid split works and retains all grid combinations
identical(grid_search[, 1], list.rbind(grid_search_subsets_list)[, 1])
identical(grid_search[, 2], list.rbind(grid_search_subsets_list)[, 2])
identical(grid_search[, 3], list.rbind(grid_search_subsets_list)[, 3])
identical(grid_search[, 4], list.rbind(grid_search_subsets_list)[, 4])


# Save --------------------------------------------------------------------
save(train_dat_df, file = "Data/Sim_Dat_Training.RData")
save(val_dat_df, file = "Data/Sim_Dat_Validation.RData")

beta_vector <- c(beta0, beta1, beta2, beta12)
save(beta_vector, file = "Data/True_beta.RData")

save(grid_search_subsets_list, file = "Data/Grid_search_subsets.RData")

