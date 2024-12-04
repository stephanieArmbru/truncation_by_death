############################## SET-UP FOR CLUSTER ##############################

## Stephanie Armbruster
## 10/20/2024

## This is the script which will be run on the cluster based on inputs.

rm(list=ls())
# setwd("Log_Reg_Simulation")

# Parameter input ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
g <- as.integer(as.numeric(args[1], na.rm = TRUE))

g <- 1

# Library -----------------------------------------------------------------
library(tidyverse)
library(pROC)

library(foreach)      # For parallel looping
library(doParallel)

library(ggpubr)
library(rlist)

# load data set
load("Data/Sim_Dat_Training.RData")
load("Data/Grid_search_subsets.RData")

load("Data/Grid_search_subsets.RData")

# source functions
source("Functions.R")

# parallel
# cores <- 10
cores <- 10


# Preparation -------------------------------------------------------------
# design matrix for training data set
Xdesign <- data.frame(Intercept = 1,
                      X1 = train_dat_df$X1,
                      X2 = train_dat_df$X2,
                      Interaction = train_dat_df$X1 * train_dat_df$X2) %>%
  as.matrix()

# for model misspecification
Xdesign_IC <- Xdesign[, -4]


# use grid subset
grid_search_subset <- grid_search_subsets_list[[g]]
n_steps <- nrow(grid_search_subset)

grid_search_subset_IC <- grid_search_subset[, -4] %>% unique()
n_steps_IC <- nrow(grid_search_subset_IC)

# AUC ---------------------------------------------------------------------
# correct model specification
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform grid search with parallelization
auc_results <- foreach(step = 1:n_steps,
                       .combine = rbind,
                       .packages = c("pROC", "dplyr")) %dopar% {
                         auc_value <- AUC_optimization(data = train_dat_df,
                                                       step = step,
                                                       grid_search = grid_search_subset,
                                                       Xdesign = Xdesign)
                         c(step, grid_search_subset[step, ], auc_value)  # Return step and corresponding AUC value
                       }

# Stop the parallel backend
stopCluster(cl)

# incorrect model specification
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform grid search with parallelization
auc_results_IC <- foreach(step = 1:n_steps_IC,
                          .combine = rbind,
                          .packages = c("pROC", "dplyr")) %dopar% {
                            auc_value <- AUC_optimization(data = train_dat_df,
                                                          step = step,
                                                          grid_search = grid_search_subset_IC,
                                                          Xdesign = Xdesign_IC)
                            c(step, grid_search_subset_IC[step, ], auc_value)  # Return step and corresponding AUC value
                          }

# Stop the parallel backend
stopCluster(cl)


# ECE ---------------------------------------------------------------------

# correct model specification
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform grid search with parallelization
ece_results <- foreach(step = 1:n_steps,
                       .combine = rbind,
                       .packages = c("dplyr")) %dopar% {
                         ece_value <- ECE_optimization(data = train_dat_df,
                                                       step = step,
                                                       grid_search = grid_search_subset,
                                                       Xdesign = Xdesign)
                         c(step, grid_search_subset[step, ], ece_value)  # Return step and corresponding AUC value
                       }

# Stop the parallel backend
stopCluster(cl)

# incorrect model specification
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform grid search with parallelization
ece_results_IC <- foreach(step = 1:n_steps_IC,
                          .combine = rbind,
                          .packages = c("dplyr")) %dopar% {
                            ece_value <- ECE_optimization(data = train_dat_df,
                                                          step = step,
                                                          grid_search = grid_search_subset_IC,
                                                          Xdesign = Xdesign_IC)
                            c(step, grid_search_subset_IC[step, ], ece_value)  # Return step and corresponding AUC value
                          }

# Stop the parallel backend
stopCluster(cl)



# Brier score -------------------------------------------------------------
# correct model specification
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform grid search with parallelization
brier_results <- foreach(step = 1:n_steps,
                         .combine = rbind,
                         .packages = c("dplyr")) %dopar% {
                           brier_value <- Brier_optimization(data = train_dat_df,
                                                             step = step,
                                                             grid_search = grid_search_subset,
                                                             Xdesign = Xdesign)
                           c(step, grid_search_subset[step, ], brier_value)  # Return step and corresponding AUC value
                         }

# Stop the parallel backend
stopCluster(cl)

# incorrect model specification
cl <- makeCluster(cores)
registerDoParallel(cl)

# Perform grid search with parallelization
brier_results_IC <- foreach(step = 1:n_steps_IC,
                            .combine = rbind,
                            .packages = c("dplyr")) %dopar% {
                              brier_value <- Brier_optimization(data = train_dat_df,
                                                                step = step,
                                                                grid_search = grid_search_subset_IC,
                                                                Xdesign = Xdesign_IC)
                              c(step, grid_search_subset_IC[step, ], brier_value)  # Return step and corresponding AUC value
                            }

# Stop the parallel backend
stopCluster(cl)



# Post-processing ---------------------------------------------------------
# find optimal among grid search subset
AUC_results <- auc_results[auc_results[, ncol(auc_results)] == max(auc_results[, ncol(auc_results)]), ]
AUC_results_IC <- auc_results_IC[auc_results_IC[, ncol(auc_results_IC)] == max(auc_results_IC[, ncol(auc_results_IC)]), ]

ECE_results <- ece_results[ece_results[, ncol(ece_results)] == max(ece_results[, ncol(ece_results)]), ]
ECE_results_IC <- ece_results_IC[ece_results_IC[, ncol(ece_results_IC)] == min(ece_results_IC[, ncol(ece_results_IC)]), ]

Brier_results <- brier_results[brier_results[, ncol(brier_results)] == min(brier_results[, ncol(brier_results)]), ]
Brier_results_IC <- brier_results_IC[brier_results_IC[, ncol(brier_results_IC)] == min(brier_results_IC[, ncol(brier_results_IC)]), ]

# save additional information for variability of results: min/max and median
AUC_results_info <- rbind(auc_results[auc_results[, ncol(auc_results)] == min(auc_results[, ncol(auc_results)]), ],
                          auc_results[auc_results[, ncol(auc_results)] == median(auc_results[, ncol(auc_results)]), ])
AUC_results_IC_info <- rbind(auc_results_IC[auc_results_IC[, ncol(auc_results_IC)] == min(auc_results_IC[, ncol(auc_results_IC)]), ],
                             auc_results_IC[auc_results_IC[, ncol(auc_results_IC)] == median(auc_results_IC[, ncol(auc_results_IC)]), ])

ECE_results_info <- rbind(ece_results[ece_results[, ncol(ece_results)] == max(ece_results[, ncol(ece_results)]), ],
                          ece_results[ece_results[, ncol(ece_results)] == median(ece_results[, ncol(ece_results)]), ])
ECE_results_IC_info <- rbind(ece_results_IC[ece_results_IC[, ncol(ece_results_IC)] == max(ece_results_IC[, ncol(ece_results_IC)]), ],
                             ece_results_IC[ece_results_IC[, ncol(ece_results_IC)] == median(ece_results_IC[, ncol(ece_results_IC)]), ])

Brier_results_info <- rbind(brier_results[brier_results[, ncol(brier_results)] == max(brier_results[, ncol(brier_results)]), ],
                            brier_results[brier_results[, ncol(brier_results)] == median(brier_results[, ncol(brier_results)]), ])
Brier_results_IC_info <- rbind(brier_results_IC[brier_results_IC[, ncol(brier_results_IC)] == max(brier_results_IC[, ncol(brier_results_IC)]), ],
                               brier_results_IC[brier_results_IC[, ncol(brier_results_IC)] == median(brier_results_IC[, ncol(brier_results_IC)]), ])


# Save --------------------------------------------------------------------
save(AUC_results, file = paste("Cluster_results/AUC_results_subset_", g, ".RData"))
save(AUC_results_IC, file = paste("Cluster_results/AUC_IC_results_subset_", g, ".RData"))

save(ECE_results, file = paste("Cluster_results/ECE_results_subset_", g, ".RData"))
save(ECE_results_IC, file = paste("Cluster_results/ECE_IC_results_subset_", g, ".RData"))

save(Brier_results, file = paste("Cluster_results/Brier_results_subset_", g, ".RData"))
save(Brier_results_IC, file = paste("Cluster_results/Brier_IC_results_subset_", g, ".RData"))

# save info overviews
save(AUC_results_info, file = "Cluster_results/Add_info/AUC_results_info_subset_", g, ".RData")
save(AUC_results_IC_info, file = "Cluster_results/Add_info/AUC_results_IC_info_subset_", g, ".RData")
save(ECE_results_info, file = "Cluster_results/Add_info/ECE_results_info_subset_", g, ".RData")
save(ECE_results_IC_info, file = "Cluster_results/Add_info/ECE_results_IC_info_subset_", g, ".RData")
save(Brier_results_info, file = "Cluster_results/Add_info/Brier_results_info_subset_", g, ".RData")
save(Brier_results_IC_info, file = "Cluster_results/Add_info/Brier_results_IC_info_subset_", g, ".RData")




