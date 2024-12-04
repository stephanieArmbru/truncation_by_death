############################ ANALYZE CLUSTER OUTPUT ############################

## Stephanie Armbruster
## 10/20/2024


## This script fits the ML estimation for correct and incorrect model specification
## and uses the cluster output to find the optimum across the criteria.

# Load data ---------------------------------------------------------------
load("Data/Sim_Dat_Training.RData")
load("Data/Sim_Dat_Validation.RData")


load("Data/True_beta.RData")
load("Data/Grid_search_subsets.RData")

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

read_cluster_result <- function(Method, split_size) {
  data <- list()
  for (g in seq(1, split_size)) {
    data[[g]] <- loadRData(paste("Cluster_results/", Methods, "_results_subset_",
                    g, ".RData", sep = ""))
  }


  }

# design matrix
Xdesign_val <- data.frame(Intercept = 1,
                          X1 = val_dat_df$X1,
                          X2 = val_dat_df$X2,
                          Interaction = val_dat_df$X1 * val_dat_df$X2) %>%
  as.matrix()

# design matrix for misspecified model
Xdesign_val_IC <- Xdesign_val[, -4]


# ML correct model specification ------------------------------------------

# model fitting according to IWLS
ML_mod <- glm(Y ~ 1 + X1 + X2 + X1*X2,
              data = train_dat_df,
              family = "binomial")

# model coefficients
ML_mod_coef <- summary(ML_mod)$coef %>%
  as.data.frame() %>%
  mutate(True = c(beta0, beta1, beta2 , beta12),
         Coef = rownames(summary(ML_mod)$coef))

list_coefficients[["ML_CMS"]] <- ggplot(data = ML_mod_coef,
                                        aes(x = Coef,
                                            y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`,
                    ymax = Estimate + 1.96 * `Std. Error`)) +
  geom_point(aes(y = True,
                 color = "True coefficient"),
             shape = 4, size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Model coefficient",
       color = "",
       caption = "ML estimation, correct") +
  coord_flip()

# values for training data
AUC_train_ML <- roc(train_dat_df$Y, ML_mod$fitted.values)$auc
ECE_train_ML <- mean(abs(train_dat_df$Y - ML_mod$fitted.value))
Brier_train_ML <- mean((train_dat_df$Y - ML_mod$fitted.value)^2)

# predict probabilities for validation data
val_dat_df$pred <- predict(ML_mod,
                           newdata = val_dat_df,
                           type = "response")


# calibration plot
list_calibration_plots[["ML_CMS"]] <- calibrationPlot(data = val_dat_df,
                                                      caption = "ML estimation, correct")


# ROC
list_ROC[["ML_CMS"]] <- ROCplot(data = val_dat_df,
                                caption = "ML estimation, correct")

# ECE / Brier score
list_ECE[["ML_CMS"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["ML_CMS"]] <- Brier(val_dat_df$Y, val_dat_df$pred)


# ML incorrect model specification ----------------------------------------
# model fitting according to IWLS
ML_mod_IC <- glm(Y ~ 1 + X1 + X2,
                 data = train_dat_df,
                 family = "binomial")

# model coefficients
ML_mod_coef_IC <- summary(ML_mod_IC)$coef %>%
  as.data.frame() %>%
  mutate(True = c(beta0, beta1, beta2),
         Coef = rownames(summary(ML_mod_IC)$coef))

list_coefficients[["ML_ICMS"]]  <- ggplot(data = ML_mod_coef_IC,
                                          aes(x = Coef,
                                              y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`,
                    ymax = Estimate + 1.96 * `Std. Error`)) +
  geom_point(aes(y = True,
                 color = "True coefficient"),
             shape = 4, size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Model coefficient",
       color = "",
       caption = "ML estimation, incorrect") +
  coord_flip()

# values for training data
AUC_train_ICML <- roc(train_dat_df$Y, ML_mod_IC$fitted.values)$auc
ECE_train_ICML <- mean(abs(train_dat_df$Y - ML_mod_IC$fitted.value))
Brier_train_ICML <- mean((train_dat_df$Y - ML_mod_IC$fitted.value)^2)

# predict probabilities for validation data
val_dat_df$pred <- predict(ML_mod_IC,
                           newdata = val_dat_df,
                           type = "response")



# calibration plot
list_calibration_plots[["ML_ICMS"]] <- calibrationPlot(data = val_dat_df,
                                                       caption = "ML estimation, incorrect")


# ROC
list_ROC[["ML_ICMS"]] <- ROCplot(data = val_dat_df,
                                 caption = "ML estimation, incorrect")

# ECE / Brier score
list_ECE[["ML_ICMS"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["ML_ICMS"]] <- Brier(val_dat_df$Y, val_dat_df$pred)


# AUC correct model specification -----------------------------------------
auc_results <- read_cluster_results(method = "AUC")

AUC_results_df <- data.frame(step = auc_results[, 1],
                             AUC = auc_results[, 2]) %>%
  mutate(beta0 = grid_search[step, 1],
         beta1 = grid_search[step, 2],
         beta2 = grid_search[step, 3],
         beta12 = grid_search[step, 4])

AUC_results_df_long <- AUC_results_df %>%
  gather(coef, value, -c(step, AUC))


# visualize results
ggplot(data = AUC_results_df_long,
       aes(x = value,
           y = AUC,
           color = coef)) +
  geom_point() +
  geom_hline(aes(yintercept = AUC_train_ML)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ coef)

# Find the best parameters (maximizing AUC)
best_AUC <- AUC_results_df %>% filter(AUC == max(AUC))


# best beta
AUC_beta <- c(0, best_AUC$beta1,
              best_AUC$beta2, best_AUC$beta12) %>%
  unique()

# predict for best beta
val_dat_df$pred <- predict_for_beta(Xdesign = Xdesign_val,
                                    beta = AUC_beta)


# calibration plot
list_calibration_plots[["AUC"]] <- calibrationPlot(data = val_dat_df,
                                                   caption = "AUC optimization")


# ROC
list_ROC[["AUC"]] <- ROCplot(data = val_dat_df,
                             caption = "AUC optimization")

# ECE / Brier score
list_ECE[["AUC"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["AUC"]] <- Brier(val_dat_df$Y, val_dat_df$pred)


# AUC incorrect model specification ---------------------------------------
AUC_results_IC_df <- data.frame(step = auc_results_IC[, 1],
                                AUC = auc_results_IC[, 2]) %>%
  mutate(beta0 = grid_search_IC[step, 1],
         beta1 = grid_search_IC[step, 2],
         beta2 = grid_search_IC[step, 3])

AUC_results_IC_df_long <- AUC_results_IC_df %>%
  gather(coef, value, -c(step, AUC))


# visualize results
ggplot(data = AUC_results_IC_df_long,
       aes(x = value,
           y = AUC,
           color = coef)) +
  geom_point() +
  geom_hline(aes(yintercept = AUC_train_ICML)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ coef)

# Find the best parameters (maximizing AUC)
best_AUC_IC <- AUC_results_IC_df %>% filter(AUC == max(AUC))


# best beta
AUC_beta_IC <- c(0, best_AUC_IC$beta1, best_AUC_IC$beta2) %>%
  unique()

# predict for best beta
val_dat_df$pred <- predict_for_beta(Xdesign = Xdesign_val_IC,
                                    beta = AUC_beta_IC)


# calibration plot
list_calibration_plots[["AUC_IC"]] <- calibrationPlot(data = val_dat_df,
                                                      caption = "AUC optimization, incorrect")


# ROC
list_ROC[["AUC_IC"]] <- ROCplot(data = val_dat_df,
                                caption = "AUC optimization, incorrect")

# ECE / Brier score
list_ECE[["AUC_IC"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["AUC_IC"]] <- Brier(val_dat_df$Y, val_dat_df$pred)


# ECE correct model specification -----------------------------------------
ece_results <- read_cluster_results(method = "ECE")

ECE_results_df <- data.frame(step = ece_results[, 1],
                             ECE = ece_results[, 2]) %>%
  mutate(beta0 = grid_search[step, 1],
         beta1 = grid_search[step, 2],
         beta2 = grid_search[step, 3],
         beta12 = grid_search[step, 4])

ECE_results_df_long <- ECE_results_df %>%
  gather(coef, value, -c(step, ECE))


# visualize results
ggplot(data = ECE_results_df_long,
       aes(x = value,
           y = ECE,
           color = coef)) +
  geom_point() +
  geom_hline(aes(yintercept = ECE_train_ML)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ coef)

# Find the best parameters (maximizing AUC)
best_ECE <- ECE_results_df %>% filter(ECE == min(ECE))


# best beta
ECE_beta <- c(best_ECE$beta0, best_ECE$beta1,
              best_ECE$beta2, best_ECE$beta12)

# predict for best beta
val_dat_df$pred <- predict_for_beta(Xdesign = Xdesign_val,
                                    beta = ECE_beta)


# calibration plot
list_calibration_plots[["ECE"]] <- calibrationPlot(data = val_dat_df,
                                                   caption = "ECE optimization")


# ROC
list_ROC[["ECE"]] <- ROCplot(data = val_dat_df,
                             caption = "ECE optimization")

# ECE / Brier score
list_ECE[["ECE"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["ECE"]] <- Brier(val_dat_df$Y, val_dat_df$pred)


# ECE incorrect model specification ---------------------------------------
ECE_results_IC_df <- data.frame(step = ece_results_IC[, 1],
                                ECE = ece_results_IC[, 2]) %>%
  mutate(beta0 = grid_search_IC[step, 1],
         beta1 = grid_search_IC[step, 2],
         beta2 = grid_search_IC[step, 3])

ECE_results_IC_df_long <- ECE_results_IC_df %>%
  gather(coef, value, -c(step, ECE))


# visualize results
ggplot(data = ECE_results_IC_df_long,
       aes(x = value,
           y = ECE,
           color = coef)) +
  geom_point() +
  geom_hline(aes(yintercept = ECE_train_ICML)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ coef)

# Find the best parameters (maximizing AUC)
best_ECE_IC <- ECE_results_IC_df %>% filter(ECE == min(ECE))


# best beta
ECE_beta_IC <- c(best_ECE_IC$beta0, best_ECE_IC$beta1, best_ECE_IC$beta2)

# predict for best beta
val_dat_df$pred <- predict_for_beta(Xdesign = Xdesign_val_IC,
                                    beta = ECE_beta_IC)

# calibration plot
list_calibration_plots[["ECE_IC"]] <- calibrationPlot(data = val_dat_df,
                                                      caption = "ECE optimization, incorrect")


# ROC
list_ROC[["ECE_IC"]] <- ROCplot(data = val_dat_df,
                                caption = "ECE optimization, incorrect")

# ECE / Brier score
list_ECE[["ECE_IC"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["ECE_IC"]] <- Brier(val_dat_df$Y, val_dat_df$pred)


# Brier correct model -----------------------------------------------------
brier_results <- read_cluster_results(method = "Brier")


Brier_results_df <- data.frame(step = brier_results[, 1],
                               Brier = brier_results[, 2]) %>%
  mutate(beta0 = grid_search[step, 1],
         beta1 = grid_search[step, 2],
         beta2 = grid_search[step, 3],
         beta12 = grid_search[step, 4])

Brier_results_df_long <- Brier_results_df %>%
  gather(coef, value, -c(step, Brier))


# visualize results
ggplot(data = Brier_results_df_long,
       aes(x = value,
           y = Brier,
           color = coef)) +
  geom_point() +
  geom_hline(aes(yintercept = Brier_train_ML)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ coef)

# Find the best parameters (maximizing AUC)
best_Brier <- Brier_results_df %>% filter(Brier == min(Brier))


# best beta
Brier_beta <- c(best_Brier$beta0, best_Brier$beta1,
                best_Brier$beta2, best_Brier$beta12)

# predict for best beta
val_dat_df$pred <- predict_for_beta(Xdesign = Xdesign_val,
                                    beta = Brier_beta)


# calibration plot
list_calibration_plots[["Brier"]] <- calibrationPlot(data = val_dat_df,
                                                     caption = "Brier score optimization")


# ROC
list_ROC[["Brier"]] <- ROCplot(data = val_dat_df,
                               caption = "Brier score optimization")

# ECE / Brier score
list_ECE[["Brier"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["Brier"]] <- Brier(val_dat_df$Y, val_dat_df$pred)



# Brier incorrect model specification -------------------------------------
Brier_results_IC_df <- data.frame(step = brier_results_IC[, 1],
                                  Brier = brier_results[, 2]) %>%
  mutate(beta0 = grid_search_IC[step, 1],
         beta1 = grid_search_IC[step, 2],
         beta2 = grid_search_IC[step, 3])

Brier_results_IC_df_long <- Brier_results_IC_df %>%
  gather(coef, value, -c(step, Brier))


# visualize results
ggplot(data = Brier_results_IC_df_long,
       aes(x = value,
           y = Brier,
           color = coef)) +
  geom_point() +
  geom_hline(aes(yintercept = Brier_train_ICML)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(color = "") +
  facet_grid(~ coef)

# Find the best parameters (maximizing AUC)
best_Brier_IC <- Brier_results_IC_df %>% filter(Brier == min(Brier))


# best beta
Brier_beta_IC <- c(best_Brier_IC$beta0, best_Brier_IC$beta1, best_Brier_IC$beta2, best_Brier_IC$beta12)

# predict for best beta
val_dat_df$pred <- predict_for_beta(Xdesign = Xdesign_val_IC,
                                    beta = Brier_beta_IC)


# calibration plot
list_calibration_plots[["Brier_IC"]] <- calibrationPlot(data = val_dat_df,
                                                        caption = "Brier score optimization, incorrect")


# ROC
list_ROC[["Brier_IC"]] <- ROCplot(data = val_dat_df,
                                  caption = "Brier score optimization, incorrect")

# ECE / Brier score
list_ECE[["Brier_IC"]] <- ECE(val_dat_df$Y, val_dat_df$pred)
list_Brier[["Brier_IC"]] <- Brier(val_dat_df$Y, val_dat_df$pred)



# Comparison --------------------------------------------------------------
# coefficient estimates for correct and incorrect model specification
coef_results <- data.frame(AUC = AUC_beta,
                           ECE = ECE_beta,
                           Brier = Brier_beta,
                           ML = ML_mod_coef[, 1],

                           AUC_IC = c(AUC_beta_IC, NA),
                           ECE_IC = c(ECE_beta_IC, NA),
                           Brier_IC = c(Brier_beta_IC, NA),
                           ML_IC = c(ML_mod_coef_IC[, 1], NA),

                           True = c(beta0, beta1, beta2, beta12),
                           Coef = c("(Intercept)", "X1", "X2", "X1:X2")) %>%
  gather(Method, Estimate, -c(Coef, True)) %>%
  mutate(is_ML = ifelse(Method == "ML" | Method == "ML_IC", "ML", "no-ML"))

# plot coefficient estimation for all methods
# for correct and incorrect model specification
list_coefficients[["overall"]] <- ggplot(data = coef_results,
                                               aes(x = Coef,
                                                   y = Estimate,
                                                   color = Method,
                                                   shape = is_ML)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  geom_point(aes(y = True,
                 color = "True coefficient"),
             shape = 4, size = 3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Model coefficient",
       color = "", shape = "",
       caption = "AUC / Brier / ECE with correct / incorrect model specification") +
  coord_flip()

# show coefficients and confidence interval for ML estimation
ggarrange(list_coefficients$ML_CMS,
          list_coefficients$ML_ICMS,
          # list_coefficients$AUC_ECE_Brier,
          nrow = 1, ncol = 2)

list_coefficients$AUC_ECE_Brier


# calibration plots
ggarrange(list_calibration_plots$ML_CMS,
          list_calibration_plots$ML_ICMS,
          list_calibration_plots$AUC,
          list_calibration_plots$ECE,
          list_calibration_plots$Brier,
          list_calibration_plots$AUC_IC,
          list_calibration_plots$ECE_IC,
          list_calibration_plots$Brier_IC,

          nrow = 3, ncol = 3)

# AUC
ggarrange(list_ROC$ML_CMS,
          list_ROC$ML_ICMS,
          list_ROC$AUC,
          list_ROC$ECE,
          list_ROC$Brier,
          list_ROC$AUC_IC,
          list_ROC$ECE_IC,
          list_ROC$Brier_IC)

# ECE
list_ECE %>% list.cbind()

# Brier score
list_Brier %>% list.cbind()
