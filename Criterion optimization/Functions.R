######################### FUNCTIONS LOGISTIC SIMULATION ########################

## Stephanie Armbruster
## 10/20/2024


## This script defined functions require for the simulation of alternative
## estimation metric for a logistics regression based on optimizing criteria

# Functions ---------------------------------------------------------------
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
  log(x / (1 - x))
}

simDat <- function(N, beta0, beta1, beta2, beta12) {
  # Gaussian covariates
  X1 <- rnorm(n = N, mean = 0, sd = 1)
  X2 <- rnorm(n = N, mean = 0, sd = 1)

  # linear predictor
  lin_pred <- beta0 + X1 * beta1 + X2 * beta2 + X1 * X2 * beta12

  # outcome
  Y <- rbinom(n = N, size = 1, prob = expit(lin_pred))


  list(X = cbind(X1, X2),
       lin_pred = lin_pred,
       Y = Y
  ) %>%
    return()
}

# calibration plot
calibrationPlot <- function(data, caption) {

  ggplot(data = data,
         aes(x = pred,
             y = Y)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, bs = "cs")) +
    theme_bw() +
    labs(x = "Predicted probabilities",
         y = "Observation",
         caption = caption)
}

ROCplot <- function(data, caption) {
  roc_temp <- roc(data$Y, data$pred)


  ggplot(data = data.frame(sens = roc_temp$sensitivities,
                           spec = roc_temp$specificities),
         aes(y = sens,
             x = 1 - spec)) +
    geom_line() +
    theme_bw() +
    labs(y = "Sensitivity",
         x = "1 - Specificity",
         title = paste("AUC = ", round(roc_temp$auc, 4)),
         caption = caption) +
    geom_abline(aes(intercept = 0,
                    slope = 1),
                linetype = "dashed",
                color = "grey") %>%
    return()
}

ECE <- function(obs, pred) {
  mean(abs(obs - pred))
}

Brier <- function(obs, pred) {
  mean((obs - pred)^2)
}

predict_for_beta <- function(Xdesign, beta) {
  # linear predictor
  lin_pred <- Xdesign %*% beta

  # prediction
  expit(lin_pred)[, 1] %>%
    return()
}

# function to optimize AUC
AUC_optimization <- function(data, step,
                             grid_search,
                             Xdesign) {
  beta <- grid_search[step, ]

  # compute prediction
  lin_pred <- Xdesign %*% beta
  pred <- expit(lin_pred)[, 1]

  # estimate AUC
  roc(data$Y, pred)$auc %>%
    return()
}

ECE_optimization <- function(data, step,
                             grid_search,
                             Xdesign) {
  beta <- grid_search[step, ]

  # compute prediction
  lin_pred <- Xdesign %*% beta
  pred <- expit(lin_pred)[, 1]

  # estimate AUC
  ECE(obs = data$Y, pred = pred) %>%
    return()
}


Brier_optimization <- function(data, step,
                               grid_search,
                               Xdesign) {
  beta <- grid_search[step, ]

  # compute prediction
  lin_pred <- Xdesign %*% beta
  pred <- expit(lin_pred)[, 1]

  # estimate AUC
  Brier(obs = data$Y, pred = pred) %>%
    return()
}
