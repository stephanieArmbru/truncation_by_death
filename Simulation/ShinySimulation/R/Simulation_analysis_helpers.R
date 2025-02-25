######## Simulation: Functions to help analyze simulation results ########

## Author: Stephanie Armbruster
## Date: 12/16/2024


# Visualize ---------------------------------------------------------------
# plot baseline death probability
plot_baseline_mortality_risk <- function(times,
                                         lambda.t) {
  ggplot(data = data.frame(times = times[-1],
                           lambda.t = lambda.t),
         aes(x = times,
             y = expit(lambda.t))) +
    geom_point() +
    theme_bw() +
    labs(y = "",
         title = "Baseline mortality (partition-specific)",
         x = "Time") %>%
    return()
}

# plot longitudinal mean
plot_longitudinal_trend <- function(times, zeta.long) {
  ggplot(data = data.frame(times = times,
                           mean = zeta.long),
         aes(x = times,
             y = mean)) +
    geom_point() +
    theme_bw() +
    labs(y = "",
         title = "Baseline mean longitudinal marker trajectory",
         x = "Time") %>%
    return()
}


# Heterogeneity for a chosen subset of patients over time
plot_heterogeneity <- function(data,
                               sample_IDs = NULL) {

  if (is.null(sample_IDs)) {
    sample_IDs <- sample(x = seq(1, length(unique(data$ID))),
                         size = 10)
  }

  ID_rank_t0 <- data %>%
    dplyr::filter(Time == 0) %>%
    pull(Marker) %>%
    rank()


  ggplot(data = data %>%
           dplyr::filter(ID %in% sample_IDs),
         aes(x = Time,
             y = Marker,
             color = ID %>% as.factor())) +
    geom_point() +
    geom_line() +
    geom_point(data = data %>%
                 group_by(ID) %>%
                 filter(ID %in% sample_IDs & Death == 1) %>%
                 filter(Time == min(Time)),
               aes(x = Time,
                   y = Marker,
                   color = "Death"),
               shape = 4,
               size = 3) +
    # geom_smooth() +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(color = "",
         x = "Time") %>%
    return()
}

# Patient-specific mortality probabilities over time; conditional on partition-specific survival
plot_mortality_prob <- function(data,
                                sample_IDs = NULL) {
  if (is.null(sample_IDs)) {
    sample_IDs <- sample(x = seq(1, length(unique(data$ID))),
                         size = 10)
  }

  ggplot(data = data %>% filter(ID %in% sample_IDs),
         aes(x = Time,
             y = Mort_Prob,
             color = ID %>% as.factor())) +
    geom_point() +
    geom_line() +
    geom_point(data = data %>%
                 group_by(ID) %>%
                 filter(ID %in% sample_IDs & Death == 1) %>%
                 filter(Time == min(Time)),
               aes(x = Time,
                   y = Mort_Prob,
                   color = "Death"),
               shape = 4,
               size = 3) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(color = "",
         x = "Time",
         y = "Conditional Death Probability") %>%
    return()
}

# Cumulative patient-specific mortality probability
plot_cum_mortality_prob <- function(data,
                                    sample_IDs = NULL) {

  if (is.null(sample_IDs)) {
    sample_IDs <- sample(x = seq(1, length(unique(data$ID))),
                         size = 10)
  }

  data_cum <- data %>%
    filter(ID %in% sample_IDs) %>%
    group_by(ID) %>%
    arrange(Time) %>%
    mutate(lag_Mort_Prob = c(0, na.omit(lag(Mort_Prob))),
           Surv = 1 - lag_Mort_Prob,
           cum_Surv = cumprod(Surv),
           marg_Mort_Prob = cum_Surv * Mort_Prob,
           cum_Mort_Prob = cumsum(marg_Mort_Prob))


  ggplot(data = data_cum,
         aes(x = Time,
             y = cum_Mort_Prob,
             color = ID %>% as.factor())) +
    geom_point() +
    geom_line() +
    # might through warning message when no death happened, but that is okay
    geom_point(data = data_cum %>%
                 group_by(ID) %>%
                 filter(ID %in% sample_IDs & Death == 1) %>%
                 filter(Time == min(Time)),
               aes(x = Time,
                   y = cum_Mort_Prob,
                   color = "Death"),
               shape = 4,
               size = 3) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(color = "",
         x = "Time",
         y = "Cumulative Death Probability") %>%
    return()
}



plot_patients_at_risk <- function(data) {
  # patients at risk
  n.risk <- data %>%
    group_by(Time) %>%
    reframe(n_risk = sum(Death == 0)) %>%
    dplyr::select(Time, n_risk)

  ggplot(data = n.risk,
         aes(x = Time,
             y = n_risk)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    labs(y = "Patients at risk",
         x = "Time") %>%
    return()
}

plot_immortal_mortal_longitudinal_trajectory <- function(data) {
  # Longitudinal marker
  # mean for all (beyond death)
  mean.Y.traj.all <- data %>%
    group_by(Time) %>%
    reframe(Time = unique(Time),
            mean_all = mean(na.omit(Marker))) %>%
    filter(Time != 0)
  # mean for survivors only
  mean.Y.traj.surv <- data %>%
    dplyr::filter(Death == 0) %>%
    group_by(Time) %>%
    reframe(Time = unique(Time),
            mean_surv = mean(na.omit(Marker))) %>%
    filter(Time != 0)

  ggplot() +
    geom_line(data = mean.Y.traj.all,
              aes(x = Time,
                  y = mean_all,
                  color = "Immortal cohort")) +
    geom_line(data = mean.Y.traj.surv,
              aes(x = Time,
                  y = mean_surv,
                  color = "Mortal cohort")) +
    # geom_smooth(method='lm', formula= y~x) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = "Longitudinal Marker",
         color = "",
         x = "Time") %>%
    return()
}

plot_diff_moral_immortal <- function(data) {
  # Longitudinal marker
  # mean for all (beyond death)
  mean.Y.traj.all <- data %>%
    group_by(Time) %>%
    reframe(Time = unique(Time),
            mean_all = mean(na.omit(Marker)))
  # mean for survivors only
  mean.Y.traj.surv <- data %>%
    dplyr::filter(Death == 0) %>%
    group_by(Time) %>%
    reframe(Time = unique(Time),
            mean_surv = mean(na.omit(Marker)))

  diff.mean.Y.traj <- inner_join(mean.Y.traj.all, mean.Y.traj.surv, by = "Time") %>%
    mutate(diff = mean_all - mean_surv) %>%
    filter(Time != 0)

  ggplot() +
    geom_line(data = diff.mean.Y.traj,
              aes(x = Time,
                  y = diff)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = "Difference in longitudinal Marker",
         color = "",
         x = "Time") %>%
    return()
}

plot_calibration <- function(data) {

  data_cum <- data %>%
    group_by(ID) %>%
    dplyr::mutate(lag_Mort_Prob = c(0, na.omit(dplyr::lag(Mort_Prob))),
           Surv_Prob = 1 - lag_Mort_Prob,
           cum_Surv = cumprod(Surv_Prob),
           marg_Mort_Prob = cum_Surv * Mort_Prob,
           cum_Mort_Prob = cumsum(marg_Mort_Prob)) %>%
    filter(Death == 1) %>%
    filter(Time == min(Time))

  ggplot(data = data_cum,
         aes(x = cum_Surv,
             y = Death)) +
    # geom_point(alpha = 0.1) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_point() +
    theme_bw() +
    labs(x = "Cumulative survival probability",
         y = "Death indicator") +
    theme(legend.position = "bottom")
}




# Calibration -------------------------------------------------------------
D_calibration <- function(data, B = 10) {

  data.id <- data %>%
    group_by(ID) %>%
    arrange(Time) %>%
    dplyr::mutate(lag_Mort_Prob = c(0, na.omit(lag(Mort_Prob))),
           Surv = 1 - lag_Mort_Prob,
           cum_Surv = cumprod(Surv)) %>%
    filter(Death == 1) %>%
    filter(Time == min(Time))

  quantiles <- 1 - seq(0,1,length.out = B + 1)
  int_B <- seq(0, 1, length.out = B + 1)

  binIndex <- seq_along(quantiles)[-1]

  uncensoredProbabilities <- data.id %>%
    # only uncensored patients
    filter(Death == 1) %>%
    pull(cum_Surv)

  uncensoredBinning <- array(0, B)
  for (b in seq(1, B)) {

    uncensoredBinning[B - b + 1] <- (((data.id$cum_Surv) >= int_B[b]) & ((data.id$cum_Surv) < int_B[b + 1])) %>%
      sum()
  }

  censoredProbabilities <- data.id %>%
    # only censored patients
    filter(Death == 0) %>%
    pull(cum_Surv)

  if(length(censoredProbabilities) > 0){
    censoredBinPositions <- prodlim::sindex(quantiles, censoredProbabilities,
                                  comp = "greater", strict = T)
    #Sometimes the probability will be 1 in which case we just want to put them in the first bin.
    censoredBinPositions <- ifelse(censoredBinPositions == 0, 1,
                                  censoredBinPositions)
    quantileWidth <- 1 / B
    firstBin <- ifelse(censoredBinPositions == B, 1,
                      (censoredProbabilities - quantiles[censoredBinPositions + 1]) / censoredProbabilities)
    restOfBins <- ifelse(censoredProbabilities == 0, 1, 1 / (B * censoredProbabilities))

    listOfContributions <- lapply(seq_along(censoredBinPositions),
                                  function(x) c(rep(0, censoredBinPositions[x] - 1),
                                                rep(firstBin[x],1),
                                                rep(restOfBins[x], B - censoredBinPositions[x])))
    censoredBinning <- colSums(plyr::ldply(listOfContributions, rbind))
  } else censoredBinning <- 0

  combinedBins <- uncensoredBinning + censoredBinning
  names(combinedBins) <- quantiles[-length(quantiles)]

  # p-value for chi squared test
  pvalue <- chisq.test(combinedBins)$p.value

  combinedBins_df <- data.frame(N = combinedBins) %>%
    mutate(Prop = (N / sum(N)) %>% round(digits = 4)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(p = pvalue)

  return(list(contingency_table = combinedBins_df,
              p = pvalue))
}



# Mortal and immortal mean ------------------------------------------------
MorImmorMean <- function(data) {
  mean.Y.traj.all <- data %>%
    group_by(Time) %>%
    reframe(Time = unique(Time),
            mean_all = mean(na.omit(Marker))) %>%
    filter(Time != 0)
  # mean for survivors only
  mean.Y.traj.surv <- data %>%
    dplyr::filter(Death == 0) %>%
    group_by(Time) %>%
    reframe(Time = unique(Time),
            mean_surv = mean(na.omit(Marker))) %>%
    filter(Time != 0)

  inner_join(mean.Y.traj.all, mean.Y.traj.surv,
             by = "Time") %>%
    return()
}



# Prediction --------------------------------------------------------------
# function to predict longitudinal marker trajectory for a new patient
# based on baseline value only
# iteratively using marker predictions over time
predict_baseline <- function(mod.coef.fixed,
                             newdata) {

  # available information at baseline
  data.t0 <- data %>%
    filter(Time == 0)
  # sample size
  n.sample <- data.t0 %>% nrow()

  # design matrix
  X.design <- data.t0[, grepl("X", colnames(data))] %>%
    as.matrix()

  # format coefficient data frame
  mod.coef.fixed.df <- data.frame(coef = names(mod.coef.fixed),
                                  value = mod.coef.fixed)

  # Intercept
  intercept <- mod.coef.fixed.df %>%
    filter(grepl("Intercept", coef)) %>%
    pull(value)

  # time-independent baseline covariate coefficients
  # while making sure they are in the right order
  beta <- mod.coef.fixed.df %>%
    filter(grepl("X", coef)) %>%
    mutate(cov_numb = as.numeric(sub("X", "", coef))) %>%
    arrange(cov_numb) %>%
    pull(value)

  # time dummy
  time.trends <- mod.coef.fixed.df %>%
    filter(grepl("Time", coef)) %>%
    mutate(Time = as.numeric(sub("Time_Dummy", "", coef)))

  # lagged marker effect
  marker.effect <- mod.coef.fixed.df %>%
    filter(grepl("lagged_Marker", coef)) %>%
    pull(value)

  # discrete time partition to predict
  times <- time.trends$Time %>% unique()

  # to save predictions
  pred.marker <- array(NA, dim = c(n.sample, length(times) + 1))

  # marker value at baseline
  pred.marker[, 1] <- data.t0$Marker

  for (t in seq_along(times)) {
    time <- times[t]

    time.trend <- time.trends %>%
      filter(Time == time) %>%
      pull(value)

    pred.marker[, t + 1] <- intercept + X.design %*% beta + time.trend + marker.effect * pred.marker[, t]
  }

  pred.marker %>%
    as.data.frame() %>%
    gather("Time",
           "Pred") %>%
    pull(Pred) %>%
    return()
}





# Stationarity ------------------------------------------------------------
estimate_ts_covariance <- function(data,
                                   k = 1
                                   ) {
  # cov_over_time <- data %>%
  #   group_by(ID) %>%
  #   arrange(Time) %>%
  #   mutate(Marker_lag = lead(Marker, k)) %>%
  #   ungroup() %>%
  #   group_by(Time) %>%
  #   reframe(cov = cov(Marker, Marker_lag))

cov_over_time <- data %>%
    group_by(ID) %>%
    reframe(cov = acf(Marker,
                      type = "covariance",
                      plot = FALSE)$acf[k, 1, 1])


  ggplot(data = cov_over_time,
         aes(x = ID,
             y = cov)) +
    geom_point() +
    theme_bw() %>%
    return()
}
