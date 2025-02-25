### SIMULATION acc. ROUANET PAPER ###

library(tidyverse)
library(MASS)
library(rlist)

# packages for generalized linear mixed models
library(nlme)
library(lmerTest)
library(lme4)

# packages for GEE
library(gee)
library(geepack)
library(geeM)

# package for joint models
library(JM)


RouanetSim <- function(n, tij,
                       eta, gamma,
                       DNAR) {

  # covariate
  X <- rbinom(n = n, size = 1, prob = 0.5)

  Y.long <- array(NA, dim = c(n, length(tij)))

  # random effects
  u.i <- mvrnorm(n = n, mu = c(0, 0),
                 Sigma = matrix(c(0.3, -0.1, -0.1, 0.1), byrow = T,
                                nrow = 2))

  # longitudinal marker
  for (i in seq(1, n)) {
    Y.long[i, ] <- 20 - 0.3 * tij + u.i[i, 1] + u.i[i, 2] * tij + rnorm(1, mean = 0, sd = 0.9)
  }

  # non-terminal event / terminal event and censoring
  risk.T <- array(1, dim = c(n, length(tij)))
  risk.C <- array(1, dim = c(n, length(tij)))

  Y.T <- array(0, dim = c(n, length(tij)))
  Y.C <- array(0, dim = c(n, length(tij)))

  weight.C <- array(1, dim = c(n, length(tij)))
  weight.CT <- array(1, dim = c(n, length(tij)))


  for (i in seq(1, n)) {

    # generate time of event for DNAR scenario
    if (DNAR) {
      hazard <- gamma[1] * exp(gamma[2] * X[i] + gamma[3] * u.i[i, 2])
      EventTime <- -log(runif(1)) / hazard
    }

    for (t in seq_along(tij[-1]) + 1) {
      # current time
      time <- tij[t]

      # is patient currently at risk?
      it.risk.C <- risk.C[i, t - 1]
      it.risk.T <- risk.T[i, t - 1]

      # estimate linear predictor for event and censoring probability
      eventProb <- gamma[1] + gamma[2] * X[i] + gamma[3] * Y.long[i, t - 1]
      censProb <- eta[1] + eta[2] * X[i] + eta[3] * Y.long[i, t - 1]

      weight.C[i, t] <- exp(censProb) / (1 + exp(censProb))
      weight.CT[i, t] <- (exp(censProb) / (1 + exp(censProb))) * (exp(eventProb) / (1 + exp(eventProb)))

      # if currently at risk
      if (it.risk.T == 1) {
        if (DNAR == F) {
          # event this time; model defined in terms of survival probability
          Y.T[i, t] <- rbinom(n = 1, size = 1,
                              prob = 1 - exp(eventProb) / (1 + exp(eventProb)))
        } else {
          Y.T[i, t] <- (EventTime <= time)
        }

        # not at risk for event or censoring in future, once event happened
        if (Y.T[i, t] == 1) {
          risk.T[i, t] <- 0
          risk.C[i, t] <- 0
        }

        # if event not happened and at risk for censoring
        if (Y.T[i, t] == 0 & it.risk.C == 1) {
          # censoring this time; model defined in terms of survival probability
          Y.C[i, t] <- rbinom(n = 1, size = 1,
                               prob = 1 - exp(censProb) / (1 + exp(censProb)))
          if (Y.C[i, t] == 1) {
            risk.C[i, t] <- 0
          }
        } else if (it.risk.C == 0) { # if not at risk for censoring
          Y.C[i, t] <- 1
          risk.C[i, t] <- 0
        }

        # if not at risk for death
      } else if (it.risk.T == 0) {
        Y.T[i, t] <- 1

        if (Y.C[i, t - 1] == 1) {
          Y.C[i, t] <- 1
        }

        risk.C[i, t] <- 0
        risk.T[i, t] <- 0
      }
    }
  }

  X_df <- data.frame(X = X,
                     id = seq(1, n))

  Y.long.df <- Y.long %>%
    as.data.frame()
  colnames(Y.long.df) <- tij
  Y.long.df <- Y.long.df %>%
    mutate(id = seq(1, n)) %>%
    gather(key = "time",
           value = "Marker",
           -id)

  Y.T.df <- Y.T %>%
    as.data.frame()
  colnames(Y.T.df) <- tij
  Y.T.df <- Y.T.df %>%
    mutate(id = seq(1, n)) %>%
    gather(key = "time",
           value = "Event",
           -id)

  Y.C.df <- Y.C %>%
    as.data.frame()
  colnames(Y.C.df) <- tij
  Y.C.df <- Y.C.df %>%
    mutate(id = seq(1, n)) %>%
    gather(key = "time",
           value = "Cens",
           -id)

  weight.C.df <- weight.C %>%
    as.data.frame()
  colnames(weight.C.df) <- tij
  weight.C.df <- weight.C.df %>%
    mutate(id = seq(1, n)) %>%
    gather(key = "time",
           value = "Weight.C",
           -id)

  weight.CT.df <- weight.CT %>%
    as.data.frame()
  colnames(weight.CT.df) <- tij
  weight.CT.df <- weight.CT.df %>%
    mutate(id = seq(1, n)) %>%
    gather(key = "time",
           value = "Weight.CT",
           -id)




  inner_join(Y.long.df,
             Y.T.df,
             by = c("id", "time")) %>%
    inner_join(Y.C.df,
               by = c("id", "time")) %>%
    inner_join(X_df,
               by = c("id")) %>%
    inner_join(weight.C.df,
               by = c("time", "id")) %>%
    inner_join(weight.CT.df,
               by = c("time", "id")) %>%
    mutate(time = as.numeric(time)) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(Weight.CT = cumprod(Weight.CT),
           Weight.C = cumprod(Weight.C)) %>%
    return()
}


# save predictions
DARMAR.RES <- list()
DNARMAR.RES <- list()


j <- 1
r <- 1

while (r <= 500) {

  # Data simulation ---------------------------------------------------------
  # set seed
  set.seed(j)

  dataDARMAR <- RouanetSim(n = 500,
                           tij = c(0, 4, 12, 16, 20),
                           eta = c(-0.5, 0.5, 0.15),
                           gamma = c(0.5, 0.5, 0.1),
                           DNAR = F) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(MarkerLag = c(0, na.omit(lag(Marker)))) %>%
    ungroup()

  # observed data: observed as long as not truncated by event or censoring
  dataDARMAR.CT <- dataDARMAR %>%
    group_by(id) %>%
    mutate(CTInd = ifelse(Event == 1 | Cens == 1, 1, 0)) %>%
    filter(cumsum(CTInd) <= 1) %>%
    ungroup() %>%
    arrange(id, time) %>%
    as.data.frame()



  dataDNARMAR <- RouanetSim(n = 1000,
                            tij = c(0, 4, 12, 16, 20),
                            eta = c(-1, 1, 0.2),
                            gamma = c(0.05, -1, -1),
                            DNAR = T) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(MarkerLag = c(0, na.omit(lag(Marker)))) %>%
    ungroup()

  # observed data: observed as long as not truncated by event or censoring
  dataDNARMAR.CT <- dataDNARMAR %>%
    group_by(id) %>%
    mutate(CTInd = ifelse(Event == 1 | Cens == 1, 1, 0)) %>%
    filter(cumsum(CTInd) <= 1) %>%
    ungroup() %>%
    arrange(id, time) %>%
    as.data.frame()



  # Models ------------------------------------------------------------------
  # linear mixed models
  LMM.DARMAR <- lmer(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2) + (1 | id),
                     data = dataDARMAR.CT)

  LMM.DNARMAR <- lmer(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2) + (1 | id),
                      data = dataDNARMAR.CT)

  # IEE
  IEE.DARMAR <- geeglm(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                       id = id,
                       data = dataDARMAR.CT)
  IEE.DNARMAR <- geeglm(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                        id = id,
                        data = dataDNARMAR.CT)

  # pooled regression for weights
  CT.Weights.DARMAR <- glm(CTInd ~ X + MarkerLag,
                           data = dataDARMAR.CT,
                           family = binomial)
  dataDARMAR.CT <- dataDARMAR.CT %>%
    mutate(prob.CT = predict(CT.Weights.DARMAR,
                             type = "response")) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(weight.CT = cumprod(1 - prob.CT)) %>%
    ungroup()

  CT.Weights.DNARMAR <- glm(CTInd ~ X + MarkerLag,
                            data = dataDNARMAR.CT,
                            family = binomial)
  dataDNARMAR.CT <- dataDNARMAR.CT %>%
    mutate(prob.CT = predict(CT.Weights.DNARMAR,
                             type = "response")) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(weight.CT = cumprod(1 - prob.CT)) %>%
    ungroup()

  C.Weights.DARMAR <- glm(Cens ~ X + MarkerLag,
                          data = dataDARMAR.CT,
                          family = binomial)

  dataDARMAR.CT <- dataDARMAR.CT %>%
    mutate(prob.C = predict(C.Weights.DARMAR,
                            type = "response")) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(weight.C = cumprod(1 - prob.C)) %>%
    ungroup()

  C.Weights.DNARMAR <- glm(Cens ~ X + MarkerLag,
                           data = dataDNARMAR.CT,
                           family = binomial)
  dataDNARMAR.CT <- dataDNARMAR.CT %>%
    mutate(prob.C = predict(C.Weights.DNARMAR,
                            type = "response")) %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(weight.C = cumprod(1 - prob.C)) %>%
    ungroup()

  # GEE with drop-out and event probability
  GEE.CT.DARMAR <- geeglm(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                          id = id,
                          data = dataDARMAR.CT,
                          weights = 1 / weight.CT)
  GEE.CT.DNARMAR <- geeglm(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                           id = id,
                           data = dataDNARMAR.CT,
                           weights = 1 / weight.CT)


  # GEE with drop-out
  GEE.C.DARMAR <- geeglm(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                         id = id,
                         data = dataDARMAR.CT,
                         weights = 1 / weight.C)
  GEE.C.DNARMAR <- geeglm(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                          id = id,
                          data = dataDNARMAR.CT,
                          weights = 1 / weight.C)


  # Joint model
  LME.DNARMAR <- lme(Marker ~ X + time + I(time^2) + I(X * time) + I(X * time^2),
                     random = ~ 1 | id,
                     data = dataDNARMAR.CT,
                     method = "REML")

  dataDNARMAR.CT.UNQ <- rbind(dataDNARMAR.CT %>%
                                group_by(id) %>%
                                filter(all(Event == 0)) %>%
                                filter(time == max(time)),
                              dataDNARMAR.CT %>%
                                group_by(id) %>%
                                filter(Event == 1))

  SURV.DNARMAR <- coxph(Surv(time, Event) ~ X + MarkerLag,
                        data = dataDNARMAR.CT.UNQ,
                        x = TRUE)
  dForm <- list(fixed = ~ 1 + I(2 * time) + X + I(2 * X * time),
                indFixed = c(3, 4, 5, 6),
                random = ~ 1,
                indRandom = 1)

  # JM.DNARMAR <- jointModel(LME.DNARMAR, SURV.DNARMAR,
  #                          timeVar = "time",
  #                          parameterization = "slope",
  #                          derivForm = dForm)
  JM.DNARMAR <- tryCatch({
    JM.DNARMAR <- jointModel(LME.DNARMAR, SURV.DNARMAR,
                             timeVar = "time",
                             parameterization = "slope",
                             derivForm = dForm)
    JM.DNARMAR  # Return the successful result
  }, error = function(e) {
    message("Error encountered at index ", i, ": ", e$message)
    return(NULL)  # Return NULL and move to the next iteration
  })

  if (is.null(JM.DNARMAR)) {
    print(j)
    j <- j + 1
    next  # Skip to the next iteration of the loop
  }

  # Estimation --------------------------------------------------------------
  DARMAR.RES[[r]] <- data.frame(time = dataDARMAR$time,
                                X = dataDARMAR$X,

                                LMM = predict(LMM.DARMAR, newdata = dataDARMAR),
                                IEE = predict(IEE.DARMAR, newdata = dataDARMAR),
                                GEE.CT = predict(GEE.CT.DARMAR, newdata = dataDARMAR),
                                GEE.C = predict(GEE.C.DARMAR, newdata = dataDARMAR))

  DNARMAR.RES[[r]] <- data.frame(time = dataDNARMAR$time,
                                 X = dataDNARMAR$X,

                                 LMM = predict(LMM.DNARMAR, newdata = dataDNARMAR),
                                 IEE = predict(IEE.DNARMAR, newdata = dataDNARMAR),
                                 GEE.CT = predict(GEE.CT.DNARMAR, newdata = dataDNARMAR),
                                 GEE.C = predict(GEE.C.DNARMAR, newdata = dataDNARMAR),
                                 JM = predict(JM.DNARMAR, newdata = dataDNARMAR))

  # move index
  print(r)
  r <- r + 1

}



# Averaging ---------------------------------------------------------------
# aggregate over repetitions for each subject
patients.DARMAR <- DARMAR.RES %>%
  list.rbind() %>%
  group_by(time, X) %>%
  reframe(LMM = mean(LMM),
          IEE = mean(IEE),
          GEE.CT = mean(GEE.CT),
          GEE.C = mean(GEE.C))

patients.DNARMAR <- DNARMAR.RES %>%
  list.rbind() %>%
  group_by(time, X) %>%
  reframe(LMM = mean(LMM),
          IEE = mean(IEE),
          GEE.CT = mean(GEE.CT),
          GEE.C = mean(GEE.C),
          JM = mean(JM))


# Observed data
patients.DARMAR <- inner_join(patients.DARMAR,
                              dataDARMAR.CT %>%
                                group_by(time, X) %>%
                                reframe(mean.CT = mean(Marker)),
                              by = c("time", "X"))

patients.DNARMAR <- inner_join(patients.DNARMAR,
                               dataDNARMAR.CT %>%
                                 group_by(time, X) %>%
                                 reframe(mean.CT = mean(Marker)),
                               by = c("time", "X"))

# Alive cohort
patients.DARMAR <- inner_join(patients.DARMAR,
                              dataDARMAR %>%
                                group_by(id) %>%
                                filter(Event != 1) %>%
                                ungroup() %>%
                                group_by(time, X) %>%
                                reframe(mean.T = mean(Marker)),
                              by = c("time", "X"))

patients.DNARMAR <- inner_join(patients.DNARMAR,
                               dataDNARMAR.CT %>%
                                 group_by(id) %>%
                                 filter(Event != 1) %>%
                                 ungroup() %>%
                                 group_by(time, X) %>%
                                 reframe(mean.T = mean(Marker)),
                               by = c("time", "X"))


# Immortal cohort
patients.DARMAR <- inner_join(patients.DARMAR,
                              dataDARMAR %>%
                                group_by(time, X) %>%
                                reframe(mean = mean(Marker)),
                              by = c("time", "X"))

patients.DNARMAR <- inner_join(patients.DNARMAR,
                               dataDNARMAR %>%
                                 group_by(time, X) %>%
                                 reframe(mean = mean(Marker)),
                               by = c("time", "X"))

ggplot(data = patients.DARMAR,
       aes(x = time,
           group = X)) +
  geom_line(aes(y = LMM,
                color = "LMM",
                linetype = "LMM")) +
  geom_line(aes(y = IEE,
                color = "IEE",
                linetype = "IEE")) +
  geom_line(aes(y = GEE.CT,
                color = "GEE, adj. for drop-out and event",
                linetype = "GEE, adj. for drop-out and event")) +
  geom_line(aes(y = GEE.C,
                color = "GEE, adj. for drop-out",
                linetype = "GEE, adj. for drop-out")) +
  # geom_line(aes(y = JM,
  #               color = "JM",
  #               linetype = "JM")) +
  geom_point(aes(y = mean.CT,
                 color = "Mean traj. among observed subjects",
                 shape = "Mean traj. among observed subjects")) +
  geom_point(aes(y = mean.T,
                 color = "Mean traj. among survivors",
                 shape = "Mean traj. among survivors")) +
  geom_point(aes(y = mean,
                 color = "Mean traj. among immortal cohort",
                 shape = "Mean traj. among immortal cohort")) +
  facet_grid(~ X) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Time",
       y = "Marker",
       color = "", shape = "", linetype = "",
       title = "DAR-MAR") +
  scale_color_manual(values = c("LMM" = "red",
                                "IEE" = "blue",
                                "GEE, adj. for drop-out and event" = "green",
                                "GEE, adj. for drop-out" = "purple",
                                # "JM" = "orange",
                                "Mean traj. among observed subjects" = "black",
                                "Mean traj. among survivors" = "grey",
                                "Mean traj. among immortal cohort" = "brown")) +
  scale_shape_manual(values = c("Mean traj. among observed subjects" = 16,
                                "Mean traj. among survivors" = 17,
                                "Mean traj. among immortal cohort" = 18),
                     guide = "none") +
  scale_linetype_manual(values = c("LMM" = "solid",
                                   "IEE" = "dashed",
                                   "GEE, adj. for drop-out and event" = "dotted",
                                   "GEE, adj. for drop-out" = "dotdash"
                                   # "JM" = "twodash"
                                   ),
                        guide = "none")
ggsave("Simulation/Standard Methods/DAR-MAR.pdf",
       width = 21,
       height = 15,
       unit = "cm")




ggplot(data = patients.DNARMAR,
       aes(x = time,
           group = X)) +
  geom_line(aes(y = LMM,
                color = "LMM",
                linetype = "LMM")) +
  geom_line(aes(y = IEE,
                color = "IEE",
                linetype = "IEE")) +
  geom_line(aes(y = GEE.CT,
                color = "GEE, adj. for drop-out and event",
                linetype = "GEE, adj. for drop-out and event")) +
  geom_line(aes(y = GEE.C,
                color = "GEE, adj. for drop-out",
                linetype = "GEE, adj. for drop-out")) +
  geom_line(aes(y = JM,
                color = "JM",
                linetype = "JM")) +
  geom_point(aes(y = mean.CT,
                 color = "Mean traj. among observed subjects",
                 shape = "Mean traj. among observed subjects")) +
  geom_point(aes(y = mean.T,
                 color = "Mean traj. among survivors",
                 shape = "Mean traj. among survivors")) +
  geom_point(aes(y = mean,
                 color = "Mean traj. among immortal cohort",
                 shape = "Mean traj. among immortal cohort")) +
  facet_grid(~ X) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Time",
       y = "Marker",
       color = "", shape = "", linetype = "",
       title = "DNAR-MAR") +
  scale_color_manual(values = c("LMM" = "red",
                                "IEE" = "blue",
                                "GEE, adj. for drop-out and event" = "green",
                                "GEE, adj. for drop-out" = "purple",
                                "JM" = "orange",
                                "Mean traj. among observed subjects" = "black",
                                "Mean traj. among survivors" = "grey",
                                "Mean traj. among immortal cohort" = "brown")) +
  scale_shape_manual(values = c("Mean traj. among observed subjects" = 16,
                                "Mean traj. among survivors" = 17,
                                "Mean traj. among immortal cohort" = 18),
                     guide = "none") +
  scale_linetype_manual(values = c("LMM" = "solid",
                                   "IEE" = "dashed",
                                   "GEE, adj. for drop-out and event" = "dotted",
                                   "GEE, adj. for drop-out" = "dotdash",
                                   "JM" = "twodash"),
                        guide = "none")

ggsave("Simulation/Standard Methods/DNAR-MAR.pdf",
       width = 21,
       height = 15,
       unit = "cm")




