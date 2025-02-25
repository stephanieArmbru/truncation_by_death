### AUTOREGRESSIVE EXPLORATION




SimDatARLongitTerminal <- function(n.sample = 100,
                                  # must include zero as start value
                                  times = seq(0, 100),
                                  # baseline covariates (no intercept)
                                  p = 10,

                                  # partition-specific probability of death
                                  lambda.t,

                                  # slope threshold for additional marker effect
                                  slope.t.threshold,
                                  # marker threshold of additional marker effect
                                  long.t.threshold,

                                  # global dependence
                                  theta.t,
                                  # partition-specific dependence
                                  varphi.t,
                                  # partition-specific effect for marker slope
                                  varphi.slope.t,
                                  # partition-specific effect for marker distance to threshold
                                  varphi.ReLU.t,
                                  # time-independent effect of baseline covariates
                                  xi.t,

                                  # partition-specific mean
                                  zeta.long,

                                  # partition-specific effect of autoregressive term
                                  eta.long,

                                  # time-independent effect of baseline covariates
                                  beta.long,

                                  # individual frailty --> later
                                  frailty = FALSE,

                                  # settings --> later; when modeling semi-comp + long
                                  biv.binary = FALSE,
                                  long.trajectory = TRUE,

                                  # force of mortality (multiplicative effect on death)
                                  force.mortality = 1,

                                  # standard deviation for longitudinal trajectory
                                  # for now, homogeneous over time and patients
                                  sd.long.trajectory = 1,

                                  # censoring
                                  # cens.poten.rate is not really the censrate --> later
                                  cens.poten.rate = 0)
{
  times_incl_zero <- times
  times <- times[-1]
  # check input dimension for parameters for death probability
  if (length(lambda.t) != length(times)) {stop("lambda.t should have length of times - 1")}
  if (length(varphi.t) != length(times)) {stop("varphi.t should have length of times - 1")}

  if (length(theta.t) != 1) {stop("theta.t should have dimension 1")}
  if(length(xi.t) != p) {stop("xi.t should have the same length as the number of covariates")}


  # check input dimensions for parameters for longitudinal marker
  if (length(zeta.long) != length(times_incl_zero)) {stop("zeta.long should have the same length as times")}

  if (length(eta.long) != length(times)) {stop("eta.long should have length of times - 1")}
  if (length(beta.long) != p) {stop("beta.long should have the same length as number of covariates")}


  # check logical input
  if (!is.logical(frailty)) {stop("frailty must have type logical")}
  if (!is.logical(biv.binary)) {stop("biv.binary must have type logical")}
  if (!is.logical(long.trajectory)) {stop("long.trajectory must have type logical")}

  # multiplicative effect for force of mortality must be positive
  if (force.mortality < 0) {stop("force.mortality must be non-negative")}
  # censoring rate must be positive
  if (cens.poten.rate < 0) {stop("cens.poten.rate must be non-negative")}

  # number of partition intervals
  K <- length(times)

  # baseline covariates: continuous Gaussian and discrete Binomial (coin toss)
  X.cont <- replicate(floor(p / 2),
                      {rnorm(n.sample, mean = 0, sd = 1)},
                      simplify = "array")

  X.bin <- replicate(ceiling(p / 2),
                     {rbinom(n = n.sample, size = 1, prob = 0.5)},
                     simplify = "array")

  # format design matrix
  X <- matrix(cbind(X.cont,
                    X.bin),
              nrow = n.sample, ncol = p,
              byrow = FALSE)

  # longitudinal covariate
  Y.traj <- array(dim = c(n.sample,
                          length(times_incl_zero)))

  # baseline longitudinal marker Y_0
  Y.traj[, 1] <- as.matrix(rnorm(n = n.sample,
                                 mean = zeta.long[1] + X %*% beta.long,
                                 sd = sd.long.trajectory))

  # mean value
  mu.k <- array(dim = c(n.sample,
                        length(times_incl_zero)))
  mu.k[, 1] <- zeta.long[1] + X %*% beta.long


  # death indicator array over time
  # baseline: all patients are at risk (I = 0)
  N.T <- array(0,
               dim = c(n.sample,
                       length(times_incl_zero)))
  # # risk indicator over time
  # at.risk.T <- array(0,
  #                    dim = c(n.sample,
  #                            length(times_incl_zero)))
  # # baseline: all patients are at risk (I = 0)
  # at.risk.T[, 1] <- 1

  # mortality probability array over time
  Prob.T <- array(0,
                  dim = c(n.sample,
                          length(times_incl_zero)))

  # time to terminal event
  Time.T <- array(max(times_incl_zero), dim = n.sample)

  # run through all partitions (+1 indexing due to zero time point)
  for(k in seq(2, K + 1)) {

    # who is at risk? (alive)
    at.risk.T <- N.T[, k - 1] == 0


    # death probability
    present.slope <- array(0, dim = c(n.sample, 1))
    if (k > 2) {
      present.slope[Y.traj[, k - 1] - Y.traj[, k - 2] <= slope.t.threshold, 1] <- 1
    }

    ## ReLU function: additive slope if marker value below a certain threshold
    ReLU.t <- array(0, dim = c(n.sample, 1))
    ReLU.t[Y.traj[, k - 1] <= long.t.threshold, 1] <- 1


    ## marker term
    ## VERSION 2:
    marker.k <- - ((theta.t +
                      varphi.t[k - 1]) * (Y.traj[, k-1])) +
      ReLU.t * (long.t.threshold - Y.traj[, k - 1]) * varphi.ReLU.t[k - 1] +
      present.slope * varphi.slope.t[k - 1]


    # linear predictor
    lin.pred.k <- lambda.t[k - 1] + marker.k + X %*% xi.t

    # avoid NA in probability
    probs.T <- expit( lin.pred.k )
    probs.T[lin.pred.k >= 5] <- 1
    probs.T[lin.pred.k <= -5] <- 0

    # save patient-specific probability for those at risk only
    # dead patients have probability zero to die again
    Prob.T[at.risk.T, k] <- probs.T[at.risk.T]

    # if people still at risk
    if (!is.na(sum(at.risk.T)) & sum(at.risk.T) > 0) {
      # draw binomial random variable for those still at risk (alive)
      N.T[at.risk.T, k] <- rbinom(sum(at.risk.T), 1, Prob.T[at.risk.T, k])

      # which patient experienced the terminal event in the current partition
      happened.T <- N.T[at.risk.T, k] == 1

      if (!is.na(sum(happened.T)) & (sum(happened.T) > 0)) {
        # impute time of death within interval
        Time.T[at.risk.T][happened.T] <- runif(sum(happened.T),
                                               min = times_incl_zero[k - 1],
                                               max = times_incl_zero[k])
      }
    }

    # compute mean for longitudinal trajectory

    ## VERSION 2:
    # linear predictor
    mu.k[, k] <- eta.long[k - 1] * mu.k[, k - 1] + zeta.long[k] + X %*% beta.long
    Y.traj[, k] <- as.matrix(rnorm(n = n.sample,
                                   mean = mu.k[, k],
                                   sd = sd.long.trajectory))

    # maintain 1 if ever one
    if (k < K + 1)  {
      N.T[N.T[, k] == 1, k + 1] <- 1
    }
  }

  # transform to data frame
  Y.traj.df <- Y.traj %>%
    as.data.frame() %>%
    mutate(ID = seq(1, n.sample))
  colnames(Y.traj.df) <- c(times_incl_zero, "ID")
  N.T.df <- N.T %>%
    as.data.frame() %>%
    mutate(ID = seq(1, n.sample))
  colnames(N.T.df) <- c(times_incl_zero, "ID")
  Prob.T.df <- Prob.T %>%
    as.data.frame() %>%
    mutate(ID = seq(1, n.sample))
  colnames(Prob.T.df) <- c(times_incl_zero, "ID")

  # transform to long format
  Y.traj.df.long <- Y.traj.df %>%
    gather(key = "Time",
           value = "Marker", -ID)

  N.T.df.long <- N.T.df %>%
    gather(key = "Time",
           value = "Death", -ID)

  Prob.T.df.long <- Prob.T.df %>%
    gather(key = "Time",
           value = "Mort_Prob", -ID)

  Time.T.df <- data.frame(ID = seq(1, n.sample),
                          Mort_Time = Time.T)

  X_df <- X %>%
    as.data.frame()
  colnames(X_df) <- paste0("X", seq(1, p))
  X_df$ID <- seq(1, n.sample)

  # for immortal cohort
  df.return.prelim <- inner_join(Y.traj.df.long,
                                 N.T.df.long,
                                 by = c("ID", "Time"))
  df.return.prelim <- inner_join(df.return.prelim,
                                 Prob.T.df.long,
                                 by = c("ID", "Time"))

  df.return.prelim <- inner_join(df.return.prelim,
                                 Time.T.df,
                                 by = "ID")
  df.return.IM <- inner_join(df.return.prelim,
                             X_df,
                             by = "ID") %>%
    mutate(Time = as.numeric(Time)) %>%
    group_by(ID) %>%
    mutate(TM = rank(Time)) %>%
    ungroup() %>%
    arrange(Time, ID)


  # for mortal cohort
  df.return.M <- df.return.IM %>%
    group_by(ID) %>%
    arrange(Time) %>%
    mutate(obs_post_death = cumsum(Death)) %>%
    filter(obs_post_death <= 1) %>%
    ungroup() %>%
    arrange(Time, ID)





  return(list(df.return.IM = df.return.IM,
              df.return.M = df.return.M,

              # formatting for estimation method
              X = X_df

              # censoring
              # cens = cens
  )
  )
}



