######## Simulation: Data Generating Mechanism ########

## Author: Stephanie Armbruster
## Date: 12/10/2024



# Helper ------------------------------------------------------------------
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
log(x / (1 - x))
}

# compute likelihood contributions for bivariate binary covariate
# based on marginal likelihood and OR
MargORtoJoint <- function(p1marg, p2marg, OR) {
  # check for independence assumption with tolerance
  # both events happen
  p12 <- ifelse(OR > 0.999 & OR < 1.001,
                p1marg*p2marg,
                (1 + (p1marg + p2marg)*(OR - 1) -
                   sqrt((1 + (p1marg + p2marg)*(OR - 1))^2 -
                          4*OR*(OR -1)*p1marg*p2marg))/(2*(OR - 1))
                )
  # non-terminal event only happens
  p1 <- p1marg - p12
  # terminal event only happens
  p2 <- p2marg - p12
  # no event happens
  p0 <- 1 - p12 - p1 - p2

  prob <- cbind(p0, p1, p2, p12)

  # sanity check probabilities
  if (any(prob<0) | any(prob>1)) {
    stop(paste("Problems with probabilities. prob =", prob))  }
  return(prob)
}


# Longitudinal Marker + Death ---------------------------------------------

## This function simulates a longitudinal marker and mortality according to
## a discrete time partition (user-defined).
## For now, it simulates p-1 time-independent covariates
## (Gaussian, fair Bernoulli).
## The Baseline covariates have a time-independent effect on marker and
## death probability.
## It includes a factor by which the force of mortality is manipulated.
## The longitudinal marker and death probability are globally and locally
## associated.
## It includes uninformative censoring at a user-defined rate.

SimulateJointLongData <- function(n.sample = 100,
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

                                  # threshold for longitudinal marker
                                  long.threshold,

                                  # partition-specific mean
                                  zeta.long,
                                  # add-on for threshold
                                  zeta.ReLU.long,

                                  # partition-specific effect of autoregressive term
                                  eta.long,
                                  # slope threshold
                                  slope.threshold,
                                  # partition-specific effect of autoregressive slope,
                                  eta.slope.long,

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


  # check thresholds
  if (length(long.threshold) != 1) {stop("long.threshold should have dimension 1")}
  if (length(long.t.threshold) != 1) {stop("long.t.threshold should have dimension 1")}
  if (length(slope.t.threshold) != 1) {stop("slope.t.threshold should have dimension 1")}


  # check input dimensions for parameters for longitudinal marker
  if (length(zeta.long) != length(times_incl_zero)) {stop("zeta.long should have the same length as times")}
  if (length(zeta.ReLU.long) != length(times_incl_zero)) {stop("zeta.ReLU.long should have the same length as times")}

  if (length(eta.long) != length(times)) {stop("eta.long should have length of times - 1")}
  if (length(eta.slope.long) != length(times)) {stop("eta.slope.long should have length of times - 1")}
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
  X_cont <- replicate(floor(p / 2),
                      {rnorm(n.sample, mean = 0, sd = 1)},
                      simplify = "array")

  X_bin <- replicate(ceiling(p / 2),
                     {rbinom(n = n.sample, size = 1, prob = 0.5)},
                     simplify = "array")

  # format design matrix
  X <- matrix(cbind(X_cont,
                    X_bin),
              nrow = n.sample, ncol = p,
              byrow = FALSE)

  # longitudinal covariate
  Y.traj <- array(dim = c(n.sample,
                          length(times_incl_zero)))

  # baseline longitudinal marker
  Y.traj[, 1] <- as.matrix(rnorm(n = n.sample,
                                 mean = zeta.long[1] + X %*% beta.long,
                                 sd = sd.long.trajectory))

  # risk indicator array over time
  # baseline: all patients are at risk (I = 0)
  N.T <- array(0,
               dim = c(n.sample,
                      length(times_incl_zero)))

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

    # compute mean for longitudinal trajectory
    ## depends on prior slope
    prior.slope <- array(0, dim = c(n.sample, 1))
    if (k > 2) {
      prior.slope[Y.traj[, k - 1] - Y.traj[, k - 2] <= slope.threshold, 1] <- 1
    }

    ## ReLU function: additional intercept if prior Y value lies below certain threshold
    ReLU <- array(0, dim = c(n.sample, 1))
    ReLU[Y.traj[, k - 1] <= long.threshold, 1] <- 1

    ## intercept term
    mu_k_intercept <- zeta.long[k] + zeta.ReLU.long[k] * ReLU


    ## autoregressive term
    mu_k_autoregressive <- eta.long[k - 1] * (Y.traj[, k - 1] - X %*% beta.long - zeta.long[k - 1]) + eta.slope.long[k - 1] * prior.slope

    # piece it all together
    mu_k <- mu_k_intercept + X %*% beta.long + mu_k_autoregressive

    # longitudinal trajectory
    Y.traj[, k] <- as.matrix(rnorm(n = n.sample,
                                   mean = mu_k,
                                   sd = sd.long.trajectory))

    # death probability (different indexing since no one can die at time zero)
    ## depends on current slope in marker trajectory
    present.slope <- array(0, dim = c(n.sample, 1))
    present.slope[Y.traj[, k] - Y.traj[, k - 1] <= slope.t.threshold, 1] <- 1

    ## ReLU function: additive slope if marker value below a certain threshold
    ReLU.t <- array(0, dim = c(n.sample, 1))
    ReLU.t[Y.traj[, k] <= long.t.threshold, 1] <- 1

    ## intercept term
    baseline_k <- lambda.t[k - 1]

    ## marker term
    marker_k <- (theta.t +
                   varphi.t[k - 1]) * sign(max(zeta.long) - Y.traj[, k]) * (max(zeta.long) - Y.traj[, k])^2 +
      ReLU.t * (long.t.threshold - Y.traj[, k])^2 * varphi.ReLU.t[k - 1] +
      present.slope * varphi.slope.t[k - 1]

    probs.T <- expit( baseline_k + marker_k + X %*% xi.t )

    # save patient-specific probability
    Prob.T[, k] <- probs.T

    # if people still at risk
    if (sum(at.risk.T) > 0) {
      # draw binomial random variable for those still at risk (alive)
      N.T[at.risk.T, k] <- rbinom(sum(at.risk.T), 1, probs.T)

      # which patient experienced the terminal event in the current partition
      happened.T <- N.T[at.risk.T, k] == 1

      if (sum(happened.T) > 0) {
        # impute time of death within interval
        Time.T[at.risk.T][happened.T] <- runif(sum(happened.T),
                                               min = times_incl_zero[k - 1],
                                               max = times_incl_zero[k])
      }
    }

    # maintain 1 if ever one
    if (k < K + 1)  {
      N.T[N.T[, k] == 1, k + 1] <- 1
    }
  }


  # add uninformative censoring (always observed at baseline)
  C <- sample(x = seq(2, K), replace = T, size = n.sample)
  somecens <- rbinom(n.sample, 1, cens.poten.rate)
  cens <- rep(0, n.sample)

  for (i in seq(1, n.sample)) {
    if (somecens[i]==1) {
      Time.T[i] <- min(C[i], Time.T[i])
      N.T[i, C[i]:K] <- NA
      Y.traj[i, C[i]:K] <- NA
      # generate censoring indicator (no death before)
      if (N.T[i, C[i] - 1] == 0) {cens[i] = 1}
    }
  }

  # transform to data frame
  Y.traj.df <- Y.traj %>% as.data.frame()
  colnames(Y.traj.df) <- times_incl_zero
  N.T.df <- N.T %>% as.data.frame()
  colnames(N.T.df) <- times_incl_zero
  Prob.T.df <- Prob.T %>% as.data.frame()
  colnames(Prob.T.df) <- times_incl_zero

  # transform to long format
  Y.traj.df.long <- Y.traj.df %>%
    gather(key = "Time",
           value = "Marker") %>%
    mutate(ID = rep(seq(1, n.sample), length(times_incl_zero)))

  N.T.df.long <- N.T.df %>%
    gather(key = "Time",
           value = "Death") %>%
    mutate(ID = rep(seq(1, n.sample), length(times_incl_zero)))

  Prob.T.df.long <- Prob.T.df %>%
    gather(key = "Time",
           value = "Mort_Prob") %>%
    mutate(ID = rep(seq(1, n.sample), length(times_incl_zero)))

  Time.T.df <- data.frame(ID = seq(1, n.sample),
                          Mort_Time = Time.T)


  X_df <- X %>%
    as.data.frame()
  colnames(X_df) <- paste0("X", seq(1, p))
  X_df$ID <- seq(1, n.sample)

  # merge all components
  df.return.prelim <- inner_join(Y.traj.df.long,
                                 N.T.df.long,
                                 by = c("ID", "Time"))
  df.return.prelim <- inner_join(df.return.prelim,
                                 Prob.T.df.long,
                                 by = c("ID", "Time"))

  df.return.prelim <- inner_join(df.return.prelim,
                                 Time.T.df,
                                 by = "ID")
  df.return <- inner_join(df.return.prelim,
                          X_df,
                          by = "ID") %>%
    mutate(Time = as.numeric(Time))

  return(list(df.return = df.return,
              cens = cens)
  )
}










# Semi-Competing risks + Longitudinal Marker ------------------------------

# SimulateBivariateLongData <- function(n.sample = 100,
#                                       times = seq(1, 100),
#                                       # baseline hazard
#                                       alpha.nt, alpha.t, lambda.or,
#                                       # baseline mean
#                                       zeta.long,
#                                       # coefficients for predictor
#                                       beta.nt,
#                                       beta.ntr, beta.t,
#                                       nu.or,
#                                       eta.long,
#
#                                       # individual frailty --> later
#                                       frailty = FALSE,
#
#                                       # settings
#                                       biv.binary,
#                                       long.trajectory,
#
#                                       # force of mortality (multiplicative effect on death)
#                                       force.mortality = 0,
#
#                                       # standard deviation for longitudinal trajectory
#                                       sd.long.trajectory = 1,
#
#                                       # cens.poten.rate is not really the censrate
#                                       cens.poten.rate = 0)
# {
#   # check input dimension for baseline hazard and baseline mean
#   if (length(alpha.nt) != length(times)) {stop("alpha.nt should have the same length as times")}
#   if (length(alpha.t) != length(times)) {stop("alpha.t should have the same length as times")}
#   if (length(lambda.or) != length(times)) {stop("lambda.or should have the same length as times")}
#   if (length(zeta.long) != length(times)) {stop("zeta.long should have the same length as times")}
#
#   # check logical input
#   if (!is.logical(frailty)) {stop("frailty must have type logical")}
#   if (!is.logical(biv.binary)) {stop("biv.binary must have type logical")}
#   if (!is.logical(long.trajectory)) {stop("long.trajectory must have type logical")}
#
#   # multiplicative effect for force of mortality must be positive
#   if (force.mortality < 0) {stop("force.mortality must be non-negative")}
#
#
#   K <- length(times)
#   p <- length(beta.nt)
#
#   # draw covariates
#   # time independent
#   X.time.fixed <- as.matrix(rnorm(n.sample))
#   # time dependent
#   X.time.dep <- matrix(nrow = n.sample, ncol = K, 0)
#   # at baseline k = 1, X_2 \sim Ber(0.6)
#   X.time.dep[, 1] <- rbinom(n.sample, 1, 0.6)
#
#   # create matrices to save risk sets and event indicators
#   risk.NT <- risk.T <- YNT <- YT <- matrix(nrow = n.sample, ncol = K, 0)
#   # create matrix for longitudinal covariate
#   Y.traj <- matrix(nrow = n.sample, ncol = K, NA)
#   Y.traj.all <- matrix(nrow = n.sample, ncol = K, 0)
#   # everyone is at risk at baseline
#   risk.NT[,1] <- 1
#   risk.T[,1] <- 1
#
#   # events / longitudinal covariate at baseline
#   # design matrix
#   Xfirst <- cbind(X.time.fixed, X.time.dep[,1])
#   # longitudinal covariate
#   Y.traj[, 1] <- as.matrix(rnorm(n = n.sample, mean = zeta.long[1], sd = sd.long.trajectory))
#   Y.traj.all[, 1] <- Y.traj[, 1]
#
#   # merge design matrix
#   YXfirst <- cbind(Y.traj[, 1], Xfirst)
#   # probability for non-terminal event
#   p.nt.first <- expit(alpha.nt[1] + YXfirst%*%beta.nt)
#   # probability for terminal event
#   p.t.first <- force.mortality * expit(alpha.t[1] + YXfirst%*%beta.t)
#   # OR
#   OR.first <- exp(lambda.or[1]+ YXfirst%*%nu.or)
#
#   # compute probabilities for each of four events that can happen (0, 0), (0, 1), (1, 0), (1, 1)
#   first.probs <- MargORtoJoint(p1marg = p.nt.first,
#                                p2marg = p.t.first,
#                                OR = OR.first)
#
#   # sample events
#   crude.data <- apply(X = first.probs, MARGIN = 1,
#                       FUN = sample,
#                       x = seq(1, 4),
#                       size = 1,
#                       replace = T)
#   # save indicator if non-terminal event happened (2, 4) or terminal-event happened (3, 4)
#   YNT[,1] <- crude.data %in% c(2,4)
#   YT[,1] <- crude.data %in% c(3,4)
#
#   YNT[YNT[,1]==1, 2] <- 1
#   YT[YT[,1]==1, 2] <- 1
#
#   # run through all partitions
#   for(k in seq(2, K)) {
#
#     # time-dependent covariate
#     Xnow <- X.time.dep[ , k - 1]
#     Xnow[Xnow==1] <- rbinom(sum(Xnow), 1, 0.9)
#     X.time.dep[ ,k] <- Xnow
#
#     # design matrix
#     Xtemp <- cbind(X.time.fixed, Xnow)
#
#     # determine risk sets
#     # for non-terminal event: neither non-terminal nor terminal event happened
#     risk.NT[,k] <- (YNT[, k - 1]==0 & YT[, k - 1]==0)
#     # for terminal event: terminal event happened
#     risk.T[,k] <- (YT[, k - 1]==0)
#     survivors <- (YT[, k - 1]==0)
#     # at risk for terminal event only, non-terminal event already happened
#     at.risk.T.only <- (risk.NT[, k]==0 & risk.T[, k]==1)
#     # at risk for both events
#     at.risk.both <- (risk.NT[, k]==1 & risk.T[, k]==1)
#
#     # design matrix
#     YXtemp <- cbind(Y.traj[, k-1], Xtemp)
#     YXtemp.all <- cbind(Y.traj.all[, k-1], Xtemp)
#
#     # longitudinal trajectory for survivors
#     Y.traj[survivors, k] <- as.matrix(rnorm(n = sum(survivors),
#                                    mean = zeta.long[k] + YXtemp[survivors, ]%*%eta.long,
#                                    sd = sd.long.trajectory))
#     # general longitudinal trajectory
#     Y.traj.all[, k] <- as.matrix(rnorm(n = sum(n.sample),
#                                        mean = zeta.long[k] + YXtemp.all%*%eta.long,
#                                        sd = sd.long.trajectory))
#
#     # design matrix
#     YXtemp <- cbind(Y.traj[, k], Xtemp)
#
#     # if patient at risk for terminal event only
#     if (sum(at.risk.T.only) > 0) {
#       probs.T.only <- force.mortality * expit(alpha.t[k] + YXtemp[at.risk.T.only, ]%*%beta.t + beta.ntr)
#       YT[at.risk.T.only, k] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
#     }
#     # if patient at risk for both events
#     if  (sum(at.risk.both) > 0) {
#       # marginal probabilities
#       p.nt.both <- expit(alpha.nt[k] + YXtemp[at.risk.both, ]%*%beta.nt)
#       p.t.both <- force.mortality * expit(alpha.t[k] + YXtemp[at.risk.both, ]%*%beta.t)
#       OR.both <- exp(lambda.or[k]+ YXtemp[at.risk.both, ]%*%nu.or)
#       # joint probabilities
#       probs.both <- MargORtoJoint(p1marg = p.nt.both, p2marg = p.t.both, OR = OR.both)
#       # sample events
#       crude.data.both <- apply(X = probs.both,
#                                MARGIN = 1,
#                                FUN = sample,
#                                x = seq(1, 4),
#                                size = 1, replace = T)
#       # save indicator
#       YNT[at.risk.both, k] <- crude.data.both %in% c(2, 4)
#       YT[at.risk.both, k] <- crude.data.both %in% c(3, 4)
#     }
#     if(k < K)  {
#       YNT[YNT[, k]==1, k + 1] <- 1
#       YT[YT[, k]==1, k + 1] <- 1
#     }
#   }
#
#   # add uninformative censoring (always observed at baseline)
#   C <- sample(x = seq(2, K), replace = T, size = n.sample)
#   somecens <- rbinom(n.sample, 1, cens.poten.rate)
#   cens <- rep(0, n.sample)
#
#   for (i in seq(1, n.sample)) {
#     if (somecens[i]==1) {
#       YNT[i, C[i]:K] <- YT[i, C[i]:K] <- Y.traj[i, C[i]:K] <- Y.traj.all[i, C[i]:K] <- NA
#       risk.NT[i, C[i]:K] <- risk.T[i, C[i]:K] <- 0
#       # generate censoring indicator (no death before)
#       if (risk.T[i, C[i] - 1]==0) {cens[i] = 1}
#     }
#   }
#
#   obs.per.person <- apply(risk.T,1, function(x) sum(x==1, na.rm = T))
#   ID <- rep(seq(1, n.sample), times = obs.per.person)
#
#   # turn data into long format
#   Xcln <- matrix(ncol = 2, nrow = sum(obs.per.person))
#   TM <- YNTcln <- YTcln <- TIMEcln <- vector(length = sum(obs.per.person))
#   Y.traj.cln <- matrix(nrow = sum(obs.per.person), ncol = 1, NA)
#
#   # trajectories
#   Y.traj.all.cln <- matrix(nrow = n.sample * K, ncol = 4, NA)
#
#   temp.ind <- 1
#   for (i in seq(1, n.sample)) {
#     nobs.i <- obs.per.person[i]
#     indicesY <- seq(temp.ind, (temp.ind + nobs.i - 1))
#
#     # event indicators
#     YNTcln[indicesY] <- YNT[i, seq(1, nobs.i)]
#     YTcln[indicesY] <- YT[i, seq(1, nobs.i)]
#     # count of observations
#     TM[indicesY] <- seq(1, nobs.i)
#     # times
#     TIMEcln[indicesY] <- times[seq(1, nobs.i)]
#
#     # longitudinal trajectory
#     Y.traj.cln[indicesY, 1] <- Y.traj[i, seq(1, nobs.i)]
#
#     # design matrix
#     Xcln[indicesY, 1] <- rep(X.time.fixed[i], nobs.i)
#     Xcln[indicesY, 2] <- X.time.dep[i, seq(1, nobs.i)]
#
#     # iterative indexing
#     temp.ind <- temp.ind + nobs.i
#
#     # trajectory for all, beyond death
#     Y.traj.all.cln[seq((i-1)*K + 1, i*K), 1] <- i
#     Y.traj.all.cln[seq((i-1)*K + 1, i*K), 2] <- YT[i, ]
#     Y.traj.all.cln[seq((i-1)*K + 1, i*K), 3] <- Y.traj.all[i, ]
#     Y.traj.all.cln[seq((i-1)*K + 1, i*K), 4] <- times
#   }
#
#   df.return <- data.frame(X = Xcln,
#                           YNT = YNTcln,
#                           YT = YTcln,
#                           Ytraj = Y.traj.cln[, 1],
#                           ID = ID,
#                           TM = TM,
#                           TIME = TIMEcln)
#   return(list(df.return = df.return,
#               cens = cens,
#               Y.traj.all = data.frame(ID = Y.traj.all.cln[, 1],
#                                       YT = Y.traj.all.cln[, 2],
#                                       Ytraj = Y.traj.all.cln[, 3],
#                                       times = Y.traj.all.cln[, 4]
#                                       )
#               )
#          )
# }









