
# Semi-competing risk scenario --------------------------------------------

expit <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
log(x / (1 - x))
}

# compute likehood contributions for bivariate binary covariate
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

SimulateBivariateLongData <- function(n.sample = 100,
                                      times = seq(1, 100),
                                      # baseline hazard
                                      alpha.nt, alpha.t, lambda.or,
                                      # baseline mean
                                      zeta.long,
                                      # coefficients for predictor
                                      beta.nt,
                                      beta.ntr, beta.t,
                                      nu.or,
                                      eta.long,

                                      # individual frailty --> later
                                      frailty = FALSE,

                                      # settings
                                      biv.binary,
                                      long.trajectory,

                                      # force of mortality (multiplicative effect on death)
                                      force.mortality = 0,

                                      # standard deviation for longitudinal trajectory
                                      sd.long.trajectory = 1,

                                      # cens.poten.rate is not really the censrate
                                      cens.poten.rate = 0)
{
  # check input dimension for baseline hazard and baseline mean
  if (length(alpha.nt) != length(times)) {stop("alpha.nt should have the same length as times")}
  if (length(alpha.t) != length(times)) {stop("alpha.t should have the same length as times")}
  if (length(lambda.or) != length(times)) {stop("lambda.or should have the same length as times")}
  if (length(zeta.long) != length(times)) {stop("zeta.long should have the same length as times")}

  # check logical input
  if (!is.logical(frailty)) {stop("frailty must have type logical")}
  if (!is.logical(biv.binary)) {stop("biv.binary must have type logical")}
  if (!is.logical(long.trajectory)) {stop("long.trajectory must have type logical")}

  # multiplicative effect for force of mortality must be positive
  if (force.mortality < 0) {stop("force.mortality must be non-negative")}


  K <- length(times)
  p <- length(beta.nt)

  # draw covariates
  # time independent
  X.time.fixed <- as.matrix(rnorm(n.sample))
  # time dependent
  X.time.dep <- matrix(nrow = n.sample, ncol = K, 0)
  # at baseline k = 1, X_2 \sim Ber(0.6)
  X.time.dep[, 1] <- rbinom(n.sample, 1, 0.6)

  # create matrices to save risk sets and event indicators
  risk.NT <- risk.T <- YNT <- YT <- matrix(nrow = n.sample, ncol = K, 0)
  # create matrix for longitudinal covariate
  Y.traj <- matrix(nrow = n.sample, ncol = K, NA)
  Y.traj.all <- matrix(nrow = n.sample, ncol = K, 0)
  # everyone is at risk at baseline
  risk.NT[,1] <- 1
  risk.T[,1] <- 1

  # events / longitudinal covariate at baseline
  # design matrix
  Xfirst <- cbind(X.time.fixed, X.time.dep[,1])
  # longitudinal covariate
  Y.traj[, 1] <- as.matrix(rnorm(n = n.sample, mean = zeta.long[1], sd = sd.long.trajectory))
  Y.traj.all[, 1] <- Y.traj[, 1]

  # merge design matrix
  YXfirst <- cbind(Y.traj[, 1], Xfirst)
  # probability for non-terminal event
  p.nt.first <- expit(alpha.nt[1] + YXfirst%*%beta.nt)
  # probability for terminal event
  p.t.first <- force.mortality * expit(alpha.t[1] + YXfirst%*%beta.t)
  # OR
  OR.first <- exp(lambda.or[1]+ YXfirst%*%nu.or)

  # compute probabilities for each of four events that can happen (0, 0), (0, 1), (1, 0), (1, 1)
  first.probs <- MargORtoJoint(p1marg = p.nt.first,
                               p2marg = p.t.first,
                               OR = OR.first)

  # sample events
  crude.data <- apply(X = first.probs, MARGIN = 1,
                      FUN = sample,
                      x = seq(1, 4),
                      size = 1,
                      replace = T)
  # save indicator if non-terminal event happened (2, 4) or terminal-event happened (3, 4)
  YNT[,1] <- crude.data %in% c(2,4)
  YT[,1] <- crude.data %in% c(3,4)

  YNT[YNT[,1]==1, 2] <- 1
  YT[YT[,1]==1, 2] <- 1

  # run through all partitions
  for(k in seq(2, K)) {

    # time-dependent covariate
    Xnow <- X.time.dep[ , k - 1]
    Xnow[Xnow==1] <- rbinom(sum(Xnow), 1, 0.9)
    X.time.dep[ ,k] <- Xnow

    # design matrix
    Xtemp <- cbind(X.time.fixed, Xnow)

    # determine risk sets
    # for non-terminal event: neither non-terminal nor terminal event happened
    risk.NT[,k] <- (YNT[, k - 1]==0 & YT[, k - 1]==0)
    # for terminal event: terminal event happened
    risk.T[,k] <- (YT[, k - 1]==0)
    survivors <- (YT[, k - 1]==0)
    # at risk for terminal event only, non-terminal event already happened
    at.risk.T.only <- (risk.NT[, k]==0 & risk.T[, k]==1)
    # at risk for both events
    at.risk.both <- (risk.NT[, k]==1 & risk.T[, k]==1)

    # design matrix
    YXtemp <- cbind(Y.traj[, k-1], Xtemp)
    YXtemp.all <- cbind(Y.traj.all[, k-1], Xtemp)

    # longitudinal trajectory for survivors
    Y.traj[survivors, k] <- as.matrix(rnorm(n = sum(survivors),
                                   mean = zeta.long[k] + YXtemp[survivors, ]%*%eta.long,
                                   sd = sd.long.trajectory))
    # general longitudinal trajectory
    Y.traj.all[, k] <- as.matrix(rnorm(n = sum(n.sample),
                                       mean = zeta.long[k] + YXtemp.all%*%eta.long,
                                       sd = sd.long.trajectory))

    # design matrix
    YXtemp <- cbind(Y.traj[, k], Xtemp)

    # if patient at risk for terminal event only
    if (sum(at.risk.T.only) > 0) {
      probs.T.only <- force.mortality * expit(alpha.t[k] + YXtemp[at.risk.T.only, ]%*%beta.t + beta.ntr)
      YT[at.risk.T.only, k] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
    }
    # if patient at risk for both events
    if  (sum(at.risk.both) > 0) {
      # marginal probabilities
      p.nt.both <- expit(alpha.nt[k] + YXtemp[at.risk.both, ]%*%beta.nt)
      p.t.both <- force.mortality * expit(alpha.t[k] + YXtemp[at.risk.both, ]%*%beta.t)
      OR.both <- exp(lambda.or[k]+ YXtemp[at.risk.both, ]%*%nu.or)
      # joint probabilities
      probs.both <- MargORtoJoint(p1marg = p.nt.both, p2marg = p.t.both, OR = OR.both)
      # sample events
      crude.data.both <- apply(X = probs.both,
                               MARGIN = 1,
                               FUN = sample,
                               x = seq(1, 4),
                               size = 1, replace = T)
      # save indicator
      YNT[at.risk.both, k] <- crude.data.both %in% c(2, 4)
      YT[at.risk.both, k] <- crude.data.both %in% c(3, 4)
    }
    if(k < K)  {
      YNT[YNT[, k]==1, k + 1] <- 1
      YT[YT[, k]==1, k + 1] <- 1
    }
  }

  # add uninformative censoring (always observed at baseline)
  C <- sample(x = seq(2, K), replace = T, size = n.sample)
  somecens <- rbinom(n.sample, 1, cens.poten.rate)
  cens <- rep(0, n.sample)

  for (i in seq(1, n.sample)) {
    if (somecens[i]==1) {
      YNT[i, C[i]:K] <- YT[i, C[i]:K] <- Y.traj[i, C[i]:K] <- Y.traj.all[i, C[i]:K] <- NA
      risk.NT[i, C[i]:K] <- risk.T[i, C[i]:K] <- 0
      # generate censoring indicator (no death before)
      if (risk.T[i, C[i] - 1]==0) {cens[i] = 1}
    }
  }

  obs.per.person <- apply(risk.T,1, function(x) sum(x==1, na.rm = T))
  ID <- rep(seq(1, n.sample), times = obs.per.person)

  # turn data into long format
  Xcln <- matrix(ncol = 2, nrow = sum(obs.per.person))
  TM <- YNTcln <- YTcln <- TIMEcln <- vector(length = sum(obs.per.person))
  Y.traj.cln <- matrix(nrow = sum(obs.per.person), ncol = 1, NA)

  # trajectories
  Y.traj.all.cln <- matrix(nrow = n.sample * K, ncol = 4, NA)

  temp.ind <- 1
  for (i in seq(1, n.sample)) {
    nobs.i <- obs.per.person[i]
    indicesY <- seq(temp.ind, (temp.ind + nobs.i - 1))

    # event indicators
    YNTcln[indicesY] <- YNT[i, seq(1, nobs.i)]
    YTcln[indicesY] <- YT[i, seq(1, nobs.i)]
    # count of observations
    TM[indicesY] <- seq(1, nobs.i)
    # times
    TIMEcln[indicesY] <- times[seq(1, nobs.i)]

    # longitudinal trajectory
    Y.traj.cln[indicesY, 1] <- Y.traj[i, seq(1, nobs.i)]

    # design matrix
    Xcln[indicesY, 1] <- rep(X.time.fixed[i], nobs.i)
    Xcln[indicesY, 2] <- X.time.dep[i, seq(1, nobs.i)]

    # iterative indexing
    temp.ind <- temp.ind + nobs.i

    # trajectory for all, beyond death
    Y.traj.all.cln[seq((i-1)*K + 1, i*K), 1] <- i
    Y.traj.all.cln[seq((i-1)*K + 1, i*K), 2] <- YT[i, ]
    Y.traj.all.cln[seq((i-1)*K + 1, i*K), 3] <- Y.traj.all[i, ]
    Y.traj.all.cln[seq((i-1)*K + 1, i*K), 4] <- times
  }

  df.return <- data.frame(X = Xcln,
                          YNT = YNTcln,
                          YT = YTcln,
                          Ytraj = Y.traj.cln[, 1],
                          ID = ID,
                          TM = TM,
                          TIME = TIMEcln)
  return(list(df.return = df.return,
              cens = cens,
              Y.traj.all = data.frame(ID = Y.traj.all.cln[, 1],
                                      YT = Y.traj.all.cln[, 2],
                                      Ytraj = Y.traj.all.cln[, 3],
                                      times = Y.traj.all.cln[, 4]
                                      )
              )
         )
}









