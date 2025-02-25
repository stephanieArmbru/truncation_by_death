##### HERLPER #####


expit <- function(x) {
  exp(x) / (1 + exp(x))
}

logit <- function(x) {
  log(x / (1 - x))
}


# prediction - Landmarking style
predictLongitMTparam.Landmark <- function(par,
                                 longit.data,
                                 formula.long,
                                 formula.long.AR,
                                 formula.long.trend,

                                 formula.T,
                                 formula.glob.dep,
                                 formula.loc.dep) {

  ID <- longit.data$ID
  TM <- longit.data$TM

  n.sample <- ID %>% unique() %>% length()

  # Longitudinal model
  if (!is.null(formula.long)) {
    # [-1] to avoid intercept column
    Xlong <- as.matrix(model.matrix(formula.long, data = longit.data)[, -1])
    p.long <- ncol(Xlong)
  } else {
    Xlong <- NULL
    p.long <- 0
  }

  if(!is.null(formula.long.AR)) {
    # remove baseline values
    XlongAR <- as.matrix(model.matrix(formula.long.AR, data = longit.data)[, -1])
    # [-seq(1, n.sample), ]
    p.long.AR <- ncol(XlongAR)
  } else {
    XlongAR <- NULL
    p.long.AR <- 0
  }

  if (!is.null(formula.long.trend)) {
    XlongTrend <- as.matrix(model.matrix(formula.long.trend, data = longit.data))
    p.long.Trend <- ncol(XlongTrend)
  } else {
    XlongTrend <- NULL
    p.long.Trend <- 0
  }

  # Terminal event model
  if (!is.null(formula.T)) {
    XT <- as.matrix(model.matrix(formula.T, data = longit.data)[, -1])
    p.T <- ncol(XT)
  } else {
    XT <- NULL
    p.T <- 0
  }
  if (!is.null(formula.glob.dep)) {
    XTglobdep <- as.matrix(model.matrix(formula.glob.dep, data = longit.data)[, -1])
    # [-seq(1, n.sample), ]
    p.T.globdep <- ncol(XTglobdep)
  } else {
    XTglobdep <- NULL
    p.T.globdep <- 0
  }
  if (!is.null(formula.loc.dep)) {
    XTlocdep <- as.matrix(model.matrix(formula.loc.dep, data = longit.data)[, -1])
    # [-seq(1, n.sample), ]
    p.T.locdep <- ncol(XTlocdep)
  } else {
    XTlocdep <- NULL
    p.T.locdep <- 0
  }

  n.par.long <- p.long + p.long.AR + p.long.Trend + 1
  n.par.T <- 1 + p.T + p.T.globdep + p.T.locdep

  # Extract parameters
  zeta <- par[1:p.long.Trend]
  eta <- par[(p.long.Trend + 1):(p.long.Trend + p.long.AR)]
  beta <- par[(p.long.Trend + p.long.AR + 1):(p.long.Trend + p.long.AR + p.long)]
  # sigma.sq <- par[(p.long.Trend + p.long.AR + p.long + 1):n.par.long]
  lambda <- par[(n.par.long + 1)]
  theta <- par[(n.par.long + 2):((n.par.long + 1) + p.T.globdep)]
  varphi <- par[(n.par.long + p.T.globdep + 2):(n.par.long + p.T.globdep + 1 + p.T.locdep)]
  xi <- par[(n.par.long + 1 + p.T.globdep + p.T.locdep + 1):(n.par.long + n.par.T)]

  # Linear Predictors
  LinPred.long <- XlongTrend %*% zeta + XlongAR %*% eta + Xlong %*% beta
  LinPred.T <- lambda + XTglobdep %*% theta + XTlocdep %*% varphi + XT %*% xi

  Prob.T <- exp(LinPred.T) / (1 + exp(LinPred.T))

  return(list(MarkerPred = LinPred.long,
              PredProb = Prob.T))
}


predictLongitMTparam.Baseline <- function(par,
                                          # baseline covariate design matrices
                                          Xlong, XT,
                                          Ylong, YT,
                                          # number of parameters
                                          p.long.AR,
                                          p.long.Trend,
                                          p.T.globdep,
                                          p.T.locdep,

                                          # time and ID information
                                          ID, TM, time
                                          ) {

  K <- length(unique(TM))

  n.sample <- ID %>% unique() %>% length()

  p.long <- ncol(Xlong)
  p.T <- ncol(XT)

  n.par.long <- p.long + p.long.AR + p.long.Trend + 1
  n.par.T <- 1 + p.T + p.T.globdep + p.T.locdep

  # Extract parameters
  zeta <- par[1:p.long.Trend]
  eta <- par[(p.long.Trend + 1):(p.long.Trend + p.long.AR)]
  beta <- par[(p.long.Trend + p.long.AR + 1):(p.long.Trend + p.long.AR + p.long)]
  # sigma.sq <- par[(p.long.Trend + p.long.AR + p.long + 1):n.par.long]
  lambda <- par[(n.par.long + 1)]
  theta <- par[(n.par.long + 2):((n.par.long + 1) + p.T.globdep)]
  varphi <- par[(n.par.long + p.T.globdep + 2):(n.par.long + p.T.globdep + 1 + p.T.locdep)]
  xi <- par[(n.par.long + 1 + p.T.globdep + p.T.locdep + 1):(n.par.long + n.par.T)]


  Xbeta <- as.matrix(Xlong) %*% matrix(beta, ncol = 1)
  Xxi <- as.matrix(XT) %*% matrix(xi,  ncol = 1)

  Pred.long.K <- array(NA, dim = c(n.sample, K))
  idx <- ID[TM == 1]
  Pred.long.K[idx, 1] <- Ylong[TM == 1]
  idx <- ID[TM == 2]
  Pred.long.K[idx, 2] <- Ylong[TM == 2]
  Pred.T.K <- array(NA, dim = c(n.sample, K))
  idx <- ID[TM == 1]
  Pred.T.K[idx, 1] <- YT[TM == 1]
  idx <- ID[TM == 2]
  Pred.T.K[idx, 2] <- YT[TM == 2]

  alive.at.2 <- ID[TM == 2][(YT[TM == 2] == 0)]

  for (k.now in seq(3, K)) {
    Pred.long.K[alive.at.2, k.now] <- zeta[1] +
      zeta[2] * time[k.now] +
      zeta[3] * time[k.now]^2 +
      eta[1] * Pred.long.K[alive.at.2, k.now - 1] +
      eta[2] * Pred.long.K[alive.at.2, k.now - 2] +
      eta[3] * time[k.now] * Pred.long.K[alive.at.2, k.now - 1] +
      Xbeta[alive.at.2, ]

    LinPred.T <- lambda +
      theta[1] * Pred.long.K[alive.at.2, k.now - 1] +
      theta[2] * Pred.long.K[alive.at.2, k.now - 2] +
      varphi[1] * time[k.now] +
      varphi[2] * time[k.now]^2 +
      varphi[3] * time[k.now] * Pred.long.K[alive.at.2, k.now - 1] +
      Xxi[alive.at.2, ]

    Pred.T.K[alive.at.2, k.now] <- exp(LinPred.T) / (1 + exp(LinPred.T))
    }


  Pred.Long.df <- Pred.long.K %>%
    as.data.frame() %>%
    mutate(ID = seq(1, n.sample))
  colnames(Pred.Long.df) <- c(times, "ID")

  Pred.Long.df <- Pred.Long.df %>%
    gather("Time",
           "Marker",
           -ID)


  Pred.T.df <- Pred.T.K %>%
    as.data.frame() %>%
    mutate(ID = seq(1, n.sample))
  colnames(Pred.T.df) <- c(times, "ID")

  Pred.T.df <- Pred.T.df %>%
    gather("Time",
           "EventProb",
           -ID)

  Pred.df <- inner_join(Pred.Long.df,
                        Pred.T.df,
                        by = c("ID", "Time")) %>%
    mutate(Time = as.numeric(Time))

  return(Pred.df)
}


