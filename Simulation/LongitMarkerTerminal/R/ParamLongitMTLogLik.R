# Log likelihood for parametric longitudinal discrete time model for longitudinal
# marker and terminal event;
# low-dimensional example simulation


#' @title The log likelihood for a longitudinal discrete time model for a longitudinal marker and a terminal event data with time-fixed covariates
#' and using unstructured parametric autoregressive function.
#' @param par A vector of pareter estimates, in the ordering marker pareters, marker variance and terminal event pareters.
#' @param ID A vector indicating to which patient a particular row belongs.
#' @param TM A vector indicating in which partition interval a particular row belongs (integer counting of intervals).
#' @param Xlong A matrix with baseline covariates (column-wise) for marker model. Each row denotes a particular subject.
#' @param XlongAR A matrix with autoregressive components (column-wise) for marker model; each row denotes a measurement for a particular subject at a particular time.
#' @param XlongTrend A matrix with time components for marker trend.
#' @param XT A matrix with baseline covariates (column-wise) for terminal event model. Each row denotes a particular subject.
#' @param XTglobdep A matrix with marker history for global dependence between marker and terminal event; each row denotes a measurement of a particular subject at a particular time.
#' @param XTlocdep A matrix with marker history for local dependence between marker and terminal event; each row denotes a measurement of a particular subject at a particular time.
#' @param Ylong A vector of marker values; each row denotes a particular subject at a particular time.
#' @param YT A vector of terminal event indicators; each row denotes a particular subject at a particular time.
#' @return The function returns an object of class \code{LongitMT} including estimates and confidence intervals for
#' the time-varying functions and coefficients.
#' @note Used in \code{LongitMTparam.R}.
#' @author Stephanie Armbruster
#' @export


ParamLongitMTLogLik <- function(par,
                                ID,
                                TM,
                                Xlong,
                                XlongAR,
                                XlongTrend,
                                XT,
                                XTglobdep,
                                XTlocdep,
                                Ylong,
                                YT

                                # sigma.sq
                                ) {

  # Unique times and subjects
  times <- unique(TM)
  K <- length(times) - 1
  ID.uqe <- unique(ID)
  n.sample <- length(ID.uqe)

  # Parameter dimensions
  p.long <- ncol(Xlong)
  p.long.AR <- ncol(XlongAR)
  p.long.Trend <- ncol(XlongTrend)
  p.T <- ncol(XT)
  p.glob.dep <- ncol(XTglobdep)
  p.loc.dep <- ncol(XTlocdep)

  n.par.long <- p.long + p.long.AR + p.long.Trend + 1
  n.par.T <- 1 + p.T + p.glob.dep + p.loc.dep

  # Extract parameters
  zeta <- par[1:p.long.Trend]
  eta <- par[(p.long.Trend + 1):(p.long.Trend + p.long.AR)]
  beta <- par[(p.long.Trend + p.long.AR + 1):(p.long.Trend + p.long.AR + p.long)]
  sigma.sq <- par[n.par.long]
  lambda <- par[(n.par.long + 1)]
  theta <- par[(n.par.long + 2):((n.par.long + 1) + p.glob.dep)]
  varphi <- par[(n.par.long + p.glob.dep + 2):(n.par.long + p.glob.dep + 1 + p.loc.dep)]
  xi <- par[(n.par.long + 1 + p.glob.dep + p.loc.dep + 1):(n.par.long + n.par.T)]

  # Linear Predictors
  LinPred.long <- XlongTrend %*% zeta + XlongAR %*% eta + Xlong %*% beta
  LinPred.T <- lambda + XTglobdep %*% theta + XTlocdep %*% varphi + XT %*% xi

  Prob.T <- exp(LinPred.T) / (1 + exp(LinPred.T))

  # Compute log-likelihood using vectorized operations
  LogLik <- (1 - YT) * ((-(1 / 2) * log(2 * pi * sigma.sq)) +
    (-1 / (2 * sigma.sq)) * (Ylong - LinPred.long)^2) +
    (1 - YT) * log(1 - Prob.T) +
    YT * log(Prob.T)

  if (any(LinPred.T <= -7)) {
    correct.ID <- LinPred.T <= -7
    LogLik[correct.ID] <- (1 - YT[correct.ID]) * ((-(1 / 2) * log(2 * pi * sigma.sq)) +
                            (-1 / (2 * sigma.sq)) * (Ylong[correct.ID] - LinPred.long[correct.ID])^2) +
      (1 - YT[correct.ID]) * log(1 - 0.0001) +
      YT[correct.ID] * log(0.0001)
  }

  if (any(LinPred.T >= 17)) {
    correct.ID <- LinPred.T >= 17
    LogLik[correct.ID] <- (1 - YT[correct.ID]) * ((-(1 / 2) * log(2 * pi * sigma.sq)) +
                                                    (-1 / (2 * sigma.sq)) * (Ylong[correct.ID] - LinPred.long[correct.ID])^2) +
      (1 - YT[correct.ID]) * log(1 - 0.9999) +
      YT[correct.ID] * log(0.9999)
  }

  if (is.nan(-sum(LogLik))) {
    LogLik <- -3871745
  }

  return(-sum(LogLik))  # Negative log-likelihood for minimization
}



##### ARCHIVE #####

# ParamLongitMTLogLik <- function(par,
#
#                                 ID,
#                                 TM,
#
#                                 # input
#                                 Xlong,
#                                 XlongAR,
#                                 XlongTrend,
#
#                                 XT,
#                                 XTglobdep,
#                                 XTlocdep,
#                                 # output
#                                 Ylong,
#                                 YT,
#                                 riskT
#
#                                 # sigma.sq = 0.25^2
# ) {
#
#   # times and ID
#   times <- unique(TM)
#   K <- length(times) - 1
#
#   ID.uqe <- unique(ID)
#   n.sample <- length(ID.uqe)
#
#   # parameter
#   p.long <- ncol(Xlong)
#   p.long.AR <- ncol(XlongAR)
#   p.long.trend <- ncol(XlongTrend)
#
#   p.T <- ncol(XT)
#   p.glob.dep <- ncol(XTglobdep)
#   p.loc.dep <- ncol(XTlocdep)
#
#   n.par.long <- p.long + p.long.AR + p.long.Trend + 1
#   n.par.T <- 1 + p.T + p.T.globdep + p.T.locdep
#
#   # longitudinal trend
#   zeta <- par[1 : p.long.Trend]
#   # AR component
#   eta <- par[(p.long.Trend + 1) : (p.long.Trend + p.long.AR)]
#   # covariate coefficients
#   beta <- par[(p.long.Trend + p.long.AR + 1) : (p.long.Trend + p.long.AR + p.long)]
#
#   # variance
#   sigma.sq <- par[(p.long.Trend + p.long.AR + p.long + 1) : (n.par.long)]
#
#   # baseline event rate
#   lambda <- par[(n.par.long + 1)]
#   # global dependence
#   theta <- par[(n.par.long + 2) : ((n.par.long + 1) + p.glob.dep)]
#   # local dependence
#   varphi <- par[(n.par.long + p.glob.dep + 2) : (n.par.long + p.glob.dep + 1 + p.loc.dep)]
#   # covariate coefficients
#   xi <- par[(n.par.long + 1 + p.glob.dep + p.loc.dep + 1) : (n.par.long + n.par.T)]
#
#   # baseline covariates
#   XBeta.long <- Xlong %*% beta
#   XBeta.T <- XT %*% xi
#
#   # history for global dependence
#   YTheta.T <- XTglobdep %*% theta
#
#   # save elements
#   LogLik <- array(0, dim = n.sample)
#
#   for (i in seq_len(n.sample)) {
#     # indices for pooled data
#     IDX <- ID == ID.uqe[i]
#     # -1 without baseline; -1 to obtain interval index
#     iTM <- TM[IDX][-1] - 1
#
#     # condition on baseline values
#     iY.long <- Ylong[IDX][-1]
#     iY.T <- YT[IDX][-1]
#
#     # condition on baseline values
#     iXlongAR <- XlongAR[IDX, ][-1, ]
#     iXTtrend <- XlongTrend[IDX, ][-1, ]
#     iYTheta.T <- YTheta.T[IDX][-1]
#     iXTlocdep <- XTlocdep[IDX, ][-1, ]
#
#
#     iLinPred.long <- iXTtrend %*% zeta +
#       XBeta.long[ID.uqe[i]] +
#       apply(as.matrix(iXlongAR * eta, drop = F),
#             MARGIN = 1,
#             FUN = function(x) sum(x))
#
#     if (is(tryCatch(iLinPred.long, error=function(e) e, warning=function(w) w),
#         "warning")) {print(i)}
#
#     iLinPred.T <- lambda +
#       XBeta.T[ID.uqe[i]] +
#       iYTheta.T +
#       apply(as.matrix(iXTlocdep * varphi, drop = F),
#             MARGIN = 1,
#             FUN = function(x) sum(x))
#
#
#     iLogLik <- (1 - iY.T) * (-1/2) * log(2 * pi * sigma.sq) +
#       (-1/(2 * sigma.sq)) * (iY.long - iLinPred.long)^2 +
#       (1 - iY.T) * (1 - exp(iLinPred.T) / (1 + exp(iLinPred.T))) +
#       iY.T * exp(iLinPred.T) / (1 + exp(iLinPred.T))
#
#     # save summed contribution
#     LogLik[i] <- sum(iLogLik)
#   }
#
#   return(-sum(LogLik))
# }
#
#
#
