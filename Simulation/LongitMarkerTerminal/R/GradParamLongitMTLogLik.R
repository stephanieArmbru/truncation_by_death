# Gradient of log likelihood for parametric longitudinal discrete time model for longitudinal
# marker and terminal event;
# low-dimensional example simulation


#' @title The gradient log likelihood for a longitudinal discrete time model for a longitudinal marker and a terminal event data with time-fixed covariates
#' and using unstructured partition-specific autoregressive function.
#' @param par A vector of parameter estimates, in the ordering marker parameters, marker variance and terminal event parameters.
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
#' @note Used in \code{LongitMTpw.R}.
#' @author Stephanie Armbruster
#' @export


GradParamLongitMTLogLik <- function(par,
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
  # Extract unique times and IDs
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
  # sigma.sq <- par[(p.long.Trend + p.long.AR + p.long + 1):n.par.long]
  lambda <- par[(n.par.long + 1)]
  theta <- par[(n.par.long + 2):((n.par.long + 1) + p.glob.dep)]
  varphi <- par[(n.par.long + p.glob.dep + 2):(n.par.long + p.glob.dep + 1 + p.loc.dep)]
  xi <- par[(n.par.long + 1 + p.glob.dep + p.loc.dep + 1):(n.par.long + n.par.T)]

  # Precompute terms
  XBeta.Long <- Xlong %*% beta
  XBeta.T <- XT %*% xi
  YTheta.T <- XTglobdep %*% theta

  # Linear Predictors
  LinPred.long <- XlongTrend %*% zeta + XlongAR %*% eta + Xlong %*% beta
  LinPred.T <- lambda + XTglobdep %*% theta + XTlocdep %*% varphi + XT %*% xi

  Prob.T <- exp(LinPred.T) / (1 + exp(LinPred.T))


  # variance dependent factor
  VarFac <- 1 / (sigma.sq)

  GradLogLik <- array(0, dim = n.par.long + n.par.T)

  GradLogLik[1:p.long.Trend] <- VarFac * ((1 - YT) * (Ylong - LinPred.long))[, 1] %*% XlongTrend
  GradLogLik[(p.long.Trend + 1):(p.long.Trend + p.long.AR)] <- VarFac * ((1 - YT) * (Ylong - LinPred.long))[, 1] %*% XlongAR
  GradLogLik[(p.long.Trend + p.long.AR + 1):(p.long.Trend + p.long.AR + p.long)] <- VarFac * ((1 - YT) * (Ylong - LinPred.long))[, 1] %*% Xlong
  # GradLogLik[(p.long.Trend + p.long.AR + p.long + 1):n.par.long] <- -VarFac / 2 + 1 / (2 * sigma.sq^2) * sum((1 - YT) * (Ylong - LinPred.long)^2)
  GradLogLik[(n.par.long + 1)] <- sum(YT - Prob.T)
  GradLogLik[(n.par.long + 2):((n.par.long + 1) + p.glob.dep)] <- (YT - Prob.T)[, 1] %*% XTglobdep
  GradLogLik[(n.par.long + p.glob.dep + 2):(n.par.long + p.glob.dep + 1 + p.loc.dep)] <- (YT - Prob.T)[, 1] %*% XTlocdep
  GradLogLik[(n.par.long + 1 + p.glob.dep + p.loc.dep + 1):(n.par.long + n.par.T)] <- (YT - Prob.T)[, 1] %*% XT

  return(-GradLogLik)
}


##### ARCHIVE #####
# GradParamLongitMTLogLik <- function(par,
#
#                                  ID,
#                                  TM,
#
#                                  # input
#                                  Xlong,
#                                  XlongAR,
#                                  XlongTrend,
#
#                                  XT,
#                                  XTglobdep,
#                                  XTlocdep,
#                                  # output
#                                  Ylong,
#                                  YT,
#                                  riskT
#
#                                  # sigma.sq = 0.25^2
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
#   p.long.Trend <- ncol(XlongTrend)
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
#
#   # baseline covariates
#   XBeta.Long <- Xlong %*% beta
#   XBeta.T <- XT %*% xi
#
#   # history for global dependence
#   YTheta.T <- XTglobdep %*% theta
#
#   # save elements
#   GradLogLik <- array(0, dim = c(n.sample, n.par.long + n.par.T))
#
#   for (i in seq_len(n.sample)) {
#     # indices for pooled data
#     IDX <- ID == ID.uqe[i]
#     # -1 without baseline; -1 to obtain interval index
#     iTM <- TM[IDX][-1] - 1
#
#     # person-time contribution to partition intervals
#     ik.l <- min(iTM)
#     ik.r <- max(iTM)
#
#     # condition on baseline values
#     iY.long <- Ylong[IDX][-1]
#     iY.T <- YT[IDX][-1]
#     irisk.T <- riskT[IDX][-1]
#
#     # condition on baseline values
#     iXlongAR <- matrix(XlongAR[IDX, ][-1, ], ncol = ncol(XlongAR))
#     iXlongTrend <- matrix(XlongTrend[IDX, ][-1, ], ncol = ncol(XlongTrend))
#     iYTheta.T <- matrix(YTheta.T[IDX][-1], ncol = ncol(YTheta.T))
#     iXTlocdep <- matrix(XTlocdep[IDX, ][-1, ], ncol = ncol(XTlocdep))
#     iXTglobdep <- matrix(XTglobdep[IDX, ][-1, ], ncol = ncol(XTglobdep))
#
#
#     iLinPred.long <- iXlongTrend %*% zeta +
#       XBeta.long[ID.uqe[i]] +
#       apply(as.matrix(iXlongAR * eta, drop = F),
#             MARGIN = 1,
#             FUN = function(x) sum(x))
#
#     iLinPred.T <- lambda +
#       XBeta.T[ID.uqe[i]] +
#       iYTheta.T +
#       apply(as.matrix(iXTlocdep * varphi, drop = F),
#             MARGIN = 1,
#             FUN = function(x) sum(x))
#
#     iProb.T <- exp(iLinPred.T) / (1 + exp(iLinPred.T))
#
#
#     iGradLogLik <- array(0, dim = c(n.par.long + n.par.T,
#                                     length(unique(iTM))))
#
#     VarFac <- 1 / (2 * sigma.sq)
#
#     for (k.now in unique(iTM)) {
#
#       IDX.now <- iTM == k.now
#
#       iY.long.now <- iY.long[IDX.now]
#       iY.T.now <- iY.T[IDX.now]
#
#       iXlongAR.now <- matrix(iXlongAR[IDX.now, ],
#                              ncol = ncol(iXlongAR))
#       iXlongTrend.now <- matrix(iXlongTrend[IDX.now, ],
#                              ncol = ncol(iXlongTrend))
#       iYTheta.T.now <- iYTheta.T[IDX.now]
#       iXTlocdep.now <- matrix(iXTlocdep[IDX.now, ],
#                               ncol = ncol(iXTlocdep))
#       iXTglobdep.now <- matrix(iXTglobdep[IDX.now, ],
#                                ncol = ncol(iXTglobdep))
#
#       iLinPred.long.now <- iLinPred.long[IDX.now]
#       iProb.T.now <- iProb.T[IDX.now]
#
#       # zeta
#       iGradLogLik[1 : p.long.Trend, k.now] <- VarFac * apply(iXlongTrend.now * (1 - iY.T.now) * (iY.long.now - iLinPred.long.now),
#                                                              MARGIN = 2,
#                                                              sum)
#
#       # eta
#       iGradLogLik[(p.long.Trend + 1) : (p.long.Trend + p.long.AR), k.now] <- VarFac * apply(iXlongAR.now * (1 - iY.T.now) * (iY.long.now - iLinPred.long.now),
#                                                                                             MARGIN = 2,
#                                                                                             sum)
#
#       # beta
#       iGradLogLik[(p.long.Trend + p.long.AR + 1) : (p.long.Trend + p.long.AR + p.long), k.now] <- VarFac * apply(matrix(Xlong[ID.uqe[i], ],
#                                                                                                                         nrow = p.long) %*% matrix((1 - iY.T.now) * (iY.long.now - iLinPred.long.now),
#                                                                                                                                                   nrow = 1),
#                                                                                                                  MARGIN = 1,
#                                                                                                                  sum)
#
#       # sigma.sq
#       iGradLogLik[(p.long.Trend + p.long.AR + p.long + 1) : (n.par.long), k.now] <- - VarFac + 1 / (2 * sigma.sq^2) * sum((1 - iY.T.now) * (iY.long.now - iLinPred.long.now)^2)
#
#       # lambda
#       iGradLogLik[(n.par.long + 1), k.now] <- sum(1  * (iY.T.now - iProb.T.now))
#
#       # theta
#       iGradLogLik[(n.par.long + 2) : ((n.par.long + 1) + p.glob.dep), k.now] <- iXTglobdep.now * (iY.T.now - iProb.T.now)
#
#       # varphi
#       iGradLogLik[(n.par.long + p.glob.dep + 2) : (n.par.long + p.glob.dep + 1 + p.loc.dep), k.now] <- iXTlocdep.now * (iY.T.now - iProb.T.now)
#
#       # xi
#       iGradLogLik[(n.par.long + 1 + p.glob.dep + p.loc.dep + 1) : (n.par.long + n.par.T), k.now] <-  apply(matrix(XT[ID.uqe[i], ],
#                                                                                                                       nrow = p.T) %*% matrix((iY.T.now - iProb.T.now),
#                                                                                                                                              nrow = 1),
#                                                                                                                MARGIN = 1,
#                                                                                                                sum)
#
#
#     }
#     # summarize to gradient
#     GradLogLik[i, ] <- apply(iGradLogLik, MARGIN = 1, sum)
#   }
#
#   SumGradLogLik <- apply(GradLogLik, MARGIN = 2, sum)
#
#   # sum over all subjects to obtain overall gradient
#
#   return(SumGradLogLik)
# }



