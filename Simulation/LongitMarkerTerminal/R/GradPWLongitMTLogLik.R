# Gradient of log likelihood for parametric longitudinal discrete time model for longitudinal
# marker and terminal event;
# partition-specific parameters, piece-wise marker trend and baseline mortality


#' @title The gradient log likelihood for a longitudinal discrete time model for a longitudinal marker and a terminal event data with time-fixed covariates
#' and using unstructured partition-specific autoregressive function.
#' @param par A vector of parameter estimates, in the ordering marker parameters, marker variance and terminal event parameters.
#' @param ID A vector indicating to which patient a particular row belongs.
#' @param TM A vector indicating in which partition interval a particular row belongs (integer counting of intervals).
#' @param Xlong A matrix with baseline covariates (column-wise) for marker model. Each row denotes a particular subject.
#' @param XlongAR A matrix with autoregressive components (column-wise) for marker model; each row denotes a measurement for a particular subject at a particular time.
#' @param XT A matrix with baseline covariates (column-wise) for terminal event model. Each row denotes a particular subject.
#' @param XTglobdep A matrix with marker history for global dependence between marker and terminal event; each row denotes a measurement of a particular subject at a particular time.
#' @param XTlocdep A matrix with marker history for local dependence between marker and terminal event; each row denotes a measurement of a particular subject at a particular time.
#'@param Ylong A vector of marker values; each row denotes a particular subject at a particular time.
#'@param YT A vector of terminal event indicators; each row denotes a particular subject at a particular time.
#'@param riskT A vector of at-risk indicator; each row denotes a particular subject at a particular time.
#' @return The function returns an object of class \code{LongitMT} including estimates and confidence intervals for
#' the time-varying functions and coefficients.
#' @note Used in \code{LongitMTpw.R}.
#' @author Stephanie Armbruster
#' @export

GradPWLongitMTLogLik <- function(par,

                                    ID,
                                    TM,

                                    # input
                                    Xlong,
                                    XlongAR,

                                    XT,
                                    XTglobdep,
                                    XTlocdep,
                                    # output
                                    Ylong,
                                    YT,
                                    riskT

                                    # sigma.sq = 0.25^2
                                    ) {

  # times and ID
  times <- unique(TM)
  K <- length(times) - 1

  ID.uqe <- unique(ID)
  n.sample <- length(ID.uqe)

  # parameter
  p.long <- ncol(Xlong)
  p.long.AR <- ncol(XlongAR)
  n.par.long <- p.long + K * p.long.AR + K + 1

  p.T <- ncol(XT)
  p.glob.dep <- ncol(XTglobdep)
  p.loc.dep <- ncol(XTlocdep)
  n.par.T <- p.T + p.glob.dep + K * p.loc.dep + K


  # longitudinal trend
  zeta <- par[1 : K]
  # AR component
  eta <- matrix(par[(K + 1) : (K + K * p.long.AR)],
                nrow = K, ncol = p.long.AR,
                byrow = T)
  # covariate coefficients
  beta <- par[(K + K * p.long.AR + 1) : (K + K * p.long.AR + p.long)]

  # variance
  sigma.sq <- par[(K + K * p.long.AR + p.long + 1) : (n.par.long)]

  # baseline event rate
  lambda <- par[(n.par.long + 1) : (n.par.long + K)]
  # global dependence
  theta <- par[(n.par.long + K + 1) : ((n.par.long + K) + p.glob.dep)]
  # local dependence
  varphi <- matrix(par[(n.par.long + K + p.glob.dep + 1) : (n.par.long + K + p.glob.dep + K * p.loc.dep)],
                   nrow = K, ncol = p.loc.dep,
                   byrow = T)
  # covariate coefficients
  xi <- par[(n.par.long + K + p.glob.dep + K * p.loc.dep + 1) : (n.par.long + n.par.T)]

  # baseline covariates
  XBeta.Long <- Xlong %*% beta
  XBeta.T <- XT %*% xi

  # history for global dependence
  YTheta.T <- XTglobdep %*% theta

  # save elements
  GradLogLik <- array(0, dim = c(n.sample, n.par.long + n.par.T))

  for (i in seq_len(n.sample)) {
    # indices for pooled data
    IDX <- ID == ID.uqe[i]
    # -1 without baseline; -1 to obtain interval index
    iTM <- TM[IDX][-1] - 1

    # person-time contribution to partition intervals
    ik.l <- min(iTM)
    ik.r <- max(iTM)

    # condition on baseline values
    iY.long <- Ylong[IDX][-1]
    iY.T <- YT[IDX][-1]
    iY.T <- riskT[IDX][-1]

    # condition on baseline values
    iXlongAR <- matrix(XlongAR[IDX, ][-1, ], ncol = ncol(XlongAR))
    iYTheta.T <- matrix(YTheta.T[IDX][-1], ncol = ncol(YTheta.T))
    iXTlocdep <- matrix(XTlocdep[IDX, ][-1, ], ncol = ncol(XTlocdep))
    iXTglobdep <- matrix(XTglobdep[IDX, ][-1, ], ncol = ncol(XTglobdep))


    iLinPred.long <- zeta[iTM] +
      XBeta.Long[ID.uqe[i]] +
      apply(as.matrix(iXlongAR * eta[iTM, ], drop = F),
            MARGIN = 1,
            FUN = function(x) sum(x))

    iLinPred.T <- lambda[iTM] +
      XBeta.T[ID.uqe[i]] +
      iYTheta.T +
      apply(as.matrix(iXTlocdep * varphi[iTM, ], drop = F),
            MARGIN = 1,
            FUN = function(x) sum(x))

    iProb.T <- exp(iLinPred.T) / (1 + exp(iLinPred.T))


    iGradLogLik <- array(0, dim = c(n.par.long + n.par.T, length(unique(iTM))))

    VarFac <- 1 / (sigma.sq)

    for (k.now in unique(iTM)) {

      IDX.now <- iTM == k.now

      iY.long.now <- iY.long[IDX.now]
      iY.T.now <- iY.T[IDX.now]

      iXlongAR.now <- matrix(iXlongAR[IDX.now, ],
                             ncol = ncol(iXlongAR))
      iYTheta.T.now <- iYTheta.T[IDX.now]
      iXTlocdep.now <- matrix(iXTlocdep[IDX.now, ],
                              ncol = ncol(iXTlocdep))
      iXTglobdep.now <- matrix(iXTglobdep[IDX.now, ],
                               ncol = ncol(iXTglobdep))

      iLinPred.long.now <- iLinPred.long[IDX.now]
      iProb.T.now <- iProb.T[IDX.now]


      # zeta_k
      iGradLogLik[k.now, k.now] <- VarFac * sum(1 * (1 - iY.T.now) * (iY.long.now - iLinPred.long.now))

      # eta_k
      iGradLogLik[(K + 2 * k.now - 1) : (K + 2 * k.now), k.now] <- VarFac * apply(iXlongAR.now * (1 - iY.T.now) * (iY.long.now - iLinPred.long.now),
                                                                  MARGIN = 2,
                                                                  sum)

      # beta
      iGradLogLik[(K + K * p.long.AR + 1) : (K + K * p.long.AR + p.long), k.now] <- VarFac * apply(matrix(Xlong[ID.uqe[i], ],
                                                                                           nrow = p.long) %*% matrix((1 - iY.T.now) * (iY.long.now - iLinPred.long.now),
                                                                                                                     nrow = 1),
                                                                                          MARGIN = 1,
                                                                                          sum)

      # sigma.sq
      iGradLogLik[(K + K * p.long.AR + p.long + 1) : (n.par.long), k.now] <- - VarFac / 2 + 1 / (2 * sigma.sq^2) * sum((1 - iY.T.now) * (iY.long.now - iLinPred.long.now)^2)

      # lambda_k
      iGradLogLik[n.par.long + k.now, k.now] <- sum(1  * (iY.T.now - iProb.T.now))

      # theta
      iGradLogLik[(n.par.long + K + 1) : ((n.par.long + K) + p.glob.dep), k.now] <- iXTglobdep.now * (iY.T.now - iProb.T.now)

      # varphi_k
      iGradLogLik[(n.par.long + K + p.glob.dep + 2 * k.now - 1) : (n.par.long + K + p.glob.dep + 2 * k.now ), k.now] <- iXTlocdep.now * (iY.T.now - iProb.T.now)

      # xi
      iGradLogLik[(n.par.long + K + p.glob.dep + K * p.loc.dep + 1) : (n.par.long + n.par.T), k.now] <-  apply(matrix(XT[ID.uqe[i], ],
                                                                                                                      nrow = p.T) %*% matrix((iY.T.now - iProb.T.now),
                                                                                                                                             nrow = 1),
                                                                                                               MARGIN = 1,
                                                                                                               sum)


    }
    # summarize to gradient
    GradLogLik[i, ] <- apply(iGradLogLik, MARGIN = 1, sum)
  }

  SumGradLogLik <- apply(GradLogLik, MARGIN = 2, sum)

  # sum over all subjects to obtain overall gradient

  return(SumGradLogLik)
}



