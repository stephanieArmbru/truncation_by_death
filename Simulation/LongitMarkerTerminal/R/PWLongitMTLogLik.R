# Log likelihood for parametric longitudinal discrete time model for longitudinal
# marker and terminal event;
# partition-specific parameters, piece-wise marker trend and baseline mortality


#' @title The log likelihood for a longitudinal discrete time model for a longitudinal marker and a terminal event data with time-fixed covariates
#' and using unstructured partition-specific autoregressive function.
#' @param par A vector of pareter estimates, in the ordering marker pareters, marker variance and terminal event pareters.
#' @param ID A vector indicating to which patient a particular row belongs.
#' @param TM A vector indicating in which partition interval a particular row belongs (integer counting of intervals).
#' @param Xlong A matrix with baseline covariates (column-wise) for marker model. Each row denotes a particular subject.
#' @param XlongAR A matrix with autoregressive components (column-wise) for marker model; each row denotes a measurement for a particular subject at a particular time.
#' @param XT A matrix with baseline covariates (column-wise) for terminal event model. Each row denotes a particular subject.
#' @param XTglobdep A matrix with marker history for global dependence between marker and terminal event; each row denotes a measurement of a particular subject at a particular time.
#' @param XTlocdep A matrix with marker history for local dependence between marker and terminal event; each row denotes a measurement of a particular subject at a particular time.
#' @param Ylong A vector of marker values; each row denotes a particular subject at a particular time.
#' @param YT A vector of terminal event indicators; each row denotes a particular subject at a particular time.
#' @param riskT A vector of at-risk indicator; each row denotes a particular subject at a particular time.
#' @return The function returns an object of class \code{LongitMT} including estimates and confidence intervals for
#' the time-varying functions and coefficients.
#' @note Used in \code{LongitMTpw.R}.
#' @author Stephanie Armbruster
#' @export

PWLongitMTLogLik <- function(par,

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

  # pareter
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
                byrow = F)
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
                   byrow = F)
  # covariate coefficients
  xi <- par[(n.par.long + K + p.glob.dep + K * p.loc.dep + 1) : (n.par.long + n.par.T)]

  # baseline covariates
  XBeta.long <- Xlong %*% beta
  XBeta.T <- XT %*% xi

  # history for global dependence
  YTheta.T <- XTglobdep %*% theta

  # save elements
  LogLik <- array(0, dim = n.sample)

  for (i in seq_len(n.sample)) {
    # indices for pooled data
    IDX <- ID == ID.uqe[i]
    # -1 without baseline; -1 to obtain interval index
    iTM <- TM[IDX][-1] - 1

    # condition on baseline values
    iY.long <- Ylong[IDX][-1]
    iY.T <- YT[IDX][-1]

    # condition on baseline values
    iXlongAR <- XlongAR[IDX, ][-1, ]
    iYTheta.T <- YTheta.T[IDX][-1]
    iXTlocdep <- XTlocdep[IDX, ][-1, ]


    iLinPred.long <- zeta[iTM] +
      XBeta.long[ID.uqe[i]] +
      apply(as.matrix(iXlongAR * eta[iTM, ], drop = F),
            MARGIN = 1,
            FUN = function(x) sum(x))

    iLinPred.T <- lambda[iTM] +
      XBeta.T[ID.uqe[i]] +
      iYTheta.T +
      apply(as.matrix(iXTlocdep * varphi[iTM, ], drop = F),
            MARGIN = 1,
            FUN = function(x) sum(x))


    iLogLik <- (1 - iY.T) * (-1/2) * log(2 * pi * sigma.sq) +
      (-1/(2 * sigma.sq)) * (iY.long - iLinPred.long)^2 +
      (1 - iY.T) * (1 - exp(iLinPred.T) / (1 + exp(iLinPred.T))) +
      iY.T * exp(iLinPred.T) / (1 + exp(iLinPred.T))

    # save summed contribution
    LogLik[i] <- sum(iLogLik)
  }

  return(-sum(LogLik))
}



