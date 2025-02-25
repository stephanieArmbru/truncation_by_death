#' @title The function to fit a longitudinal discrete time model for a longitudinal marker and a terminal event data with time-fixed covariates
#' and using parametric autoregressive function.
#' @description The function implements the proposed methodology XX (my future paper) for time-fixed covariates under possible
#' right censoring and left truncation. Data should be first converted to longitudinal bivariate representation by applying
#' \code{\link{TimesToLongitMT}} function. This function fits the model based on maximum likelihood method.
#' @param longit.data A list or data frame with entries named \code{ID} (ID column), \code{TM} (Time column), \code{risk.T},  \code{YT}, \code{Ylong}, as well as the marker transformations as detailed in
#' formula.long.AR, formula.glob.dep, formula.loc.dep. and the baseline covariates stacked over partition interval and subject (in that order).
#' See details below.
#' @param times A vector of increasing times (for example, the interval partition points \eqn{\tau_1}... \eqn{\tau_K}).
#' @param formula.long A formula of the form \code{ ~ x1 + x2 } where \code{x1} and \code{x2} are covariates to be used for the longitudinal marker sub-model.
#' @param formula.long.AR A formula of the form \code{~ y1 + y2} where \code{y1} and \code{y2} are the autoregressive components
#' (e.g. lag1 marker value; marker slop) for the longitudinal marker sub-model.
#' @param formula.long.trend A formula of the form \code{~ f(time)} where \code{f(time)} describes how the marker trend depends on time (true time variable).
#' @param formula.T A formula of the form \code{ ~ x1 + x3} where \code{x1} and \code{x3} are covariates to be used for
#' for the terminal event probability sub-model.
#' @param formula.glob.dep A formula of the form \code{ ~ y1 + y3} where \code{y1} and \code{y3} are some transformations of the marker history
#' to be used to model the time-invariant global dependence between the marker and the terminal event in the terminal event probability sub-model.
#' @param formula.loc.dep A formula of the form \code{ ~ y1 + y4} where \code{y1} and \code{y4} are some transformations of the marker history
#' to be used to model the time-variant local dependence between the marker and the terminal event in the terminal event probability sub-model.
#' @param init Initial values for the parameters.
#' @param maxit.optim For internal use of \code{optim}. Default is 50000.
#' @details Each of the matrices in \code{longit.data},  \code{risk.T}, \code{YT}, \code{Ylong} have a row for each
#' unit and a column for each interval. Then, \code{risk.T} indicate whether the unit is at risk of the terminal event in each interval.
#' The matrix \code{YT} indicates whether the terminal event was observed by the end of each interval.
#' The function \code{\link{TimesToLongitMT}} can be used to obtain this representation of marker and terminal event data.
#' @return The function returns an object of class \code{LongitMT} including estimates and confidence intervals for
#' the time-varying functions and coefficients.
#' @note For time-varying covariates use \code{\link{LongitMTtimeDep}} or \code{\link{LongitMTparamTimeDep}}.
#' @examples
#' \dontrun{}
#' @author Stephanie Armbruster
#' @export
LongitMTparam <- function(longit.data,
                          times = NULL,
                          formula.long,
                          formula.long.AR,
                          formula.long.trend,

                          formula.T,
                          formula.glob.dep,
                          formula.loc.dep,

                          data,

                          init = NULL,
                          maxit.optim = 50000)
{
  ID <- longit.data$ID
  TM <- longit.data$TM

  n.sample <- ID %>% unique() %>% length()

  # Longitudinal model
  if (!is.null(formula.long)) {
    # [-1] to avoid intercept column
    Xlongmat <- as.matrix(model.matrix(formula.long, data = longit.data)[, -1])
    p.long <- ncol(Xlongmat)
  } else {
    Xlongmat <- NULL
    p.long <- 0
  }

  if(!is.null(formula.long.AR)) {
    # remove baseline values
    XlongARmat <- as.matrix(model.matrix(formula.long.AR, data = longit.data)[, -1])
    # [-seq(1, n.sample), ]
    p.long.AR <- ncol(XlongARmat)
  } else {
    XlongARmat <- NULL
    p.long.AR <- 0
  }

  if (!is.null(formula.long.trend)) {
    XlongTrendmat <- as.matrix(model.matrix(formula.long.trend, data = longit.data))
    p.long.Trend <- ncol(XlongTrendmat)
  } else {
    XlongTrendmat <- NULL
    p.long.Trend <- 0
  }

  # Terminal event model
  if (!is.null(formula.T)) {
    XTmat <- as.matrix(model.matrix(formula.T, data = longit.data)[, -1])
    p.T <- ncol(XTmat)
  } else {
    XTmat <- NULL
    p.T <- 0
  }
  if (!is.null(formula.glob.dep)) {
    XTglobdepmat <- as.matrix(model.matrix(formula.glob.dep, data = longit.data)[, -1])
    # [-seq(1, n.sample), ]
    p.T.globdep <- ncol(XTglobdepmat)
  } else {
    XTglobdepmat <- NULL
    p.T.globdep <- 0
  }
  if (!is.null(formula.loc.dep)) {
    XTlocdepmat <- as.matrix(model.matrix(formula.loc.dep, data = longit.data)[, -1])
    # [-seq(1, n.sample), ]
    p.T.locdep <- ncol(XTlocdepmat)
  } else {
    XTlocdepmat <- NULL
    p.T.locdep <- 0
  }


  K <- length(unique(TM)) - 1
  if (is.null(times)) times <- 1:K

  # parameter
  n.param.long <- p.long + p.long.AR + p.long.Trend + 1
  n.param.T <- p.T + p.T.globdep + p.T.locdep + 1

  # + 1 for variance
  n.params <- n.param.long + n.param.T

  # initialize coefficient estimate
  if (is.null(init)) {
    init <- rep(1, n.params)
  }


  # ParamLongitMTLogLik(ID = ID,
  #                     TM = TM,
  #                     Xlong = Xlongmat,
  #                     XlongAR = XlongARmat,
  #                     XlongTrend = XlongTrendmat,
  #
  #                     XT = XTmat,
  #                     XTglobdep = XTglobdepmat,
  #                     XTlocdep = XTlocdepmat,
  #
  #                     Ylong = longit.data$Ylong,
  #                     YT = longit.data$YT,
  #                     # sigma.sq = 0.25^2,
  #                     par = init)
  # GradParamLongitMTLogLik(ID = ID,
  #                         TM = TM,
  #                         Xlong = Xlongmat,
  #                         XlongAR = XlongARmat,
  #                         XlongTrend = XlongTrendmat,
  #
  #                         XT = XTmat,
  #                         XTglobdep = XTglobdepmat,
  #                         XTlocdep = XTlocdepmat,
  #
  #                         Ylong = longit.data$Ylong,
  #                         YT = longit.data$YT,
  #                         sigma.sq = 0.25^2,
  #                         par = init)


  res.opt <- optim(par = init,
                   # observed likelihood
                   fn = ParamLongitMTLogLik,

                   # gradient for observed likelihood
                   # gr = GradParamLongitMTLogLik,

                   hessian = T,
                   control = list(maxit = maxit.optim),

                   method = "L-BFGS-B",

                   # input for function
                   ID = ID,
                   TM = TM,
                   Xlong = Xlongmat,
                   XlongAR = XlongARmat,
                   XlongTrend = XlongTrendmat,

                   XT = XTmat,
                   XTglobdep = XTglobdepmat,
                   XTlocdep = XTlocdepmat,

                   Ylong = longit.data$Ylong,
                   YT = longit.data$YT
                   # sigma.sq = 0.25^2
                   )


  # save results
  fit <- list()
  fit$formula.long <- formula.long
  fit$formula.long.AR <- formula.long.AR
  fit$formula.long.trend <- formula.long.trend

  fit$formula.T <- formula.T
  fit$formula.T.glob.dep <- formula.glob.dep
  fit$formula.T.loc.dep <- formula.loc.dep

  fit$optim.conv <- res.opt$convergence
  fit$est <- res.opt$par
  fit$lik <- -res.opt$value
  fit$hess <- res.opt$hessian

  # covariance matrix for estimates
  # fit$v.hat <- solve(res.opt$hessian)
  # fit$se <- sqrt(diag(fit$v.hat))
  # if (any(is.nan(fit$se)))
  # {
  #   fit$se[is.nan(fit$se)] <- 100
  #   warning("close to infinite SE for some parameters, set SE=100")
  # }
  # find coefficient estimates for partition-specific baseline values,
  # covariate coefficients and SE for each coefficient

  # longitudinal trend
  fit$coef.zeta <- fit$est[1 : p.long.Trend]
  # AR component
  fit$coef.eta <- fit$est[(p.long.Trend + 1):(p.long.Trend + p.long.AR)]
  # covariate coefficients
  fit$coef.beta <- fit$est[(p.long.Trend + p.long.AR + 1) : (p.long.Trend + p.long.AR + p.long)]

  # variance
  fit$coef.sigma.sq <- fit$est[n.param.long]

  # baseline event rate
  fit$coef.lambda <- fit$est[(n.param.long + 1)]
  # global dependence
  fit$coef.theta <- fit$est[(n.param.long + 2) : ((n.param.long + 1) + p.T.globdep)]
  # local dependence
  fit$coef.varphi <- fit$est[(n.param.long + p.T.globdep + 2) : (n.param.long + p.T.globdep + 1 + p.T.locdep)]
  # covariate coefficients
  fit$coef.xi <- fit$est[(n.param.long + 1 + p.T.globdep + p.T.locdep + 1) : (n.param.long + n.param.T)]

  # standard error
  # longitudinal trend
  # fit$se.zeta <- fit$se[1 : p.long.Trend]
  # # AR component
  # fit$se.eta <- matrix(fit$se[(p.long.Trend + 1) : (p.long.Trend + p.long.AR)],
  #                      nrow = K, ncol = p.long.AR,
  #                      byrow = F)
  # # covariate coefficients
  # fit$se.beta <- fit$se[(p.long.Trend + p.long.AR + 1) : (p.long.Trend + p.long.AR + p.long)]
  #
  # # variance
  # fit$se.sigma.sq <- fit$se[n.param.long]
  #
  # # baseline event rate
  # fit$se.lambda <- fit$se[n.param.long + 1]
  # # global dependence
  # fit$se.theta <- fit$se[(n.param.long + 2) : ((n.param.long + 1) + p.T.globdep)]
  # # local dependence
  # fit$se.varphi <- matrix(fit$se[(n.param.long + p.T.globdep + 2) : (n.param.long + p.T.globdep + 1 + p.T.locdep)],
  #                         nrow = K, ncol = p.T.locdep,
  #                         byrow = F)
  # # covariate coefficients
  # fit$se.xi <- fit$se[(n.param.long + 1 + p.T.globdep + p.T.locdep + 1) : (n.param.long + n.param.T)]


  # calculate ci for the baseline prob.T1, prob.T2 and OR.
  # Using the appropriate transformation of the CI for B%*%eta
  sub.vhat.long <- fit$v.hat[1:n.param.long, 1:n.param.long]
  sub.vhat.T <- fit$v.hat[(n.param.long + 1):(n.param.long + n.param.T), (n.param.long + 1):(n.param.long + n.param.T)]

  # add in confidence interval estimation (sandwich estimator)

  # add class definition
  # class(fit) <- "LongitTerminal"
  fit %>% return()
}




