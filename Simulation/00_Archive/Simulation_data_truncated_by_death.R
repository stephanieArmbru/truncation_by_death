#### SIMULATION OF SEMICOMPETING SURIVAL DATA AND LONGITUDINAL COVARIATE ####
set.seed(2024)


# according to Harrison's data generating mechanism
# based on B-spline hazard functions
# inverse sampling method

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(checkmate)
library(SemiCompRisksFreq)




# generate random walk
RW <- function(n, x0, mu, variance) {
  # cumulative random error
  z <- cumsum(rnorm(n = n,
                    mean = 0,
                    sd = sqrt(variance)))
  # discrete time
  t <- seq_along(n)

  # sum
  x <- x0 + t * mu + z
  return(x)
}



# Parameters --------------------------------------------------------------
# we are starting with a simple case: one non-terminal, one terminal event
# one longitudinal covariate
# for now no clustering

# assume Markov property for longitudinal covariate

# generate example values
n <- 100
K <- 100
Lmat <- matrix(data = RW(n = n * K, x0 = 10, mu = 0.5, variance = 0.1),
               nrow = n,
               ncol = K)
x <- matrix(rnorm(3 * n), ncol = 3)
x1 <- x2 <- x3 <- x

tau1.true <- 0.5
tau2.true <- 0.1
tau3.true <- 0.6

beta1.true <- c(0.25,0.6,-0.3)
beta2.true <- c(0.4,0.75,-0.4)
beta3.true <- c(0.7,0.9,-0.9)

alpha1.true <- exp(-0.6)
alpha2.true <- exp(-0.3)
alpha3.true <- exp(-0.4)

kappa1.true <- exp(-1.2)
kappa2.true <- exp(-2.9)
kappa3.true <- exp(-3.2)

hazard <- "Weibull"
model <- "semi-Markov"

frailty_type <- "gamma"
theta.true <- 0



h3tv_degree <- 0
knots_temp <- c(3, 8, 15.5)
h3tv_knots <- c(0,knots_temp,Inf)
beta3tv.true <- c(1.5,0,-1.5)

knots_list <- list(knots_temp)


LT_interval <- c(0,0)
cens <- c(60,60)



simSCL <- function(n,

                   # longitudinal covariate
                   Lmat,
                   tau1.true, tau2.true, tau3.true,

                   # proportional hazards components
                   x1, x2, x3,
                   beta1.true, beta2.true, beta3.true,

                   # baseline (shape / scale)
                   alpha1.true, alpha2.true, alpha3.true,
                   kappa1.true, kappa2.true, kappa3.true,

                   phi1.true, phi2.true, phi3.true,

                   # frailty
                   theta.true,
                   frailty_type = "gamma",


                   # clustering id
                   id = NULL,
                   # true covariance for cluster effect
                   SigmaV.true = NULL,

                   hazard = "weibull",
                   knots_list,

                   model = "semi-markov",

                   beta2frail.true = 1,
                   beta3frail.true = 1,
                   beta3tv.true = NULL,

                   h3tv_degree = 3,
                   h3tv_knots = NULL,

                   LT_interval = c(0, 0),
                   cens = c(0, 0)
                   ) {
  n <- nrow(x1)

  # use checkmate to perform tests on arguments
  if (!(checkMatrix(Lmat) & nrow(Lmat) == n)) {
    stop("Longitudinal covariate must be a matrix with one row for each subject.")
  }
  if (!all(checkMatrix(x1),
         checkMatrix(x2),
         checkMatrix(x3),

         checkNumeric(x1),
         checkNumeric(x2),
         checkNumeric(x3),

         nrow(x1) == nrow(x2),
         nrow(x2) == nrow(x3))) {
    stop("Covariate input X1, X2, X3 must be a numeric matrix with one row for each subject.")
  }

  if (!all(checkVector(beta1.true),
           checkVector(beta2.true),
           checkVector(beta3.true),

           checkNumeric(beta1.true),
           checkNumeric(beta2.true),
           checkNumeric(beta3.true),

           length(beta1.true) == ncol(x1),
           length(beta2.true) == ncol(x2),
           length(beta3.true) == ncol(x3)

  )) {
    stop("Beta input must be a vector and have same dimensions as corresponding design matrix.")
  }

  if (!match.arg(frailty_type, c("gamma", "lognorm", "none"))) {
    stop("Frailty term must be either Gamma or lognormal distributed.
         Other options are not implemented.")
  }

  if (frailty_type == "none") {
    if (!is.null(theta.true)) {
      stop("If no frailty should be included, the true theta should be set to NULL.")
    }
  } else if (frailty_type == "gamma") {
    if (!(checkInteger(theta.true) & theta.true != 0)) {
      stop("For gamma frailty, the theta value must be non-zero integer.")
    }
  } else if (frailty_type == "lognorm") {
    if(theta.true <= 0) {
      stop("For lognormal frailty, the theta value must be strictly positive.")
    }
  }

  if (!match.arg(hazard, c("weibull", "piecewise"))) {
    stop("The baseline hazard function can be either piecewise defined (according to knots argument) or a Weibull distribution.")
  }


  if (any(LT_interval < 0, cens < 0)) {
    stop("Left / Right censoring boundaries must be non-negative.")
  }
  if (LT_interval[1] > LT_interval[2]) {
    stop("Interval boundaries for left censoring incorrect.
         Left boundary cannot be larger than right boundary.")
  }
  if (cens[1] > cens[2]) {
    stop("Interval boundaries for right censoring incorrect.
         Left boundary cannot be larger than right boundary.")
  }

  if (hazard == "weibull") {
    if (!is.null(knots_list)) {
      stop("If the hazard is defined to be Weibull, the argument knots_list must be NULL.")
    }
  }

  if (hazard == "piecewise") {
    if (!(checkList(knots_list))) {
      stop("Knots_list must be a list.")
    }
    if (!any(knots_list %>% length() == 3, knots_list %>% length() == 1)) {
      stop("Knots must be either a list with one element containing the knots lambdas,
           or a list with three elements corresponding to hazard-specific knot vectors.")
    }
  }





  if (!is.null(id) & is.null(SigmaV.true)) {
    stop("SigmaV.true must be given to simulate correlated data.")
  }

  # establish dimensions
  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]

  # left censoring included?
  anyLT <- max(LT_interval) > 0

  # sample frailty
  if (frailty_type == "gamma") {
      gamma.true <- stats::rgamma(n, 1 / theta.true, 1 / theta.true)
    } else if (frailty_type == "lognorm") {
      gamma.true <- stats::rlnorm(n,
                                  meanlog = 0,
                                  sdlog = sqrt(theta.true))
    } else if(frailty_type == "none") {
    gamma.true <- rep(1, n)
    }

  # compute predictor
  LP1	<- if (p1 > 0) as.vector(x1 %*% beta1.true) else 0
  LP2	<- if (p2 > 0) as.vector(x2 %*% beta2.true) else 0
  LP3	<- if (p3 > 0) as.vector(x3 %*% beta3.true) else numeric(n)
  #made a vector bc it is subset automatically by yesR below

  # add longitudinal component
  long1 <- tau1.true * Lmat[, -1]
  long2 <- tau2.true * Lmat[, -1]
  long3 <- tau3.true * Lmat[, -1]


  #incorporate clustering random effects
  if (!is.null(id)) {
    J <- length(unique(id))
    nj <- as.vector(table(id))
    # sample multivariate normal cluster frailty, mean zero 0
    Vmat <- MASS::mvrnorm(J,
                          rep(0, 3),
                          SigmaV.true)

    Vmat <- cbind(V1 = rep(Vmat[,1], nj),
                  V2 = rep(Vmat[,2], nj),
                  V3 = rep(Vmat[,3], nj))
    # add to predictor
    LP1 <- LP1 + as.vector(Vmat[,1])
    LP2 <- LP2 + as.vector(Vmat[,2])
    LP3 <- LP3 + as.vector(Vmat[,3])
  } else {Vmat <- NULL}

  # for piecewise hazard
  if (hazard == "piecewise"){
    #make it so that index can be made rowwise whether there is 1 col or 3
    #function already accounts for whether or not knots vectors begin with 0
    if (length(knots_list) == 3){
        knots01 <- knots_list[[1]]
        knots02 <- knots_list[[2]]
        knots03 <- knots_list[[3]]
      } else if (length(knots_list) == 1){
        knots01 <- knots02 <- knots03 <- knots_list[[1]]
      }
  }

  # generate censoring time
  if (cens[2] == 0){
    Cen <- rep(Inf,n)
  } else {
    #fix what to do about intersection of left truncation and censoring interval!

    # if(cens[1] < LT_interval[2]){
    #   warning(paste0("Censoring distribution cannot overlap truncation distribution.",
    #                  "Setting minimum censoring time to ", LT_interval[2],"."))
    #   cens[1] <- LT_interval[2]
    # }
    Cen <- stats::runif(n,
                        min = cens[1],
                        max = cens[2])
  }

  # generate left truncation
  if (anyLT) {
    yL <- stats::runif(n,
                       min = LT_interval[1],
                       max = LT_interval[2])
    # for weibull baseline hazard
    if(hazard == "weibull"){
      # left truncation means conditioning on being above left boundary
      # sample time-to-event based on inverse probability sampling
      R_bound <- stats::pweibull(yL,
                                 lower.tail = FALSE,
                                 shape = alpha1.true,
                                 scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true)) / alpha1.true))

      D_bound <- stats::pweibull(yL,
                                 lower.tail = FALSE,
                                 shape = alpha2.true,
                                 scale = exp(-(log(kappa2.true) + LP2 + beta2frail.true * log(gamma.true))/alpha2.true))

      R_prob <- stats::runif(n,
                             min = 0,
                             max = R_bound)
      D_prob <- stats::runif(n,
                             min = 0,
                             max = D_bound)


      R <- stats::qweibull(R_prob,
                           lower.tail = FALSE,
                           shape = alpha1.true,
                           scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true)) / alpha1.true))
      D <- stats::qweibull(D_prob,
                           lower.tail = FALSE,
                           shape = alpha2.true,
                           scale = exp(-(log(kappa2.true) + LP2 + beta2frail.true * log(gamma.true)) / alpha2.true))
    } else{
      R_bound <- ppwexp(yL,
                        lower.tail = FALSE,
                        log.p=FALSE,
                        phi=phi1.true,
                        knots_vec=knots01,
                        eta = LP1 + log(gamma.true))

      D_bound <- ppwexp(yL,
                        lower.tail = FALSE,
                        log.p=FALSE,
                        phi=phi2.true,
                        knots_vec=knots02,
                        eta = LP2 + beta2frail.true * log(gamma.true))

      R_prob <- stats::runif(n,min=0, max=R_bound)
      D_prob <- stats::runif(n,min=0, max=D_bound)

      R <- qpwexp(R_prob, lower.tail = FALSE, phi=phi1.true, knots_vec=knots01,
                  eta = LP1 + log(gamma.true))
      D <- qpwexp(D_prob, lower.tail = FALSE, phi=phi2.true, knots_vec=knots02,
                  eta = LP2 + beta2frail.true * log(gamma.true))
    }
  } else {
    yL <- numeric(n)
    if(hazard == "weibull"){
      R <- stats::rweibull(n,
                           shape = alpha1.true,
                           scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true)) / alpha1.true))
      D <- stats::rweibull(n, shape = alpha2.true,
                           scale = exp(-(log(kappa2.true) + LP2 + beta2frail.true * log(gamma.true)) / alpha2.true))
    } else{
      #still generate data by inverse probability
      R <- qpwexp(stats::runif(n),
                  lower.tail = FALSE,
                  phi=phi1.true,
                  knots_vec=knots01,
                  eta = LP1 + log(gamma.true))
      D <- qpwexp(stats::runif(n),
                  lower.tail = FALSE,
                  phi=phi2.true,
                  knots_vec=knots02,
                  eta = LP2 + beta2frail.true * log(gamma.true))
    }
  }

  # generate indicator for event
  yesR <- R < D

  #now, incorporate a possibly time-varying component into LP3 (semi-markov only)
  #if we set total number of parameters to 0, then we have no time-varying component.
  if( !(tolower(model) %in% c("markov","m")) && !is.null(beta3tv.true)){
    p3tv <- length(beta3tv.true)
    if(h3tv_degree == "linear"){ #linear
      stopifnot(p3tv==1)
      x3tv <- as.matrix(pmin(R,D,Cen))
      h3tv_knots <- c(0,Inf)
    } else if(h3tv_degree == "log1p") {
      stopifnot(p3tv==1)
      x3tv <- as.matrix(log1p(pmin(R,D,Cen)))
      h3tv_knots <- c(0,Inf)
    } else if(h3tv_degree == "cs"){ #cubic spline model
      #in cubic spline model, boundary knots are set directly at min/max endpoints,
      #so no need to fix at 0
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+1)
      if(is.null(h3tv_knots)){
        h3tv_knots <- stats::quantile(pmin(R,Cen)[yesR==1 & R<Cen], h3_quantile_seq)
      }
      x3tv <- splines::ns(x = pmin(R,D,Cen), knots = h3tv_knots[-c(1,length(h3tv_knots))],
                          Boundary.knots = h3tv_knots[c(1,length(h3tv_knots))],
                          intercept = FALSE)
    } else { #if we don't use restricted cubic, then we are using a regular b-spline with specified degree
      #this also includes piecewise constant if degree is 0
      h3tv_degree <- as.numeric(h3tv_degree)
      stopifnot(p3tv>=h3tv_degree)
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+2-h3tv_degree)[-c(1,p3tv+2-h3tv_degree)]
      #fixing piecewise endpoint at maximum is ok, because splines2 prediction will extrapolate beyond it
      if(is.null(h3tv_knots)){
        h3tv_knots <- c(0,stats::quantile(pmin(R,Cen)[yesR==1 & R<Cen],
                                          h3_quantile_seq),max(pmin(R,D,Cen)))
      }
      x3tv <- splines2::bSpline(x = pmin(R,D,Cen), intercept = FALSE, degree = h3tv_degree,
                                knots = h3tv_knots[-c(1,length(h3tv_knots))],
                                Boundary.knots = h3tv_knots[c(1,length(h3tv_knots))])
    }
    colnames(x3tv) <- paste0("h3tv",1:p3tv)
    LP3 <- LP3 + x3tv %*% beta3tv.true
  } else{
    p3tv <- 0
  }

  #check the lower.bound stuff to make sure that I have the right probs and not 1-probs
  if(tolower(model) %in% c("markov","m")){
    if(tolower(hazard) %in% c("weibull","wb")){
      M_bound <- stats::pweibull(R[yesR],lower.tail = FALSE, shape = alpha3.true,
                                 scale = exp(-(log(kappa3.true) + LP3[yesR] +
                                                 beta3frail.true * log(gamma.true[yesR]))/alpha3.true),)
      M_prob <- stats::runif(sum(yesR),min=0, max=M_bound)
      D[yesR] <- stats::qweibull(p = M_prob, lower.tail = FALSE, shape = alpha3.true,
                                 scale = exp(-(log(kappa3.true) + LP3[yesR] +
                                                 beta3frail.true * log(gamma.true[yesR]))/alpha3.true))
    } else{
      M_bound <- ppwexp(R[yesR], lower.tail = FALSE, log.p=FALSE,
                        phi=phi3.true, knots_vec=knots03,
                        eta = LP3[yesR] + beta3frail.true * log(gamma.true[yesR]))
      M_prob <- stats::runif(sum(yesR),min=0, max=M_bound)
      D[yesR] <- qpwexp(M_prob, lower.tail = FALSE, phi=phi3.true,
                        knots_vec=knots03, eta = LP3[yesR] + beta3frail.true * log(gamma.true[yesR]))
    }
  } else{
    if(tolower(hazard) %in% c("weibull","wb")){
      D[yesR] <- R[yesR] + stats::rweibull(sum(yesR), shape = alpha3.true,
                                           scale = exp(-(log(kappa3.true) + LP3[yesR] +
                                                           beta3frail.true * log(gamma.true[yesR]))/alpha3.true))
    } else {
      D[yesR] <- R[yesR] + qpwexp(stats::runif(sum(yesR)), lower.tail = FALSE,
                                  phi=phi3.true, knots_vec=knots03,
                                  eta =  LP3[yesR] + beta3frail.true * log(gamma.true[yesR]))
    }
  }
  delta1 <- rep(NA, n)
  delta2 <- rep(NA, n)
  y1 <- R
  y2 <- D

  #cases where terminal occurs before non-terminal and censoring
  ind01 <- which(D < R & D < Cen)
  y1[ind01] <- D[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1

  #cases where nonterminal occurs, then censoring before terminal
  ind10 <- which(R < D & R < Cen & D >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0

  #cases where censoring occurs first
  ind00 <- which(R >= Cen & D >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0

  #cases where nonterminal occurs, then terminal, then censoring
  ind11 <- which(R < Cen & D < Cen & R < D)
  delta1[ind11] <- 1
  delta2[ind11] <- 1

  #this is a gut-check that the values I build the basis for h3tv on above
  #match the values of y1 for all observations that matter (e.g., those with delta1==1)
  stopifnot(all((y1==pmin(R,Cen))[delta1==1]))

  ret <- data.frame(cbind(y1, delta1, y2, delta2, yL, gamma.true,
                          if(p3tv > 0) x3tv,
                          id, Vmat))
  if(!is.null(beta3tv.true)){
    attr(ret,which = "p3tv") <- p3tv
    attr(ret,which = "h3tv_degree") <- h3tv_degree
    attr(ret,which = "h3tv_knots") <- h3tv_knots
  }

  # if(!(tolower(hazard) %in% c("weibull","wb"))){
  #   attr(ret,which = "knots01") <- knots01
  #   attr(ret,which = "knots02") <- knots02
  #   attr(ret,which = "knots03") <- knots03
  # }

  return(ret)
}


# Exponential simulation -------------------------------------------------
n = 100
rate=c(1, 1.5, -1)
intervals=c(1, 2, 3)
cumulative=FALSE


if(is.null(intervals)){
  if (cumulative){return(cumsum(stats::rexp(n,rate[1])))}else
    return(stats::rexp(n,rate[1]))}

k <- length(rate)
if (k==1){
  if(cumulative){return(cumsum(stats::rexp(n,rate)))}else
    return(stats::rexp(n,rate))
}
if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
tx <- 0
j <- 1
times <- array(0,n)
timex <- cumsum(intervals)
indx <- array(TRUE,n)
for(i in 1:k){
  nindx <- sum(indx)
  if (nindx==0) break
  increment <- stats::rexp(nindx,rate[i])
  # here I do not understand the cumsum() usage
  if (cumulative) times[indx] <- tx + cumsum(increment)
  else times[indx] <- tx + increment
  if (i<k){
    tx <- timex[i]
    indx <- (times > timex[i])
  }
}
# sampled times
times

