

try <- function(x, a, b) {((b - a) / 2) * x + ((b + a) / 2)}

x <- seq(-1, 1, by = 0.1)

ggplot(data = data.frame(x = x, y = 1),
       aes(x = x, y = y)) +
  geom_point()

x_transform <- try(x = x, a = 2, b = 5)


ggplot(data = data.frame(x = x_transform, y = 1),
       aes(x = x, y = y)) +
  geom_point()


# prediction function

temp_fit <- fit_cov
x_temp <- example_bl

tseq = tseq_pred
para = temp_fit$estimate
x1new = x_temp
x2new = x_temp
x3new = x_temp
frailnew = NULL
frailty = temp_fit$frailty
model = temp_fit$model
nP0 = temp_fit$nP0
nP = temp_fit$nP
p3tv = 0
h3tv_basis_func = h3_b_fun
hazard = temp_fit$hazard
knots_list = temp_fit$knots_list
n_quad = temp_fit$n_quad
quad_method = temp_fit$quad_method
Finv = temp_fit$Finv
alpha = 0.05





marg_flag <- is.null(frailnew)
if(is.null(frailnew)){
  frailnew <- 1
}

#it's crazy but I think I'll specify things related to the effect of t1 on h3 as follows:
#you tell me how many total parameters are used in estimating the h3tv effect
#and you give me the function that defines the basis for that effect.
#Could be linear, piecewise, spline, I don't care. All I know is that it's an effect on the log-hazard scale.
#also, assume that it's the last thing in the vector of parameters

#setting up some preliminary quantities
nP0_tot <- if(frailty) sum(nP0) + 1 else sum(nP0)
nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3]-p3tv) #leave off the h3tv effects!!!

eta_list <- beta_list <- vector(mode="list",3)
beta_list[[1]] <- if(!is.null(x1new)) para[nP_start[1]:nP_end[1]] else numeric(0)
beta_list[[2]] <- if(!is.null(x2new)) para[nP_start[2]:nP_end[2]] else numeric(0)
beta_list[[3]] <- if(!is.null(x3new)) para[nP_start[3]:nP_end[3]] else numeric(0)
eta_list[[1]] <- if(length(beta_list[[1]]) > 0) as.vector(x1new %*% beta_list[[1]] + log(frailnew)) else 0
eta_list[[2]] <- if(length(beta_list[[2]]) > 0) as.vector(x2new %*% beta_list[[2]] + log(frailnew)) else 0
eta_list[[3]] <- if(length(beta_list[[3]]) > 0) as.vector(x3new %*% beta_list[[3]] + log(frailnew)) else 0

h3tv_phi <- if(p3tv>0) utils::tail(para,n=p3tv) else numeric(0)

se_fit_flag <- !is.null(Finv)
if(frailty){
  theta <- exp(para[sum(nP0) + 1])
}

#prefill the list with all potential inputs
value_temp <- list()
value <- list(tseq=tseq,
              p_term_only=NULL, se_p_term_only=NULL,
              p_nonterm=NULL, se_p_nonterm=NULL,
              p_nonterm_only=NULL, se_p_nonterm_only=NULL,
              p_both=NULL, se_p_both=NULL,
              p_neither=NULL, se_p_neither=NULL, se_logp_neither=NULL,
              p_term_only_marg=NULL, se_p_term_only_marg=NULL,
              p_nonterm_marg=NULL, se_p_nonterm_marg=NULL,
              p_nonterm_only_marg=NULL, se_p_nonterm_only_marg=NULL,
              p_both_marg=NULL, se_p_both_marg=NULL,
              p_neither_marg=NULL, se_p_neither_marg=NULL,
              x1new=x1new, x2new=x2new, x3new=x3new, frailnew=frailnew, class=NULL)

quad_weights <- SemiCompRisksFreq:::get_quad_pointsweights(n_quad=n_quad,
                                                           quad_method=quad_method)$weights
quad_points <- SemiCompRisksFreq:::transform_quad_points(n_quad = n_quad,
                                                         quad_method=quad_method, a = 0,b = tseq)
H_list <-  vector(mode="list",3)
for(i in 1:3){
  xtemp <- switch(i, x1new, x2new, x3new)

  basis <- get_basis(y = tseq,knots_vec = knots_list[[i]],hazard = hazard, deriv = FALSE)
  dbasis <- get_basis(y = tseq,knots_vec = knots_list[[i]],hazard = hazard, deriv = TRUE)
  basis_quad <- get_basis(y=quad_points, knots_vec=knots_list[[i]],hazard=hazard,deriv = FALSE)

  #don't bother with h3_tv here because H_list[[3]] is only used by semi-markov setting
  H_list[[i]] <- get_haz(tseq=tseq, eta=eta_list[[i]],
                         hazard=hazard, phi=para[nP0_start[i]:nP0_end[i]],
                         basis=basis, dbasis=dbasis,
                         basis_quad=basis_quad, quad_weights=quad_weights,
                         func_type="H", log_out=FALSE)
}

# browser()

#if tseq=0, just fix these.
#Then, it can be the default for tseq to start with 0
if(tseq[1]==0){
  H_list[[1]][1,] <- 0
  H_list[[2]][1,] <- 0
  H_list[[3]][1,] <- 0
}

#next, compute conditional probabilities
#note, for the moment these are length(tseq) by length(eta) matrices!
value$p_neither <- exp(-H_list[[1]] - H_list[[2]])

#make it so that you've "started integrating" from time tseq[1]
#in practice, this is no different if tseq = 0
p_neither0 <- exp(-H_list[[1]][1,] - H_list[[2]][1,])
value$p_neither <- t(t(value$p_neither) / p_neither0)

#instead of assuming that there's an implicit 0, just directly use the first thing!
#this allows for left truncation by just beginning the integral at an arbitrary time
S1_b <- exp(-H_list[[1]][-1,,drop=FALSE])
S1_a <- exp(-H_list[[1]][-length(tseq),,drop=FALSE])
S2_b <- exp(-H_list[[2]][-1,,drop=FALSE])
S2_a <- exp(-H_list[[2]][-length(tseq),,drop=FALSE])

#and then in turn, we know that wherever the integral is starting from, it's
#gotta start at 0, so just add that back
value$p_term_only <- rbind(0, 0.5 *
                             apply( (S1_a + S1_b) * (S2_a - S2_b), MARGIN = 2, cumsum))
value$p_nonterm <- rbind(0, 0.5 *
                           apply( (S1_a - S1_b) * (S2_a + S2_b), MARGIN = 2, cumsum))

#if there is left truncation, then correct for it by rescaling by the probability
#of no event by the starting timepoint
value$p_nonterm <- t(t(value$p_nonterm) / p_neither0)
value$p_term_only <- t(t(value$p_term_only) / p_neither0)

if(frailty & marg_flag){
  #Following Xu (2010), marginal hazards reflect interplay of H1 and H2
  #this could be computed on log-scale, fyi! could be useful for some reason
  value$p_neither_marg <- exp( (-1/theta) * log1p(theta * (H_list[[1]] + H_list[[2]])) )
  #make it so that you've "started integrating" from time tseq[1]
  #in practice, this is no different if tseq = 0
  p_neither0_marg <- exp( (-1/theta) * log1p(theta * (H_list[[1]][1,] + H_list[[2]][1,])) )
  value$p_neither_marg <- t(t(value$p_neither_marg) / p_neither0_marg)

  H1_b <- H_list[[1]][-1,,drop=FALSE]
  H1_a <- H_list[[1]][-length(tseq),,drop=FALSE]
  H2_b <- H_list[[2]][-1,,drop=FALSE]
  H2_a <- H_list[[2]][-length(tseq),,drop=FALSE]

  #Here, we interchange the integral of gamma with the summation of trapezoid rule
  #to get sum of "integrated integrands"
  value$p_term_only_marg <- rbind(0, 0.5 * apply(
    (1 + theta * (H1_a + H2_a))^(-1/theta) +
      (1 + theta * (H1_b + H2_a))^(-1/theta) -
      (1 + theta * (H1_a + H2_b))^(-1/theta) -
      (1 + theta * (H1_b + H2_b))^(-1/theta), MARGIN=2, cumsum))
  value$p_nonterm_marg <- rbind(0, 0.5 * apply(
    (1 + theta * (H1_a + H2_a))^(-1/theta) +
      (1 + theta * (H1_a + H2_b))^(-1/theta) -
      (1 + theta * (H1_b + H2_a))^(-1/theta) -
      (1 + theta * (H1_b + H2_b))^(-1/theta), MARGIN=2, cumsum))

  value$p_term_only_marg <- t(t(value$p_term_only_marg) / p_neither0_marg)
  value$p_nonterm_marg <- t(t(value$p_nonterm_marg) / p_neither0_marg)

  # #Below is a log-scale implementation of these same computations,
  # #I'm sure we could do this throughout but I'm not sure how essential it is...
  # value$p_term_only_marg <- 0.5 * exp(apply(
  #   logdiffexp(
  #     logsumexp((-1/theta)*log1p(theta * (H1_a + H2_a)),
  #                  (-1/theta)*log1p(theta * (H1_b + H2_a))),
  #     logsumexp((-1/theta)*log1p(theta * (H1_a + H2_b)),
  #                  (-1/theta)*log1p(theta * (H1_b + H2_b)))
  #   ), MARGIN=2, logcumsumexp))
  # value$p_term_only_marg <- 0.5 * exp(apply(
  #   logdiffexp(
  #     logsumexp((-1/theta)*log1p(theta * (H1_a + H2_a)),
  #                  (-1/theta)*log1p(theta * (H1_a + H2_b))),
  #     logsumexp((-1/theta)*log1p(theta * (H1_b + H2_a)),
  #                  (-1/theta)*log1p(theta * (H1_b + H2_b)))
  #   ), MARGIN=2, logcumsumexp))

}

#now compute SEs for these quantities
#loop through "individuals" (I know this is woefully inefficient but idk what to say...)
if(se_fit_flag){
  Finv_sub12 <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                       nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2]),
                     c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                       nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2])]

  #while I'm here, set up the matrix for the more complicated  probabilities in the loop far below
  Finv_sub123 <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                        nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                        nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv)),
                      c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                        nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                        nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv))]

  if(frailty & marg_flag){
    Finv_sub12_marg <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                              nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                              nP0_tot), #add frailty
                            c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                              nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                              nP0_tot)] #add frailty

    #while I'm here, set up the matrix for the more complicated probabilities in the loop far below
    Finv_sub123_marg <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                               nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                               nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv),
                               nP0_tot),
                             c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                               nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                               nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv),
                               nP0_tot)]
  }
  value$se_p_neither_marg <-
    value$se_p_neither <- value$se_logp_neither <-
    value$se_p_nonterm <- value$se_p_term_only <-
    value$se_p_nonterm_marg <- value$se_p_term_only_marg <-
    matrix(data=NA, nrow=NROW(H_list[[1]]),ncol=NCOL(H_list[[1]]))

  #again set up containers for the p values to be computed in the loop far below
  value$se_p_both_marg <- value$se_p_nonterm_only_marg <-
    value$se_p_both <- value$se_p_nonterm_only <-
    matrix(data = NA, nrow = length(tseq), ncol=NCOL(value$H1))


  basis1 <- get_basis(y = tseq,knots_vec = knots_list[[1]],hazard = hazard,deriv = FALSE)
  dbasis1 <- get_basis(y = tseq,knots_vec = knots_list[[1]],hazard = hazard,deriv = TRUE)
  basis1_quad <- get_basis(y=quad_points, knots_vec=knots_list[[1]],hazard=hazard,deriv = FALSE)
  basis2 <- get_basis(y = tseq,knots_vec = knots_list[[2]],hazard = hazard,deriv = FALSE)
  dbasis2 <- get_basis(y = tseq,knots_vec = knots_list[[2]],hazard = hazard,deriv = TRUE)
  basis2_quad <- get_basis(y=quad_points, knots_vec=knots_list[[2]],hazard=hazard,deriv = FALSE)


  # here is the problem
  for(j in 1:max(1, NCOL(value$H1))){
    dH1b <- get_jac(tseq=tseq,
                    xnew=x1new[j,,drop=FALSE],
                    beta=beta_list[[1]],
                    eta=eta_list[[1]][j], #WILL THIS CAUSE A PROBLEM IF ONE X HAS NO ELEMENTS BUT ANOTHER HAS SEVERAL?
                    hazard=hazard,
                    phi=para[nP0_start[1]:nP0_end[1]],
                    basis=basis1,
                    dbasis=dbasis1,
                    basis_quad=basis1_quad,
                    quad_weights=quad_weights,
                    func_type="H",
                    H = H_list[[1]][,j],
                    log_out = FALSE )
    dH2b <- get_jac(tseq=tseq, xnew=x2new[j,,drop=FALSE], beta=beta_list[[2]],
                    eta=eta_list[[2]][j], #WILL THIS CAUSE A PROBLEM IF ONE X HAS NO ELEMENTS BUT ANOTHER HAS SEVERAL?
                    hazard=hazard, phi=para[nP0_start[2]:nP0_end[2]],
                    basis=basis2, dbasis=dbasis2,
                    basis_quad=basis2_quad, quad_weights=quad_weights,
                    func_type="H", H = H_list[[2]][,j],
                    log_out = FALSE )

    # does this cause a problem depending on left truncation??
    if(tseq[1]==0){
      dH1b[1,] <- 0
      dH2b[1,] <- 0
    }
    dH1a <- dH1b[-NROW(dH1b),,drop=FALSE]
    dH2a <- dH2b[-NROW(dH2b),,drop=FALSE]
    dH1b <- dH1b[-1,,drop=FALSE]
    dH2b <- dH2b[-1,,drop=FALSE]

    #see scratchwork for this derivation, basically it's leveraging the summation of trapezoid rule
    #for p_term_only (aka CIF_t)
    #apply(S1_a*S2_a + S1_b*S2_a - S1_a*S2_b - S1_b*S2_b, MARGIN = 2, cumsum)
    J <- 0.5 * cbind(
      dH1a*S1_a[,j]*(S2_a[,j] - S2_b[,j]) + dH1b*S1_b[,j]*(S2_a[,j] - S2_b[,j]),
      dH2a*S2_a[,j]*(S1_a[,j] + S1_b[,j]) - dH2b*S2_b[,j]*(S1_a[,j] + S1_b[,j])
    )
    value$se_p_term_only[,j] <- c(0,sqrt(diag(apply(apply(J %*% Finv_sub12 %*% t(J),1,cumsum),1,cumsum))))

    #for p_nonterm (aka CIF_nt)
    #apply(S2_a*S1_a + S2_b*S1_a - S2_a*S1_b - S2_b*S1_b, MARGIN = 2, cumsum)
    J <- 0.5 * cbind(
      dH1a*S1_a[,j]*(S2_a[,j] + S2_b[,j]) - dH1b*S1_b[,j]*(S2_a[,j] + S2_b[,j]),
      dH2a*S2_a[,j]*(S1_a[,j] - S1_b[,j]) + dH2b*S2_b[,j]*(S1_a[,j] - S1_b[,j])
    )
    value$se_p_nonterm[,j] <- c(0, sqrt(diag(apply(apply(J %*% Finv_sub12 %*% t(J),1,cumsum),1,cumsum))))

    #SE of "neither" probability
    J <- -cbind(dH1b, dH2b)
    value$se_logp_neither[,j] <- c(0, sqrt(rowSums(J * (J %*% Finv_sub12))))
    J <- value$p_neither[-1,j] * J
    value$se_p_neither[,j] <- c(0, sqrt(rowSums(J * (J %*% Finv_sub12))))

    if(frailty & marg_flag){
      #for p_term_only (aka CIF_t)
      #apply(S1_a*S2_a + S1_b*S2_a - S1_a*S2_b - S1_b*S2_b, MARGIN = 2, cumsum)
      J <- 0.5 * cbind(
        dH1a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) -
                  (1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1)) +
          dH1b * ((1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
        dH2a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) +
                  (1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1)) -
          dH2b * ((1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1) +
                    (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
        -marg_logtheta_deriv(a = H1_a[,j] + H2_a[,j], ltheta = log(theta)) +
          -marg_logtheta_deriv(a = H1_b[,j] + H2_a[,j], ltheta = log(theta)) -
          -marg_logtheta_deriv(a = H1_a[,j] + H2_b[,j], ltheta = log(theta)) -
          -marg_logtheta_deriv(a = H1_b[,j] + H2_b[,j], ltheta = log(theta)))
      value$se_p_term_only_marg[,j] <- c(0, sqrt(diag(apply(apply(J %*% Finv_sub12_marg %*% t(J),1,cumsum),1,cumsum))))

      #for p_nonterm (aka CIF_nt)
      #apply(S2_a*S1_a + S2_b*S1_a - S2_a*S1_b - S2_b*S1_b, MARGIN = 2, cumsum)
      J <- 0.5 * cbind(
        dH1a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) +
                  (1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1)) -
          dH1b * ((1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1) +
                    (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
        dH2a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) -
                  (1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1)) +
          dH2b * ((1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
        -marg_logtheta_deriv(a = H1_a[,j] + H2_a[,j], ltheta = log(theta)) +
          -marg_logtheta_deriv(a = H1_a[,j] + H2_b[,j], ltheta = log(theta)) -
          -marg_logtheta_deriv(a = H1_b[,j] + H2_a[,j], ltheta = log(theta)) -
          -marg_logtheta_deriv(a = H1_b[,j] + H2_b[,j], ltheta = log(theta)))
      value$se_p_nonterm_marg[,j] <- c(0, sqrt(diag(apply(apply(J %*% Finv_sub12_marg %*% t(J),1,cumsum),1,cumsum))))

      #SE of marginal "neither" probability
      J <- cbind(dH1b * (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1),
                 dH2b * (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1),
                 -marg_logtheta_deriv(a = H1_b[,j] + H2_b[,j], ltheta = log(theta)))
      value$se_p_neither_marg[,j] <- c(0, sqrt(rowSums(J * (J %*% Finv_sub12_marg))))
    }
  }
}

# browser()

#now, to handle the quantities that depend on h3...

#I think now, we'll have to use loops to incorporate h3 correctly
value$p_both_marg <- value$p_nonterm_only_marg <-
  value$p_both <- value$p_nonterm_only <-
  matrix(data = NA, nrow = length(tseq),
         ncol=NCOL(value$H1))

# if(abs(max(diff(tseq)) - min(diff(tseq))) > 1e-6){
#   warning("tseq points are not equally spaced. Some predicted functions may be incorrect.")
# }

value$p_nonterm_only[1,] <- 0
value$p_both[1,] <- 0
if(frailty & marg_flag){
  value$p_nonterm_only_marg[1,] <- 0
  value$p_both_marg[1,] <- 0
}

# #and further assume that t1 is not included as a covariate...
# if(tolower(model) %in% c("m","markov")){
#   H3_b <- value$H3[-1,,drop=FALSE]
#   H3_a <- value$H3[-length(tseq),,drop=FALSE]
#   value$p_nonterm_only <- exp(value$H3) * rbind(0, 0.5 *
#             apply(S2_a*exp(-H3_a)*S1_a + S2_b*exp(-H3_b)*S1_a -
#                     S2_a*exp(-H3_a)*S1_b - S2_b*exp(-H3_b)*S1_b, MARGIN = 2, cumsum))
# }



if(p3tv > 0){
  h3_tv_basis <- h3tv_basis_func(tseq)
  h3_tv_mult <- exp( as.vector(h3_tv_basis %*% h3tv_phi) )
} else{
  h3_tv_mult <- rep(1,length(tseq))
  h3_tv_basis <- NULL
}




for(i in 1:(length(tseq)-1)){

  # browser()

  #and further assume that t1 is not included as a covariate...
  if(tolower(model) %in% c("m","markov")){
    tseq_temp <- tseq[1:(i+1)]
    #last row of H3[1:i] is H3(t2i), so
    #compute matrix of H3(t2) - H3(t1) for every t1 between 0 and t1i
    #by "sweeping" out H3(t2i) from every row of -H3(t1) matrix

    # H3_a_sub <- t( -t(value$H3[1:i,,drop=FALSE]) + value$H3[i,])

    #YAY THIS UPDATE WORKS!! EVEN FOR LEFT TRUNCATED and IRREGULARLY SPACE DATA!
    H3_temp <- t( -t(H_list[[3]][1:(i+1),,drop=FALSE]) + H_list[[3]][(i+1),])
    H3_b_sub <- H3_temp[-1,,drop=FALSE]
    H3_a_sub <- H3_temp[-(i+1),,drop=FALSE]

  } else{
    tseq_temp <- tseq[(i+1)] - tseq[1:(i+1)]

    #semi-markov is trickier and I think may in fact need to be re-predicted in general
    #both to account for irregular spacing and for dependence on t_1

    #take the first i entries, and then reverse their order
    #to get a vector that increments from H3(t2i-0) to H3(t2i-t2i)

    #OK THIS WORKS FOR NON-TRUNCATED CASE! But not truncated case
    #because I think we'd need to compute new values of H3 regardless! i.e., H_3(5.5-5)
    #is not already computed if we are only computing H_3 on the range of 5 to 6
    # temp <- value$H3[1:(i+1),,drop=FALSE][(i+1):1,,drop=FALSE]
    # H3_b_sub <- temp[-1,,drop=FALSE]
    # H3_a_sub <- temp[-(i+1),,drop=FALSE]

    #reprediction solves the above issue
    #only lingering issue is how to bring in dependence on t1
    H3_temp <- pred_helper_uni(tseq = tseq_temp,
                               xnew = x3new,
                               phi = para[nP0_start[3]:nP0_end[3]],
                               beta = para[nP_start[3]:nP_end[3]],
                               hazard = hazard, knots_vec = knots_list[[3]],
                               n_quad = n_quad, quad_method = quad_method,
                               vcov_sub = NULL, offset = log(frailnew),
                               alpha=alpha)$H * h3_tv_mult[1:(i+1)] #NOTE COLUMNWISE MULTIPLICATION BY T1 EFFECT FACTOR

    #the final one should always be zero
    H3_temp[i+1,] <- 0
    H3_b_sub <- H3_temp[-1,,drop=FALSE]
    H3_a_sub <- H3_temp[-(i+1),,drop=FALSE]
  }

  S1_a_sub <- S1_a[1:i,,drop=FALSE]
  S1_b_sub <- S1_b[1:i,,drop=FALSE]
  integrand_a_sub <- S2_a[1:i,,drop=FALSE] * exp(-H3_a_sub)
  integrand_b_sub <- S2_b[1:i,,drop=FALSE] * exp(-H3_b_sub)
  value$p_nonterm_only[i+1,] <- 0.5 *
    colSums( (integrand_a_sub + integrand_b_sub) * (S1_a_sub - S1_b_sub) )

  integrand_a_sub <- S2_a[1:i,,drop=FALSE] * (1-exp(-H3_a_sub))
  integrand_b_sub <- S2_b[1:i,,drop=FALSE] * (1-exp(-H3_b_sub))
  value$p_both[i+1,] <- 0.5 *
    colSums( (integrand_a_sub + integrand_b_sub) * (S1_a_sub - S1_b_sub) )

  if(frailty & marg_flag){
    H1_a_sub <- H1_a[1:i,,drop=FALSE]
    H1_b_sub <- H1_b[1:i,,drop=FALSE]
    H2_a_sub <- H2_a[1:i,,drop=FALSE]
    H2_b_sub <- H2_b[1:i,,drop=FALSE]
    value$p_nonterm_only_marg[i+1,] <- 0.5 * colSums(
      (1 + theta * (H1_a_sub + H2_a_sub + H3_a_sub))^(-1/theta) +
        (1 + theta * (H1_a_sub + H2_b_sub + H3_b_sub))^(-1/theta) -
        (1 + theta * (H1_b_sub + H2_a_sub + H3_a_sub))^(-1/theta) -
        (1 + theta * (H1_b_sub + H2_b_sub + H3_b_sub))^(-1/theta))
  }

  #RETURN TO THESE ONCE I FINISH SOLVING THE MATTER OF THE SEMI-MARKOV!!
  if(se_fit_flag){
    #now, to generate standard errors for these predictions

    basis <- get_basis(y = tseq_temp,knots_vec = knots_list[[3]],hazard = hazard,deriv = FALSE)
    dbasis <- get_basis(y = tseq_temp,knots_vec = knots_list[[3]],hazard = hazard,deriv = TRUE)
    quad_points_sub <- transform_quad_points(n_quad = n_quad,
                                             quad_method=quad_method, a = 0,b = tseq_temp)
    basis_quad <- get_basis(y=quad_points_sub, knots_vec=knots_list[[3]], hazard = hazard,deriv = FALSE)

    dH1a_sub <- dH1a[1:i,,drop=FALSE]
    dH1b_sub <- dH1b[1:i,,drop=FALSE]
    dH2a_sub <- dH2a[1:i,,drop=FALSE]
    dH2b_sub <- dH2b[1:i,,drop=FALSE]

    # browser()

    for(j in 1:NCOL(value$H1)){
      #jacobian for H3 has to follow the same structure as H3 itself does above
      #this dH3a will be "updated" below absed on the model type
      dH3a <- get_jac(tseq=tseq_temp, xnew=x1new[j,,drop=FALSE], beta=beta_list[[3]],
                      eta=eta_list[[1]][j], #WILL THIS CAUSE A PROBLEM IF ONE X HAS NO ELEMENTS BUT ANOTHER HAS SEVERAL?
                      hazard=hazard, phi=para[nP0_start[3]:nP0_end[3]],
                      basis=basis, dbasis=dbasis,
                      basis_quad=basis_quad, quad_weights=quad_weights,
                      func_type="H", H = H3_temp[,j],
                      log_out = FALSE )
      if(tolower(model) %in% c("m","markov")){
        #last row of H3[1:i] is H3(t2i), so
        #compute matrix of H3(t2) - H3(t1) for every t1 between 0 and t1i
        #by "sweeping" out H3(t2i) from every row of H3(t1) matrix
        dH3a <- t( -t(dH3a) + dH3a[(i+1),])
      } else if(p3tv>0){
        #if semi-markov and some kind of t1-dependence, just include the corresponding columns here!
        dH3a <- cbind(dH3a, H3_temp[,j] * h3_tv_basis[1:(i+1),,drop=FALSE])
      }
      #we know that by definition, the last element of tseq_temp is supposed to be zero
      dH3a[(i+1),] <- 0

      dH3b <- dH3a[-1,,drop=FALSE]
      dH3a <- dH3a[-(i+1),,drop=FALSE]

      #for p_term_only (aka CIF_t)
      #apply(S1_a*S2_a + S1_b*S2_a - S1_a*S2_b - S1_b*S2_b, MARGIN = 2, cumsum)
      S2_a_sub <- S2_a[1:i,,drop=FALSE]
      S2_b_sub <- S2_b[1:i,,drop=FALSE]
      J <- 0.5 * cbind(
        dH1a_sub * (S1_a_sub[,j] * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) +
                                      S2_b_sub[,j] * exp(-H3_b_sub[,j]))) -
          dH1b_sub * (S1_b_sub[,j] * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) +
                                        S2_b_sub[,j] * exp(-H3_b_sub[,j]))),
        dH2a_sub * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) * (S1_a_sub[,j] - S1_b_sub[,j])) +
          dH2b_sub * (S2_b_sub[,j] * exp(-H3_b_sub[,j]) * (S1_a_sub[,j] -S1_b_sub[,j])),
        dH3a * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) * (S1_a_sub[,j] - S1_b_sub[,j])) +
          dH3b * (S2_b_sub[,j] * exp(-H3_b_sub[,j]) * (S1_a_sub[,j] - S1_b_sub[,j])))
      value$se_p_nonterm_only[i+1,j] <- sqrt(sum(J %*% Finv_sub123 %*% t(J)))

      J <- 0.5 * cbind(
        dH1a_sub * (S1_a_sub[,j] * (S2_a_sub[,j] * (1-exp(-H3_a_sub[,j])) +
                                      S2_b_sub[,j] * (1-exp(-H3_b_sub[,j])))) -
          dH1b_sub * (S1_b_sub[,j] * (S2_a_sub[,j] * (1-exp(-H3_a_sub[,j])) +
                                        S2_b_sub[,j] * (1-exp(-H3_b_sub[,j])))),
        dH2a_sub * (S2_a_sub[,j] * ((1-exp(-H3_a_sub[,j])) * S1_a_sub[,j] -
                                      (1-exp(-H3_a_sub[,j])) * S1_b_sub[,j])) +
          dH2b_sub * (S2_b_sub[,j] * ((1-exp(-H3_b_sub[,j])) * S1_a_sub[,j] -
                                        (1-exp(-H3_b_sub[,j])) * S1_b_sub[,j])),
        dH3a * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) * (S1_b_sub[,j] - S1_a_sub[,j])) +
          dH3b * (S2_b_sub[,j] * exp(-H3_b_sub[,j]) * (S1_b_sub[,j] - S1_a_sub[,j])))
      value$se_p_both[i+1,j] <- sqrt(sum(J %*% Finv_sub123 %*% t(J)))

      #last piece, se for marginal probabilities of "both" and "nonterm only"
      if(frailty & marg_flag){
        #for p_nonterm (aka CIF_nt)
        #apply(S2_a*S1_a + S2_b*S1_a - S2_a*S1_b - S2_b*S1_b, MARGIN = 2, cumsum)
        J <- 0.5 * cbind(
          dH1a_sub * ((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) +
                        (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)) -
            dH1b_sub * ((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) +
                          (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
          dH2a_sub * ((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) -
                        (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
            dH2b_sub * ((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1) -
                          (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
          dH3a * ((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
            dH3b * ((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1) -
                      (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
          -marg_logtheta_deriv(a = H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) +
            -marg_logtheta_deriv(a = H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)) -
            -marg_logtheta_deriv(a = H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) -
            -marg_logtheta_deriv(a = H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)))
        value$se_p_nonterm_only_marg[i+1,j] <- sqrt(sum(J %*% Finv_sub123_marg %*% t(J)))

        #for p_both
        J <- 0.5 * cbind(
          dH1a_sub * (((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
                        ((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                           (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1))) -
            dH1b_sub * (((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                           (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
                          ((1 + theta * (H1_b_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                             (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1))),
          dH2a_sub * (((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) -
                        ((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                           (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1))) +
            dH2b_sub * (((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                           (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)) -
                          ((1 + theta * (H1_b_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                             (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1))),
          dH3a * ((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
            dH3b * ((1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1) -
                      (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
          (marg_logtheta_deriv(a = H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) -
             marg_logtheta_deriv(a = H1_a_sub[,j] + H2_a_sub[,j], ltheta = log(theta))) +
            (marg_logtheta_deriv(a = H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)) -
               marg_logtheta_deriv(a = H1_a_sub[,j] + H2_b_sub[,j], ltheta = log(theta))) -
            (marg_logtheta_deriv(a = H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) -
               marg_logtheta_deriv(a = H1_b_sub[,j] + H2_a_sub[,j], ltheta = log(theta))) -
            (marg_logtheta_deriv(a = H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)) -
               marg_logtheta_deriv(a = H1_b_sub[,j] + H2_b_sub[,j], ltheta = log(theta))))
        value$se_p_both_marg[i+1,j] <- sqrt(sum(J %*% Finv_sub123_marg %*% t(J)))
      }
    }

  }

}

value$p_nonterm_only <- t( t(value$p_nonterm_only) / p_neither0)
value$p_both <- t( t(value$p_both) / p_neither0)
if(frailty & marg_flag){
  value$p_nonterm_only_marg <- t( t(value$p_nonterm_only_marg) / p_neither0_marg)

  #The marginal expression for "both" probability is literally
  #the difference of these two quantities...
  value$p_both_marg <- value$p_nonterm_marg - value$p_nonterm_only_marg

}

value
