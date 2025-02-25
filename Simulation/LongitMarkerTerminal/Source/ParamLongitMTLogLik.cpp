#include <RcppArmadillo.h>
#include <iostream>
#include <map>
#include <string>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]

// NOT UP TO DATE CURRENTLY

// [[Rcpp::export]]
double ParamLongitMTLogLik(const arma::vec& par,
                           const arma::vec& ID,
                           const arma::vec& TM,
                           const arma::mat& Xlong,
                           const arma::mat& XlongAR,
                           const arma::mat& XlongTrend,
                           const arma::mat& XT,
                           const arma::mat& XTglobdep,
                           const arma::mat& XTlocdep,
                           const arma::vec& Ylong,
                           const arma::vec& YT,
                           const arma::vec& riskT) {

  // Extract unique times and subjects
  arma::vec times = unique(TM);
  arma::vec ID_uqe = unique(ID);
  int n_sample = ID_uqe.n_elem;

  // Parameter dimensions
  int p_long = Xlong.n_cols;
  int p_long_AR = XlongAR.n_cols;
  int p_long_Trend = XlongTrend.n_cols;
  int p_T = XT.n_cols;
  int p_glob_dep = XTglobdep.n_cols;
  int p_loc_dep = XTlocdep.n_cols;

  int n_par_long = p_long + p_long_AR + p_long_Trend + 1;
  int n_par_T = 1 + p_T + p_glob_dep + p_loc_dep;

  // Extract parameters correctly
  arma::vec zeta = par.subvec(0, p_long_Trend - 1);
  arma::vec eta = par.subvec(p_long_Trend, p_long_Trend + p_long_AR - 1);
  arma::vec beta = par.subvec(p_long_Trend + p_long_AR, p_long_Trend + p_long_AR + p_long - 1);
  double sigma_sq = par(n_par_long - 1);
  double lambda = par(n_par_long);
  arma::vec theta = par.subvec(n_par_long + 1, n_par_long + p_glob_dep);
  arma::vec varphi = par.subvec(n_par_long + p_glob_dep + 1, n_par_long + p_glob_dep + p_loc_dep);
  arma::vec xi = par.subvec(n_par_long + p_glob_dep + p_loc_dep + 1, n_par_long + n_par_T - 1);

  // Precompute fixed terms
  arma::vec XBeta_Long = Xlong * beta;
  arma::vec XBeta_T = XT * xi;
  arma::vec YTheta_T = XTglobdep * theta;

  // Initialize LogLik vector
  arma::vec LogLik(n_sample, fill::zeros);

  for (int i = 0; i < n_sample; i++) {
    arma::uvec IDX = find(ID == ID_uqe(i));
    if (IDX.n_elem <= 1) continue; // Skip if only one observation

    arma::vec iTM = TM.elem(IDX.subvec(1, IDX.n_elem - 1)) - 1;

    // Extract subject-specific data
    arma::vec iY_long = Ylong.elem(IDX.subvec(1, IDX.n_elem - 1));
    arma::vec iY_T = YT.elem(IDX.subvec(1, IDX.n_elem - 1));

    arma::mat iXlongAR = XlongAR.rows(IDX.subvec(1, IDX.n_elem - 1));
    arma::mat iXlongTrend = XlongTrend.rows(IDX.subvec(1, IDX.n_elem - 1));
    arma::vec iYTheta_T = YTheta_T.elem(IDX.subvec(1, IDX.n_elem - 1));
    arma::mat iXTlocdep = XTlocdep.rows(IDX.subvec(1, IDX.n_elem - 1));

    // Compute linear predictors
    arma::vec iLinPred_long = iXlongTrend * zeta + XBeta_Long(ID_uqe(i)) + iXlongAR * eta;

    arma::vec iLinPred_T = lambda + XBeta_T(ID_uqe(i)) + iYTheta_T + iXTlocdep * varphi;

    arma::vec iProb_T = exp(iLinPred_T) / (1 + exp(iLinPred_T));

    // Compute log-likelihood using vectorized operations
    arma::vec iLogLik = (1 - iY_T) * (-0.5 * log(2 * datum::pi * sigma_sq)) +
      (-0.5 / sigma_sq) * pow(iY_long - iLinPred_long, 2) +
      (1 - iY_T) % (1 - iProb_T) +
      iY_T % iProb_T;


    // Sum up log-likelihood for this subject
    LogLik(i) = sum(iLogLik);
  }

  return -sum(LogLik);  // Negative log-likelihood for minimization
}


// [[Rcpp::export]]
Rcpp::List InvestigateParamLongitMTLogLik(const arma::vec& par,
                                          const arma::vec& ID,
                                          const arma::vec& TM,
                                          const arma::mat& Xlong,
                                          const arma::mat& XlongAR,
                                          const arma::mat& XlongTrend,
                                          const arma::mat& XT,
                                          const arma::mat& XTglobdep,
                                          const arma::mat& XTlocdep,
                                          const arma::vec& Ylong,
                                          const arma::vec& YT,
                                          const arma::vec& riskT) {

  // Extract unique times and subjects
  arma::vec times = unique(TM);
  arma::vec ID_uqe = unique(ID);
  int n_sample = ID_uqe.n_elem;

  // Parameter dimensions
  int p_long = Xlong.n_cols;
  int p_long_AR = XlongAR.n_cols;
  int p_long_Trend = XlongTrend.n_cols;
  int p_T = XT.n_cols;
  int p_glob_dep = XTglobdep.n_cols;
  int p_loc_dep = XTlocdep.n_cols;

  int n_par_long = p_long + p_long_AR + p_long_Trend + 1;
  int n_par_T = 1 + p_T + p_glob_dep + p_loc_dep;

  // Extract parameters correctly
  arma::vec zeta = par.subvec(0, p_long_Trend - 1);
  arma::vec eta = par.subvec(p_long_Trend, p_long_Trend + p_long_AR - 1);
  arma::vec beta = par.subvec(p_long_Trend + p_long_AR, p_long_Trend + p_long_AR + p_long - 1);
  double sigma_sq = par(n_par_long - 1);
  double lambda = par(n_par_long);
  arma::vec theta = par.subvec(n_par_long + 1, n_par_long + p_glob_dep);
  arma::vec varphi = par.subvec(n_par_long + p_glob_dep + 1, n_par_long + p_glob_dep + p_loc_dep);
  arma::vec xi = par.subvec(n_par_long + p_glob_dep + p_loc_dep + 1, n_par_long + n_par_T - 1);

  // Extract subject-specific data safely
  arma::vec XBeta_Long = Xlong * beta;
  arma::vec XBeta_T = XT * xi;
  arma::vec YTheta_T = XTglobdep * theta;

  arma::uvec IDX = find(ID == ID_uqe(1));
  if (IDX.n_elem <= 1) return Rcpp::List::create(); // Handle edge case

  arma::vec iTM = TM.elem(IDX.subvec(1, IDX.n_elem - 1)) - 1;
  arma::vec iY_long = Ylong.elem(IDX.subvec(1, IDX.n_elem - 1));
  arma::vec iY_T = YT.elem(IDX.subvec(1, IDX.n_elem - 1));
  arma::mat iXlongAR = XlongAR.rows(IDX.subvec(1, IDX.n_elem - 1));
  arma::mat iXlongTrend = XlongTrend.rows(IDX.subvec(1, IDX.n_elem - 1));
  arma::vec iYTheta_T = YT.elem(IDX.subvec(1, IDX.n_elem - 1));
  arma::mat iXTlocdep = XTlocdep.rows(IDX.subvec(1, IDX.n_elem - 1));

  // Compute linear predictors
  arma::vec iLinPred_long = iXlongTrend * zeta + XBeta_Long(ID_uqe(1)) +
    iXlongAR * eta;

  arma::vec iLinPred_T = lambda + XBeta_T(ID_uqe(1)) + iYTheta_T + iXTlocdep * varphi;


  // Store dimensions in an Rcpp::List
  Rcpp::List dim_list = Rcpp::List::create(
    _["zeta"] = zeta.n_elem,
    _["eta"] = eta.n_elem,
    _["beta"] = beta.n_elem,
    _["sigma_sq"] = 1,
    _["lambda"] = 1,
    _["theta"] = theta.n_elem,
    _["varphi"] = varphi.n_elem,
    _["xi"] = xi.n_elem,

    _["XBetaLong"] = Rcpp::IntegerVector::create(XBeta_Long.n_rows, XBeta_Long.n_cols),
    _["XBetaT"] = Rcpp::IntegerVector::create(XBeta_T.n_rows, XBeta_T.n_cols),
    _["YthetaT"] = Rcpp::IntegerVector::create(YTheta_T.n_rows, YTheta_T.n_cols),

    _["iTM"] = iTM.n_elem,
    _["iY_long"] = iY_long.n_elem,
    _["iY_T"] = iY_T.n_elem,
    _["iXlongAR"] = Rcpp::IntegerVector::create(iXlongAR.n_rows, iXlongAR.n_cols),
    _["iXlongTrend"] = Rcpp::IntegerVector::create(iXlongTrend.n_rows, iXlongTrend.n_cols),
    _["iYTheta_T"] = iYTheta_T.n_elem,
    _["iXTlocdep"] = Rcpp::IntegerVector::create(iXTlocdep.n_rows, iXTlocdep.n_cols),

    _["LinPredLong"] = Rcpp::IntegerVector::create(iLinPred_long.n_rows, iLinPred_long.n_cols),
    _["LinPredT"] = Rcpp::IntegerVector::create(iLinPred_T.n_rows, iLinPred_T.n_cols)
  );

  return dim_list;
}


