#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// NOT UP TO DATE CURRENTLY

// [[Rcpp::export]]
arma::vec GradParamLongitMTLogLik(const arma::vec& par,
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

  // Extract unique times and IDs
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

  // Precompute terms
  arma::vec XBeta_Long = Xlong * beta;
  arma::vec XBeta_T = XT * xi;
  arma::vec YTheta_T = XTglobdep * theta;

  // Initialize gradient matrix
  arma::mat GradLogLik(n_sample, n_par_long + n_par_T, fill::zeros);

  for (int i = 0; i < n_sample; i++) {
    arma::uvec IDX = find(ID == ID_uqe(i));
    if (IDX.n_elem <= 1) continue; // Skip if only one observation

    arma::vec iTM = TM.elem(IDX.subvec(1, IDX.n_elem - 1)) - 1;
    arma::vec unique_TM = unique(iTM);

    // Extract subject-specific data
    arma::vec iY_long = Ylong.elem(IDX.subvec(1, IDX.n_elem - 1));
    arma::vec iY_T = YT.elem(IDX.subvec(1, IDX.n_elem - 1));

    arma::mat iXlongAR = XlongAR.rows(IDX.subvec(1, IDX.n_elem - 1));
    arma::mat iXlongTrend = XlongTrend.rows(IDX.subvec(1, IDX.n_elem - 1));
    arma::vec iYTheta_T = YTheta_T.elem(IDX.subvec(1, IDX.n_elem - 1));
    arma::mat iXTlocdep = XTlocdep.rows(IDX.subvec(1, IDX.n_elem - 1));
    arma::mat iXTglobdep = XTglobdep.rows(IDX.subvec(1, IDX.n_elem - 1));

    // Linear predictors
    arma::vec iLinPred_long = iXlongTrend * zeta + XBeta_Long(ID_uqe(i)) +
      sum(iXlongAR.each_col() % eta, 1);

    arma::vec iLinPred_T = lambda + XBeta_T(ID_uqe(i)) + iYTheta_T +
      sum(iXTlocdep.each_col() % varphi, 1);

    arma::vec iProb_T = exp(iLinPred_T) / (1 + exp(iLinPred_T));

    double VarFac = 1 / (2 * sigma_sq);

    // Gradient computation
    arma::mat iGradLogLik(n_par_long + n_par_T, unique_TM.n_elem, fill::zeros);

    for (unsigned int k_idx = 0; k_idx < unique_TM.n_elem; k_idx++) {
      double k_now = unique_TM(k_idx);
      arma::uvec IDX_now = find(iTM == k_now);

      arma::vec iY_long_now = iY_long.elem(IDX_now);
      arma::vec iY_T_now = iY_T.elem(IDX_now);
      arma::mat iXlongAR_now = iXlongAR.rows(IDX_now);
      arma::mat iXlongTrend_now = iXlongTrend.rows(IDX_now);
      arma::mat iXTlocdep_now = iXTlocdep.rows(IDX_now);
      arma::mat iXTglobdep_now = iXTglobdep.rows(IDX_now);

      arma::vec iLinPred_long_now = iLinPred_long.elem(IDX_now);
      arma::vec iProb_T_now = iProb_T.elem(IDX_now);

      iGradLogLik.submat(0, k_idx, p_long_Trend - 1, k_idx) =
        VarFac * (iXlongTrend_now.t() * ((1 - iY_T_now) % (iY_long_now - iLinPred_long_now)));

      iGradLogLik.submat(p_long_Trend, k_idx, p_long_Trend + p_long_AR - 1, k_idx) =
        VarFac * (iXlongAR_now.t() * ((1 - iY_T_now) % (iY_long_now - iLinPred_long_now)));

      iGradLogLik(n_par_long - 1, k_idx) =
        -VarFac + (1 / (2 * pow(sigma_sq, 2))) * sum((1 - iY_T_now) % pow(iY_long_now - iLinPred_long_now, 2));

        iGradLogLik(n_par_long, k_idx) = sum(iY_T_now - iProb_T_now);

        iGradLogLik.submat(n_par_long + 1, k_idx, n_par_long + p_glob_dep, k_idx) =
          iXTglobdep_now.t() * (iY_T_now - iProb_T_now);

        iGradLogLik.submat(n_par_long + p_glob_dep + 1, k_idx, n_par_long + p_glob_dep + p_loc_dep, k_idx) =
          iXTlocdep_now.t() * (iY_T_now - iProb_T_now);
    }

    GradLogLik.row(i) = sum(iGradLogLik, 1).t();
  }

  return sum(GradLogLik, 0).t();
}
