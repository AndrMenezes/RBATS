#include <RcppArmadillo.h>

// [[Rcpp::export]]
void update_poisson_dglm(int y, arma::vec F, arma::mat G, arma::mat D,
                         arma::vec &a, arma::mat &R,
                         double *parm1, double *parm2){

  // Prior moments for the linear predictor
  double f = arma::as_scalar(F.t() * a);
  double q = arma::as_scalar(F.t() * R * F);

  // Conjugate prior parameters for the poisson (approximation)
  *parm1 = 1 / q;
  *parm2 = std::exp(-f)  / q;

  // Update the posterior moments of linear predictor
  double f_s = std::log((*parm1 + y)  / (*parm2 + 1)) + 1 / (2 * (*parm1 + y));
  double q_s = (2 * *parm1 + 2 * y - 1) / (2 * (*parm1 + y));

  // Linear Bayes update for the state parameters
  arma::vec m = a + R * F * (f_s - f) / q;
  arma::mat C = R - R * F * F.t() * R * (1 - q_s / q) / q;

  // Evolution step: priors for time t + 1 from the posteriors m, C
  a = G * m;
  R = G * C * G.t();
  // Make sure the covariance matrix is symmetric
  R = (R + R.t()) / 2;

  // Discount information
  R = D % R;

  std::cout << a << std::endl;
  std::cout << R << std::endl;
  std::cout << *parm2 << std::endl;
  std::cout << *parm2 << std::endl;
  // std::cout << *m << std::endl;
  // std::cout << *C << std::endl;

}

Rcpp::List forward_filter_poisson_dglm(arma::vec y, arma::vec F, arma::mat G,
                                       arma::mat D,
                                       arma::vec a, arma::mat R){


  const int k = y.size();
  arma::cube R_seq(R.n_rows, R.n_cols, k + 1, arma::fill::zeros);
  arma::cube C_seq(R.n_rows, R.n_cols, k, arma::fill::zeros);
  arma::mat a_seq(a.size(), k + 1, arma::fill::zeros);
  arma::mat m_seq(a.size(), k, arma::fill::zeros);
  arma::vec parm1_seq(k, arma::fill::zeros);
  arma::vec parm2_seq(k, arma::fill::zeros);

  double *parm1, *parm2;
  a_seq.col(0) = a;
  R_seq.slice(0) = R;
  arma::vec m = a;
  arma::mat C = R;

  for (int t = 0; t < k; ++t){
    update_poisson_dglm(y[t], F, G, D, a, R, m, C, *parm1, *parm2);

    parm1_seq(t) = *parm1;
    parm2_seq(t) = *parm1;
    m_seq.col(t) = m;
    C_seq.slice(t) = C;
    a_seq.col(t + 1) = a;
    R_seq.slice(t + 1) = R;
  }


}




