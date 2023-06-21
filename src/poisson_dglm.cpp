#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List update_poisson_dglm(int y, arma::vec F, arma::mat G, arma::mat D,
                               arma::vec a, arma::mat R){

  // Prior moments for the linear predictor
  double f = arma::as_scalar(F.t() * a);
  double q = arma::as_scalar(F.t() * R * F);

  // Conjugate prior parameters for the poisson (1st order approximation of digamma)
  double parm1 = 1 / q;
  double parm2 = std::exp(-f)  / q;

  // Update the posterior moments of linear predictor (2nd order approximation of digamma)
  double f_s = R::digamma(parm1 + y) - std::log(parm2 + 1);
//    std::log((parm1 + y)  / (parm2 + 1)) + 1 / (2 * (parm1 + y));
  double q_s = R::trigamma(parm1 + y);
    // (2 * parm1 + 2 * y - 1) / (2 * std::pow((parm1 + y), 2));

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

  // Output
  return(Rcpp::List::create(
      Rcpp::Named("f")=f,
      Rcpp::Named("q")=q,
      Rcpp::Named("a")=a,
      Rcpp::Named("R")=R,
      Rcpp::Named("m")=m,
      Rcpp::Named("C")=C,
      Rcpp::Named("parm1")=parm1,
      Rcpp::Named("parm2")=parm2,
      Rcpp::Named("lpl")=0));
}

// [[Rcpp::export]]
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
  arma::vec f_seq(k, arma::fill::zeros);
  arma::vec q_seq(k, arma::fill::zeros);

  a_seq.col(0) = a;
  R_seq.slice(0) = R;

  //Rcpp::List tmp;
  for (int t = 0; t < k; ++t){
    Rcpp::List tmp = update_poisson_dglm(y[t], F, G, D, a_seq.col(t), R_seq.slice(t));

    parm1_seq(t) = tmp["parm1"];
    parm2_seq(t) = tmp["parm2"];
    f_seq(t) = tmp["f"];
    q_seq(t) = tmp["q"];
    m_seq.col(t) = Rcpp::as<arma::vec>(tmp["m"]);
    C_seq.slice(t) = Rcpp::as<arma::mat>(tmp["C"]);
    a_seq.col(t + 1) = Rcpp::as<arma::vec>(tmp["a"]);
    R_seq.slice(t + 1) = Rcpp::as<arma::mat>(tmp["R"]);
  }

  return(Rcpp::List::create(
      Rcpp::Named("a")=a_seq,
      Rcpp::Named("R")=R_seq,
      Rcpp::Named("m")=m_seq,
      Rcpp::Named("C")=C_seq,
      Rcpp::Named("f")=f_seq,
      Rcpp::Named("q")=q_seq,
      Rcpp::Named("parm1")=parm1_seq,
      Rcpp::Named("parm2")=parm2_seq,
      Rcpp::Named("loglik")=0));


}




