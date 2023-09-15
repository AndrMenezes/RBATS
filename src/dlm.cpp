#include <RcppArmadillo.h>

void evolve_prior(arma::vec &a, arma::mat &R,
                  arma::vec &m, arma::mat &C,
                  arma::mat &G, arma::mat &D,
                  double &s,
                  int &n_parms,
                  int &ar_order,
                  int &tf_order,
                  Rcpp::IntegerVector &i_ar,
                  Rcpp::IntegerVector &i_tf,
                  double &x_tf
                    ) {

  a = m;

  // Check if there is linear components, then evolve its prior mean
  if (n_parms > 0) {
    // Prior evolution for the linear components
    a.subvec(0, n_parms - 1) = G.submat(0, 0, n_parms - 1, n_parms - 1) * m.subvec(0, n_parms - 1);
  }

  // Evolution for the autoregressive components
  if (ar_order > 0) {
    // First component of AR ( sum_{i=1}^p phi_{i} xi_i )
    a.subvec(i_ar(0), i_ar(0)) = (
      m.subvec(i_ar(0), i_ar(ar_order - 1)).t()
    * m.subvec(i_ar(ar_order), i_ar(2 * ar_order - 1)));

    // Changing the elements of xi_{2, ..., order} elements
    for (int j=0; j < ar_order - 1; j++) {
      a(i_ar(j + 1)) = m(i_ar(j));
    }

    // Changing the nonlinear part of G
    for (int j=0; j < ar_order; j++) {
      // xi_1 to xi_order
      G(i_ar(0), i_ar(j)) = m(i_ar(j + ar_order));
      // phi_1 to phi_order
      G(i_ar(0), i_ar(j + ar_order)) = m(i_ar(j));
    }

  }

  // Evolution for the transfer function components
  if (tf_order > 0) {
    // First component of AR ( sum_{i=1}^p lambda_{i} E_i + psi * x)
    a.subvec(i_tf(0), i_tf(0)) = (
      m.subvec(i_tf(0), i_tf(tf_order - 1)).t()
    * m.subvec(i_tf(tf_order), i_tf(2 * tf_order - 1))
    + m(i_tf(2 * tf_order)) * x_tf);

    // Changing the elements of E_{2, ..., order} elements
    for (int j=0; j < tf_order - 1; j++) {
      a(i_tf(j + 1)) = m(i_tf(j));
    }

    // Changing the nonlinear part of G
    for (int j=0; j < tf_order; j++) {
      // E_1 to E_order
      G(i_tf(0), i_tf(j)) = m(i_tf(j + tf_order));
      // lambda_1 to lambda_order
      G(i_tf(0), i_tf(j + tf_order)) = m(i_tf(j));
    }
    // derivative w.r.t psi is X_t
    G(i_tf(0), i_tf(2 * tf_order)) = x_tf;

  }

  // Prior covariance matrix
  R = G * C * G.t();

  // Discount information
  R %= D;

  if (ar_order > 0) {
    R(i_ar(0), i_ar(0)) += s;
  }

  // Make sure the covariance matrix is symmetric
  R += R.t();
  R /= 2;

}

void forecast_marginal(double &f, double &q,
                       arma::vec &F, arma::vec &a, arma::mat &R, double &s,
                       int &ar_order) {

  // 1-step ahead predictive moments
  f = arma::as_scalar(F.t() * a);
  q = arma::as_scalar(F.t() * R * F);

  if (ar_order == 0) {
    q += s;
  }

}

void update_posterior(arma::vec &m, arma::mat &C,
                      double &n, double &s,
                      double &y, double &f, double &q,
                      arma::vec &a, arma::mat &R,
                      arma::vec &F, double &df_variance) {

  m = a;
  C = R;
  // If the data is not missing update the model
  if (!std::isnan(y)) {

    double e = y - f;
    arma::vec A = (R * F) / q;
    double r = (n + e * e / q) / (n + 1);

    // Kalman filter update
    m = a + A * e;
    C = r * (R - q * A * A.t());
    s *= r;
    n++;
    n *= df_variance;
  }


}


// [[Rcpp::export]]
Rcpp::List forward_filter_dlm(arma::vec y,
                              arma::mat F,
                              arma::mat G,
                              arma::mat D,
                              arma::vec m,
                              arma::mat C,
                              double n,
                              double s,
                              double df_variance,
                              int n_parms,
                              int ar_order,
                              int tf_order,
                              Rcpp::IntegerVector i_ar,
                              Rcpp::IntegerVector i_tf,
                              Rcpp::NumericVector xreg_tf){

  const int k = y.size();

  arma::cube R_seq(C.n_rows, C.n_cols, k, arma::fill::zeros);
  arma::cube C_seq(C.n_rows, C.n_cols, k, arma::fill::zeros);
  arma::cube G_seq(C.n_rows, C.n_cols, k, arma::fill::zeros);
  arma::mat a_seq(m.size(), k, arma::fill::zeros);
  arma::mat m_seq(m.size(), k, arma::fill::zeros);
  arma::vec f_seq(k, arma::fill::zeros);
  arma::vec q_seq(k, arma::fill::zeros);
  arma::vec s_seq(k, arma::fill::zeros);
  arma::vec n_seq(k, arma::fill::zeros);

  double f = 0, q = 1, x_tf = 0;

  arma::vec a = m;
  arma::mat R = C;
  arma::vec F_t;

  for (int t = 0; t < k; ++t) {

    x_tf = xreg_tf(t);
    // Posterior at time t-1 to prior for time t
    evolve_prior(a, R,
                 m, C,
                 G, D,
                 s, n_parms,
                 ar_order,
                 tf_order,
                 i_ar,
                 i_tf,
                 x_tf);

    F_t = F.col(t);

    // y_t | D_{t-1}
    forecast_marginal(f, q,
                      F_t,
                      a, R, s,
                      ar_order);

    // (theta_t | D_{t})
    update_posterior(m, C, n, s,
                     y[t], f, q,
                     a, R, F_t, df_variance);

    // Saving the parameters
    f_seq(t) = f;
    q_seq(t) = q;
    m_seq.col(t) = m;
    C_seq.slice(t) = C;
    G_seq.slice(t) = G;
    a_seq.col(t) = a;
    R_seq.slice(t) = R;
    n_seq(t) = n;
    s_seq(t) = s;
  }

  return(Rcpp::List::create(
      Rcpp::Named("a")=a_seq,
      Rcpp::Named("R")=R_seq,
      Rcpp::Named("m")=m_seq,
      Rcpp::Named("C")=C_seq,
      Rcpp::Named("G")=G_seq,
      Rcpp::Named("n")=n_seq,
      Rcpp::Named("s")=s_seq,
      Rcpp::Named("f")=f_seq,
      Rcpp::Named("q")=q_seq));

}

// [[Rcpp::export]]
Rcpp::List backward_smoother_dlm(arma::mat F,
                                 arma::cube G_seq,
                                 arma::mat m_seq,
                                 arma::mat a_seq,
                                 arma::cube C_seq,
                                 arma::cube R_seq){

  int T_end = m_seq.n_cols;

  // Initialize a_t(0) and R_t(0)
  arma::vec ak = m_seq.col(T_end - 1);
  arma::mat Rk = C_seq.slice(T_end - 1);

  // Empty objects to keep the values
  arma::mat ak_seq(m_seq.n_rows, T_end, arma::fill::zeros);
  arma::cube Rk_seq(C_seq.n_rows, C_seq.n_cols, T_end, arma::fill::zeros);
  arma::vec fk_seq(T_end, arma::fill::zeros);
  arma::vec qk_seq(T_end, arma::fill::zeros);
  arma::mat B_t_k;

  // Starting values (we need to subtract -1, because it is the index in C++)
  ak_seq.col(T_end - 1) = ak;
  Rk_seq.slice(T_end - 1) = Rk;
  fk_seq(T_end - 1) = arma::as_scalar(F.col(0).t() * ak);
  qk_seq(T_end - 1) = arma::as_scalar(F.col(0).t() * Rk * F.col(0));

  // (theta_{t-k} | D_t) for k = 1,..., T - 1 and t = T.
  for (int k = 1; k < T_end; k++){

    // Using {inv_sympd} because is faster than {inv} and R is a covariance matrix (usually p.d.)

    // B_{t-k}
    B_t_k = C_seq.slice(T_end - k - 1) * G_seq.slice(T_end-k).t() * arma::inv_sympd(R_seq.slice(T_end-k), arma::inv_opts::allow_approx);

    // a_t(-k) and R_t(-k)
    ak = m_seq.col(T_end - k - 1) + B_t_k * (ak - a_seq.col(T_end - k));
    Rk = C_seq.slice(T_end - k - 1) + B_t_k * (Rk - R_seq.slice(T_end - k)) * B_t_k.t();

    // Saving the values
    ak_seq.col(T_end - k - 1) = ak;
    Rk_seq.slice(T_end - k - 1) = Rk;
    fk_seq(T_end - k - 1) = arma::as_scalar(F.col(T_end - k).t() * ak);
    qk_seq(T_end - k - 1) = arma::as_scalar(F.col(T_end - k).t() * Rk * F.col(T_end - k));

  }

  return(Rcpp::List::create(
      Rcpp::Named("ak")=ak_seq,
      Rcpp::Named("Rk")=Rk_seq,
      Rcpp::Named("fk")=fk_seq,
      Rcpp::Named("qk")=qk_seq));

}

