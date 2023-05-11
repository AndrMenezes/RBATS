#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List update_dlm(const double y,
                      const arma::vec F,
                      const arma::mat G,
                      const arma::mat D,
                      arma::vec a,
                      arma::mat R,
                      const double n,
                      const double s,
                      const double df_variance){

  arma::vec m = a;
  arma::mat C = R;
  double s_new = s;
  double n_new = n;
  double lpl = 0;

  // One-step ahead predictive moments
  double f = arma::as_scalar(F.t() * a);
  double q = s + arma::as_scalar(F.t() * R * F);
  // std::cout << F.t() << "\n";
  // std::cout << a << "\n";
  // std::cout << f << "\n";

  if (!std::isnan(y)) {

    double e = y - f;
    arma::vec A = (R * F) / q;
    double r = (df_variance * n + std::pow(e, 2) / q) / (df_variance * (n + 1));

    // Kalman filter update
    m = a + A * e;
    C = r * (R - q * A * A.t());
    s_new = r * s;
    n_new = n + 1;

    lpl = std::log(1/sqrt(q) * R::dt(e/sqrt(q), df_variance * n, FALSE));
  }

  // Evolution step: priors for time t + 1 from the posteriors m, C
  a = G * m;
  R = G * C * G.t();

  // Discount information
  R = D % R;


  return(Rcpp::List::create(
      Rcpp::Named("f")=f,
      Rcpp::Named("q")=q,
      Rcpp::Named("a")=a,
      Rcpp::Named("R")=R,
      Rcpp::Named("m")=m,
      Rcpp::Named("C")=C,
      Rcpp::Named("n")=df_variance * n_new,
      Rcpp::Named("s")=s_new,
      Rcpp::Named("lpl")=lpl));
}


Rcpp::List forward_filter_dlm_list(arma::vec const y,
                              arma::vec const F,
                              arma::mat const G,
                              arma::mat const D,
                              arma::vec a,
                              arma::mat R,
                              double n,
                              double s,
                              const double df_variance){

  const int k = y.size();

  Rcpp::List list_a = Rcpp::List(k + 1);
  Rcpp::List list_R = Rcpp::List(k + 1);
  Rcpp::List list_n = Rcpp::List(k + 1);
  Rcpp::List list_s = Rcpp::List(k + 1);
  Rcpp::List list_m = Rcpp::List(k);
  Rcpp::List list_C = Rcpp::List(k);
  Rcpp::List list_f = Rcpp::List(k);
  Rcpp::List list_q = Rcpp::List(k);
  list_a[0] = a;
  list_R[0] = R;
  list_n[0] = n;
  list_s[0] = s;
  double lpl = 0;

  for (int t = 0; t < k; ++t){
    Rcpp::List tmp = update_dlm(y[t], F, G, D, list_a[t], list_R[t],
                                list_n[t], list_s[t], df_variance);
    double curr_lpl = tmp["lpl"];
    lpl += curr_lpl;
    list_f(t) = tmp["f"];
    list_q(t) = tmp["q"];
    list_m(t) = tmp["m"];
    list_C(t) = tmp["C"];
    list_a(t + 1) = tmp["a"];
    list_R(t + 1) = tmp["R"];
    list_n(t + 1) = tmp["n"];
    list_s(t + 1) = tmp["s"];
  }

  return(Rcpp::List::create(
    Rcpp::Named("a")=list_a,
    Rcpp::Named("R")=list_R,
    Rcpp::Named("m")=list_m,
    Rcpp::Named("C")=list_C,
    Rcpp::Named("n")=list_n,
    Rcpp::Named("s")=list_s,
    Rcpp::Named("f")=list_f,
    Rcpp::Named("q")=list_q));
}

// [[Rcpp::export]]
Rcpp::List forward_filter_dlm(arma::vec const y,
                              arma::vec const F,
                              arma::mat const G,
                              arma::mat const D,
                              arma::vec a,
                              arma::mat R,
                              double n,
                              double s,
                              const double df_variance){

  const int k = y.size();

  arma::cube R_seq(R.n_rows, R.n_cols, k + 1, arma::fill::zeros);
  arma::cube C_seq(R.n_rows, R.n_cols, k, arma::fill::zeros);
  arma::mat a_seq(a.size(), k + 1, arma::fill::zeros);
  arma::mat m_seq(a.size(), k, arma::fill::zeros);
  arma::vec f_seq(k, arma::fill::zeros);
  arma::vec q_seq(k, arma::fill::zeros);
  arma::vec s_seq(k + 1, arma::fill::zeros);
  arma::vec n_seq(k + 1, arma::fill::zeros);

  a_seq.col(0) = a;
  R_seq.slice(0) = R;
  s_seq(0) = s;
  n_seq(0) = n;
  double lpl = 0;
  double curr_lpl = 0;

  for (int t = 0; t < k; ++t){
    Rcpp::List tmp = update_dlm(y[t], F, G, D, a_seq.col(t), R_seq.slice(t),
                                n_seq[t], s_seq[t], df_variance);

    curr_lpl = tmp["lpl"];
    lpl += curr_lpl;

    f_seq(t) = tmp["f"];
    q_seq(t) = tmp["q"];
    m_seq.col(t) = Rcpp::as<arma::vec>(tmp["m"]);
    C_seq.slice(t) = Rcpp::as<arma::mat>(tmp["C"]);
    a_seq.col(t + 1) = Rcpp::as<arma::vec>(tmp["a"]);
    R_seq.slice(t + 1) = Rcpp::as<arma::mat>(tmp["R"]);
    n_seq(t + 1) = tmp["n"];
    s_seq(t + 1) = tmp["s"];
  }

  return(Rcpp::List::create(
    Rcpp::Named("a")=a_seq,
    Rcpp::Named("R")=R_seq,
    Rcpp::Named("m")=m_seq,
    Rcpp::Named("C")=C_seq,
    Rcpp::Named("n")=n_seq,
    Rcpp::Named("s")=s_seq,
    Rcpp::Named("f")=f_seq,
    Rcpp::Named("q")=q_seq,
    Rcpp::Named("loglik")=lpl));

}

// [[Rcpp::export]]
Rcpp::List forward_filter_dlm_X(arma::vec const y,
                                arma::vec const F,
                                arma::mat const X,
                                arma::mat const G,
                                arma::mat const D,
                                arma::vec a,
                                arma::mat R,
                                double n,
                                double s,
                                const double df_variance){

  const int k = y.size();

  arma::cube R_seq(R.n_rows, R.n_cols, k + 1, arma::fill::zeros);
  arma::cube C_seq(R.n_rows, R.n_cols, k, arma::fill::zeros);
  arma::mat a_seq(a.size(), k + 1, arma::fill::zeros);
  arma::mat m_seq(a.size(), k, arma::fill::zeros);
  arma::vec f_seq(k, arma::fill::zeros);
  arma::vec q_seq(k, arma::fill::zeros);
  arma::vec s_seq(k + 1, arma::fill::zeros);
  arma::vec n_seq(k + 1, arma::fill::zeros);

  a_seq.col(0) = a;
  R_seq.slice(0) = R;
  s_seq(0) = s;
  n_seq(0) = n;
  double lpl = 0;
  double curr_lpl = 0;

  for (int t = 0; t < k; ++t){

    Rcpp::List tmp = update_dlm(y[t], arma::join_cols(F, X.col(t)), G, D,
                                a_seq.col(t), R_seq.slice(t),
                                n_seq[t], s_seq[t], df_variance);

    curr_lpl = tmp["lpl"];
    lpl += curr_lpl;

    f_seq(t) = tmp["f"];
    q_seq(t) = tmp["q"];
    m_seq.col(t) = Rcpp::as<arma::vec>(tmp["m"]);
    C_seq.slice(t) = Rcpp::as<arma::mat>(tmp["C"]);
    a_seq.col(t + 1) = Rcpp::as<arma::vec>(tmp["a"]);
    R_seq.slice(t + 1) = Rcpp::as<arma::mat>(tmp["R"]);
    n_seq(t + 1) = tmp["n"];
    s_seq(t + 1) = tmp["s"];
  }

  return(Rcpp::List::create(
    Rcpp::Named("a")=a_seq,
    Rcpp::Named("R")=R_seq,
    Rcpp::Named("m")=m_seq,
    Rcpp::Named("C")=C_seq,
    Rcpp::Named("n")=n_seq,
    Rcpp::Named("s")=s_seq,
    Rcpp::Named("f")=f_seq,
    Rcpp::Named("q")=q_seq,
    Rcpp::Named("loglik")=lpl));

}

Rcpp::List backward_smoother()



