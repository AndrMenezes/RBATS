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

//*****************************************************************************
// Forward filter (theta_t | D_{t-1})

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

// Same as above but with time-varying covariates
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

//*****************************************************************************
// Forward filter with automatic monitoring and adaptation

// [[Rcpp::export]]
double bayes_factor(double e, double mu, double tau){
  return(
    tau * std::exp((std::pow(e - mu, 2) - std::pow(e, 2)) / (2 * std::pow(tau,2)))
  );
  // return(
  //   std::exp(0.5 * (std::pow(mu, 2)) - e * mu)
  // );
}

// [[Rcpp::export]]
Rcpp::List forward_filter_dlm_monitor(arma::vec const y,
                                      arma::vec const F,
                                      arma::mat const G,
                                      arma::mat const D,
                                      arma::vec a,
                                      arma::mat R,
                                      double n,
                                      double s,
                                      const double df_variance,
                                      const double bf_threshold,
                                      const double location_shift,
                                      const double scale_shift,
                                      arma::mat const exception_D){

  const int k = y.size();
  const int p = a.size();

  // Empty objects to save the results
  arma::cube R_seq(R.n_rows, R.n_cols, k + 1, arma::fill::zeros);
  arma::cube C_seq(R.n_rows, R.n_cols, k, arma::fill::zeros);
  arma::mat a_seq(a.size(), k + 1, arma::fill::zeros);
  arma::mat m_seq(a.size(), k, arma::fill::zeros);
  arma::vec f_seq(k, arma::fill::zeros);
  arma::vec q_seq(k, arma::fill::zeros);
  arma::vec s_seq(k + 1, arma::fill::zeros);
  arma::vec n_seq(k + 1, arma::fill::zeros);
  arma::vec H_seq(k, arma::fill::zeros);
  arma::vec L_seq(k, arma::fill::zeros);
  arma::vec l_seq(k, arma::fill::zeros);

  //
  a_seq.col(0) = a;
  R_seq.slice(0) = R;
  s_seq(0) = s;
  n_seq(0) = n;
  double lpl = 0;
  double curr_lpl = 0;


  // Initialize the predictve
  double f = arma::as_scalar(F.t() * a);
  double q = s + arma::as_scalar(F.t() * R * F);
  double e = (y[0] - f) / std::sqrt(q);

  // Bayes factor and run-length
  double Ht = bayes_factor(e, location_shift, scale_shift);
  double Lt = Ht;
  int lt = 1;

  // Make the first update outside loop
  Rcpp::List tmp = update_dlm(y[0], F, G, D, a_seq.col(0), R_seq.slice(0),
                              n_seq[0], s_seq[0], df_variance);

  curr_lpl = tmp["lpl"];
  lpl += curr_lpl;

  f_seq(0) = tmp["f"];
  q_seq(0) = tmp["q"];
  H_seq(0) = Ht;
  L_seq(0) = Lt;
  l_seq(0) = lt;
  m_seq.col(0) = Rcpp::as<arma::vec>(tmp["m"]);
  C_seq.slice(0) = Rcpp::as<arma::mat>(tmp["C"]);
  a_seq.col(1) = Rcpp::as<arma::vec>(tmp["a"]);
  R_seq.slice(1) = Rcpp::as<arma::mat>(tmp["R"]);
  n_seq(1) = tmp["n"];
  s_seq(1) = tmp["s"];

  // Start the automatic monitor with intervention after time t=max(2*p, 6)
  int t1 = std::max(2 * p, 6);
  for (int t = 1; t < t1; ++t){
    // Compute the standardized error e_t ~ N(0, 1)
    f = arma::as_scalar(F.t() * a);
    q = s + arma::as_scalar(F.t() * R * F);
    e = (y[t] - f) / std::sqrt(q);

    // Bayes factor
    Ht = bayes_factor(e, location_shift, scale_shift);
    if (Lt < 1){
      lt = 1 + lt;
    } else {
      lt = 1;
    }
    Lt = Ht * std::min<double>(1, Lt);

    H_seq(t) = Ht;
    L_seq(t) = Lt;
    l_seq(t) = lt;

    // Bayes update
    Rcpp::List tmp = update_dlm(y[t], F, G, D, a_seq.col(t), R_seq.slice(t),
                                n_seq[t], s_seq[t], df_variance);

    // Update the loglike
    curr_lpl = tmp["lpl"];
    lpl += curr_lpl;

    a = Rcpp::as<arma::vec>(tmp["a"]);
    R = Rcpp::as<arma::mat>(tmp["R"]);

    // Save parameters
    f_seq(t) = f;
    q_seq(t) = q;
    m_seq.col(t) = Rcpp::as<arma::vec>(tmp["m"]);
    C_seq.slice(t) = Rcpp::as<arma::mat>(tmp["C"]);
    a_seq.col(t + 1) = a;
    R_seq.slice(t + 1) = Rcpp::as<arma::mat>(tmp["R"]);
    n_seq(t + 1) = tmp["n"];
    s_seq(t + 1) = tmp["s"];
  }

  // Start the automatic monitor with intervention after initial time
  double y_curr;
  for (int t = t1; t < k; ++t){

    y_curr = y[t];
    // Compute the standardized error e_t ~ N(0, 1)
    f = arma::as_scalar(F.t() * a);
    q = s + arma::as_scalar(F.t() * R * F);
    e = (y_curr - f) / std::sqrt(q);

    // Bayes factor
    Ht = bayes_factor(e, location_shift, scale_shift);
    if (Lt < 1){
      lt = 1 + lt;
    } else {
      lt = 1;
    }
    Lt = Ht * std::min<double>(1, Lt);

    H_seq(t) = Ht;
    L_seq(t) = Lt;
    l_seq(t) = lt;

    std::cout << "t: " << t << " f: " << f << " q: " << q;
    std::cout << " e: " << e << " H: " << Ht << " L: " << Lt <<  " l: " << lt << "\n";

    // Check for parametric change
    if ((Ht >= bf_threshold) & ((Lt < bf_threshold) | (lt > 2))){
      // Get back at time t - lt - 1 and increase the prior, R.
      int index = t - lt - 1;
      std::cout << "Parametric change at " << index << "\n";
      R_seq.slice(index) = R_seq.slice(index) % exception_D;
      for (int i = index; i < t; i++) {
        Rcpp::List tmp = update_dlm(y[i], F, G, D, a_seq.col(i), R_seq.slice(i),
                                    n_seq[i], s_seq[i], df_variance);
        f_seq(i) = tmp["f"];
        q_seq(i) = tmp["q"];
        m_seq.col(i) = Rcpp::as<arma::vec>(tmp["m"]);
        C_seq.slice(i) = Rcpp::as<arma::mat>(tmp["C"]);
        a_seq.col(i + 1) = Rcpp::as<arma::vec>(tmp["a"]);
        R_seq.slice(i + 1) = Rcpp::as<arma::mat>(tmp["R"]);
        n_seq(i + 1) = tmp["n"];
        s_seq(i + 1) = tmp["s"];
      }
      Lt = 1;
      lt = 0;
    }

    bool potential_outlier = FALSE;
    if (Ht < bf_threshold) {
      std::cout << "Potential outlier \n";
      potential_outlier = TRUE;
      y_curr = NA_REAL;
      Lt = 1;
      lt = 0;
    }

    Rcpp::List tmp = update_dlm(y_curr, F, G, D, a_seq.col(t), R_seq.slice(t),
                                n_seq[t], s_seq[t], df_variance);

    curr_lpl = tmp["lpl"];
    lpl += curr_lpl;

    a = Rcpp::as<arma::vec>(tmp["a"]);
    R = Rcpp::as<arma::mat>(tmp["R"]);

    if (potential_outlier){
      R = R % exception_D;
    }

    f_seq(t) = tmp["f"];
    q_seq(t) = tmp["q"];
    m_seq.col(t) = Rcpp::as<arma::vec>(tmp["m"]);
    C_seq.slice(t) = Rcpp::as<arma::mat>(tmp["C"]);
    a_seq.col(t + 1) = a;
    R_seq.slice(t + 1) = R;
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
      Rcpp::Named("H")=H_seq,
      Rcpp::Named("L")=L_seq,
      Rcpp::Named("l")=l_seq,
      Rcpp::Named("loglik")=lpl));

}



// Smoothing distribution based on the last observed time T_end
// [[Rcpp::export]]
Rcpp::List backward_smoother_dlm(arma::vec const F,
                                 arma::mat const G,
                                 arma::mat m_seq,
                                 arma::mat a_seq,
                                 arma::cube C_seq,
                                 arma::cube R_seq){

  double T_end = m_seq.n_cols;

  // Initialize a_t(0) and R_t(0)
  arma::vec ak = m_seq.col(T_end - 1);
  arma::mat Rk = C_seq.slice(T_end - 1);

  // Empty objects to keep the values
  arma::mat ak_seq(m_seq.n_rows, T_end, arma::fill::zeros);
  arma::cube Rk_seq(C_seq.n_rows, C_seq.n_cols, T_end, arma::fill::zeros);
  arma::vec fk_seq(T_end, arma::fill::zeros);
  arma::vec qk_seq(T_end, arma::fill::zeros);

  // Starting values (we need to subtract -1, because the index in C++)
  ak_seq.col(T_end - 1) = ak;
  Rk_seq.slice(T_end - 1) = Rk;
  fk_seq(T_end - 1) = arma::as_scalar(F.t() * ak);
  qk_seq(T_end - 1) = arma::as_scalar(F.t() * Rk * F);


  // (theta_{t-k} | D_t) for k = 1,..., T - 1 and t = T.
  for (int k = 1; k < T_end; k++){

    // Using {inv_sympd} because is faster than {inv} and R is a covariance matrix (usually p.d.)

    // B_{t-k}
    arma::mat B_t_k = C_seq.slice(T_end-k-1) * G.t() * arma::inv_sympd(R_seq.slice(T_end-k));

    // a_t(-k) and R_t(-k)
    ak = m_seq.col(T_end - k - 1) + B_t_k * (ak - a_seq.col(T_end - k));
    Rk = C_seq.slice(T_end - k - 1) + B_t_k * (Rk - R_seq.slice(T_end - k)) * B_t_k.t();

    // Saving the values
    ak_seq.col(T_end - k - 1) = ak;
    Rk_seq.slice(T_end - k - 1) = Rk;
    fk_seq(T_end - k - 1) = arma::as_scalar(F.t() * ak);
    qk_seq(T_end - k - 1) = arma::as_scalar(F.t() * Rk * F);

  }

  return(Rcpp::List::create(
      Rcpp::Named("ak")=ak_seq,
      Rcpp::Named("Rk")=Rk_seq,
      Rcpp::Named("fk")=fk_seq,
      Rcpp::Named("qk")=qk_seq));

}


// Smoothing distribution with covariates
// [[Rcpp::export]]
Rcpp::List backward_smoother_dlm_X(arma::vec const F,
                                   arma::mat const G,
                                   arma::mat const X,
                                   arma::mat m_seq,
                                   arma::mat a_seq,
                                   arma::cube C_seq,
                                   arma::cube R_seq){

  double T_end = m_seq.n_cols;

  // Initialize a_t(0) and R_t(0)
  arma::vec ak = m_seq.col(T_end - 1);
  arma::mat Rk = C_seq.slice(T_end - 1);

  // Empty objects to keep the values
  arma::mat ak_seq(m_seq.n_rows, T_end, arma::fill::zeros);
  arma::cube Rk_seq(C_seq.n_rows, C_seq.n_cols, T_end, arma::fill::zeros);
  arma::vec fk_seq(T_end, arma::fill::zeros);
  arma::vec qk_seq(T_end, arma::fill::zeros);

  // Define the regression vector for start values
  arma::vec Ft = arma::join_cols(F, X.col(T_end - 1));

  // Starting values (we need to subtract -1, because the index in C++)
  ak_seq.col(T_end - 1) = ak;
  Rk_seq.slice(T_end - 1) = Rk;
  fk_seq(T_end - 1) = arma::as_scalar(Ft.t() * ak);
  qk_seq(T_end - 1) = arma::as_scalar(Ft.t() * Rk * Ft);


  // (theta_{t-k} | D_t) for k = 1,..., T - 1 and t = T.
  for (int k = 1; k < T_end; k++){

    // Using {inv_sympd} because is faster than {inv} and R is a covariance matrix (usually p.d.)

    // B_{t-k}
    arma::mat B_t_k = C_seq.slice(T_end-k-1) * G.t() * arma::inv_sympd(R_seq.slice(T_end-k));

    // a_t(-k) and R_t(-k)
    ak = m_seq.col(T_end - k - 1) + B_t_k * (ak - a_seq.col(T_end - k));
    Rk = C_seq.slice(T_end - k - 1) + B_t_k * (Rk - R_seq.slice(T_end - k)) * B_t_k.t();

    Ft = arma::join_cols(F, X.col(T_end - k - 1));

    // Saving the values
    ak_seq.col(T_end - k - 1) = ak;
    Rk_seq.slice(T_end - k - 1) = Rk;
    fk_seq(T_end - k - 1) = arma::as_scalar(Ft.t() * ak);
    qk_seq(T_end - k - 1) = arma::as_scalar(Ft.t() * Rk * Ft);

  }

  return(Rcpp::List::create(
      Rcpp::Named("ak")=ak_seq,
      Rcpp::Named("Rk")=Rk_seq,
      Rcpp::Named("fk")=fk_seq,
      Rcpp::Named("qk")=qk_seq));

}


















