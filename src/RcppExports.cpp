// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// forward_filter_dlm
Rcpp::List forward_filter_dlm(arma::vec y, arma::mat F, arma::mat G, arma::mat D, arma::vec m, arma::mat C, double n, double s, double df_variance, int n_parms, int ar_order, int tf_order, int fixed_tf_parm, Rcpp::IntegerVector i_ar, Rcpp::IntegerVector i_tf, Rcpp::NumericVector xreg_tf);
RcppExport SEXP _RBATS_forward_filter_dlm(SEXP ySEXP, SEXP FSEXP, SEXP GSEXP, SEXP DSEXP, SEXP mSEXP, SEXP CSEXP, SEXP nSEXP, SEXP sSEXP, SEXP df_varianceSEXP, SEXP n_parmsSEXP, SEXP ar_orderSEXP, SEXP tf_orderSEXP, SEXP fixed_tf_parmSEXP, SEXP i_arSEXP, SEXP i_tfSEXP, SEXP xreg_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type df_variance(df_varianceSEXP);
    Rcpp::traits::input_parameter< int >::type n_parms(n_parmsSEXP);
    Rcpp::traits::input_parameter< int >::type ar_order(ar_orderSEXP);
    Rcpp::traits::input_parameter< int >::type tf_order(tf_orderSEXP);
    Rcpp::traits::input_parameter< int >::type fixed_tf_parm(fixed_tf_parmSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type i_ar(i_arSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type i_tf(i_tfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xreg_tf(xreg_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_filter_dlm(y, F, G, D, m, C, n, s, df_variance, n_parms, ar_order, tf_order, fixed_tf_parm, i_ar, i_tf, xreg_tf));
    return rcpp_result_gen;
END_RCPP
}
// backward_smoother_dlm
Rcpp::List backward_smoother_dlm(arma::mat F, arma::cube G_seq, arma::mat m_seq, arma::mat a_seq, arma::cube C_seq, arma::cube R_seq);
RcppExport SEXP _RBATS_backward_smoother_dlm(SEXP FSEXP, SEXP G_seqSEXP, SEXP m_seqSEXP, SEXP a_seqSEXP, SEXP C_seqSEXP, SEXP R_seqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type G_seq(G_seqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type m_seq(m_seqSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a_seq(a_seqSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type C_seq(C_seqSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type R_seq(R_seqSEXP);
    rcpp_result_gen = Rcpp::wrap(backward_smoother_dlm(F, G_seq, m_seq, a_seq, C_seq, R_seq));
    return rcpp_result_gen;
END_RCPP
}
// forecast_dlm
Rcpp::List forecast_dlm(int horizon, arma::mat F, arma::mat G, arma::mat D, arma::vec m, arma::mat C, double s, int n_parms, int ar_order, int tf_order, int fixed_tf_parm, Rcpp::IntegerVector i_ar, Rcpp::IntegerVector i_tf, Rcpp::NumericVector xreg_tf);
RcppExport SEXP _RBATS_forecast_dlm(SEXP horizonSEXP, SEXP FSEXP, SEXP GSEXP, SEXP DSEXP, SEXP mSEXP, SEXP CSEXP, SEXP sSEXP, SEXP n_parmsSEXP, SEXP ar_orderSEXP, SEXP tf_orderSEXP, SEXP fixed_tf_parmSEXP, SEXP i_arSEXP, SEXP i_tfSEXP, SEXP xreg_tfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type n_parms(n_parmsSEXP);
    Rcpp::traits::input_parameter< int >::type ar_order(ar_orderSEXP);
    Rcpp::traits::input_parameter< int >::type tf_order(tf_orderSEXP);
    Rcpp::traits::input_parameter< int >::type fixed_tf_parm(fixed_tf_parmSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type i_ar(i_arSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type i_tf(i_tfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xreg_tf(xreg_tfSEXP);
    rcpp_result_gen = Rcpp::wrap(forecast_dlm(horizon, F, G, D, m, C, s, n_parms, ar_order, tf_order, fixed_tf_parm, i_ar, i_tf, xreg_tf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RBATS_forward_filter_dlm", (DL_FUNC) &_RBATS_forward_filter_dlm, 16},
    {"_RBATS_backward_smoother_dlm", (DL_FUNC) &_RBATS_backward_smoother_dlm, 6},
    {"_RBATS_forecast_dlm", (DL_FUNC) &_RBATS_forecast_dlm, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_RBATS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
