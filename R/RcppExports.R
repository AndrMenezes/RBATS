# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

forward_filter_dlm <- function(y, F, G, D, m, C, n, s, df_variance, n_parms, ar_order, tf_order, fixed_tf_parm, i_ar, i_tf, xreg_tf) {
    .Call(`_RBATS_forward_filter_dlm`, y, F, G, D, m, C, n, s, df_variance, n_parms, ar_order, tf_order, fixed_tf_parm, i_ar, i_tf, xreg_tf)
}

backward_smoother_dlm <- function(F, G_seq, m_seq, a_seq, C_seq, R_seq) {
    .Call(`_RBATS_backward_smoother_dlm`, F, G_seq, m_seq, a_seq, C_seq, R_seq)
}

forecast_dlm <- function(horizon, F, G, D, m, C, s, n_parms, ar_order, tf_order, fixed_tf_parm, i_ar, i_tf, xreg_tf) {
    .Call(`_RBATS_forecast_dlm`, horizon, F, G, D, m, C, s, n_parms, ar_order, tf_order, fixed_tf_parm, i_ar, i_tf, xreg_tf)
}

