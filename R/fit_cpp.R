forward_filter_cpp.dlm <- function(model, y, a, R, n = 1, s = 1,
                                   monitor = FALSE,
                                   monitor_parameters = list(
                                     bf_threshold = 0.135,
                                     location_shift = 4,
                                     scale_shift = 1,
                                     discount_factors = list(
                                       level = 0.20,
                                       growth = 0.40,
                                       seasonal = 0.80,
                                       regressors = 0.80)
                                     )) {

  if (is.null(model[["X"]])) {
    out <- forward_filter_dlm(y = y, F = model[["FF"]], G = model[["GG"]],
                              D = model[["D"]], a = a, R = R, n = n,
                              s = s, df_variance = model[["df_variance"]])
  } else {
    out <- forward_filter_dlm_X(y = y, F = model[["FF"]], X = t(model[["X"]]),
                                G = model[["GG"]], D = model[["D"]],
                                a = a, R = R, n = n, s = s,
                                df_variance = model[["df_variance"]])
  }
  structure(out, class = "dlm.forward_filter")
}

backward_smoother_cpp.dlm <- function(model, filtering_parameters) {

  if (is.null(model[["X"]])) {
    out <- backward_smoother_dlm(F = model[["FF"]], G = model[["GG"]],
                                 m_seq = filtering_parameters[["m"]],
                                 a_seq = filtering_parameters[["a"]],
                                 C_seq = filtering_parameters[["C"]],
                                 R_seq = filtering_parameters[["R"]])
  } else {
    out <- backward_smoother_dlm_X(F = model[["FF"]], G = model[["GG"]],
                                   X = t(model[["X"]]),
                                   m_seq = filtering_parameters[["m"]],
                                   a_seq = filtering_parameters[["a"]],
                                   C_seq = filtering_parameters[["C"]],
                                   R_seq = filtering_parameters[["R"]])
  }
  structure(out, class = "dlm.backward_smoother")
}

fit_cpp.dlm <- function(model, y, a, R, n = 1, s = 1, smooth = TRUE) {

  filtered <- forward_filter_cpp.dlm(model = model, y = y,
                                       a = a, R = R, n = n, s = s)
  smoothed = NULL
  if (smooth) {
    smoothed <- backward_smoother_cpp.dlm(model = model,
                                          filtering_parameters = filtered)
  }

  # Update the model object
  t_end <- ncol(filtered[["m"]])
  model[["prior"]] <- list(a = filtered[["a"]][, t_end + 1],
                           R = filtered[["R"]][, , t_end + 1],
                           n = filtered[["n"]][t_end + 1],
                           s = filtered[["s"]][t_end + 1])
  model[["posterior"]] <- list(m = filtered[["m"]][, t_end],
                               C = filtered[["C"]][, , t_end],
                               n = filtered[["n"]][t_end],
                               s = filtered[["s"]][t_end])

  # Return an object of class fit.dlm
  structure(list(model = model, filtered = filtered, smoothed = smoothed,
                 call = match.call()),
            class = "dlm.fit")

}
