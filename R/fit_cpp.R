forward_filter_cpp.dlm <- function(model, y, a, R, n = 1, s = 1) {

  if (is.null(model[["X"]])) {
    out <- forward_filter_dlm(y = y, F = model[["FF"]], G = model[["GG"]],
                              D = model[["D"]], a = a, R = R, n = n,
                              s = s, df_variance = model[["df_variance"]])
  } else {
    cat("Not implemented yet.")
    out <- forward_filter_dlm_X(y = y, F = model[["FF"]], X = t(model[["X"]]),
                                G = model[["GG"]], D = model[["D"]],
                                a = a, R = R, n = n, s = s,
                                df_variance = model[["df_variance"]])
  }
  structure(out, class = "dlm.forward_filter")
}

backward_smooth_cpp.dlm <- function(model, filtering_parameters) {

  list()


}

fit_cpp.dlm <- function(model, y, a, R, n = 1, s = 1, smooth = TRUE) {

  filtered <- forward_filter_cpp.dlm(model = model, y = y,
                                       a = a, R = R, n = n, s = s)
  smoothed = NULL
  if (smooth) {
    smoothed <- list()
  }

  # Update the model object
  model$prior <- list()
  model$posterior <- list()

  # Return an object of class fit.dlm
  structure(list(model = model, filtered = filtered, smoothed = smoothed,
                 call = match.call()),
            class = "dlm.fit")

}
