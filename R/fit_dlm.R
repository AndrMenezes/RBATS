#' @title Forward filtering and Backward Smoothing for Bayesian Dynamic Linear Models
#'
#' @description Methods \code{forward_filter} and \code{backward_smoother}
#' for objects of class \code{dlm}. The implementation is written in
#' \code{C++}.
#'
#' @param model a model object of class \code{dlm}.
#' @param y vector. Observed value of time series.
#' @param a,R vector and matrix. prior state mean vector and covariance matrix.
#' @param n,s numeric. Prior sample size and mean for the variance, respectively.
#' For both parameter the default is 1.
#' @param smooth logical. Should smoothing be performed? Default is \code{TRUE}.
#' @param filtering_parameters list. Historical prior moments
#' (\code{a} and \code{R}) and historical posterior moments
#' (\code{m} and \code{C}).
#' @param ... currently not used.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' Prado, R.; West, M. Time Series Modeling, Computation, and Inference. CRC Press, 2010.

#' @rdname fit.dlm
#' @export
fit <- function(model, y, ...) {
  UseMethod("fit", model)
}

#' @rdname fit.dlm
#' @export
forward_filter <- function(model, y, ...) {
  UseMethod("forward_filter", model)
}

#' @rdname fit.dlm
#' @export
backward_smoother <- function(model, filtering_parameters, ...) {
  UseMethod("backward_smoother", model)
}

#' @rdname fit.dlm
#' @export
fit.dlm <- function(model, y, a, R, n = 1, s = 1, smooth = TRUE, ...) {

  filtered <- forward_filter.dlm(model = model, y = y,
                                     a = a, R = R, n = n, s = s)
  # Perform the smoothing
  smoothed <- NULL
  if (smooth) {
    smoothed <- backward_smoother.dlm(model = model,
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
  model["loglik"] <- filtered[["loglik"]]

  # Return an object of class fit.dlm
  structure(list(model = model, filtered = filtered, smoothed = smoothed,
                 call = match.call()),
            class = "dlm.fit")

}


#' @rdname fit.dlm
#' @export
forward_filter.dlm <- function(model, y, a, R, n = 1, s = 1, ...) {

  if (is.null(model[["X"]])) {
    # Check if it is to perform the monitor
    if (is.null(model[["monitor"]])) {
      out <- forward_filter_dlm(y = y, F = model[["FF"]], G = model[["GG"]],
                                D = model[["D"]], a = a, R = R, n = n,
                                s = s, df_variance = model[["df_variance"]])

    } else {
      # Check if the monitor is bilateral
      if (model[["monitor"]][["bilateral"]]) {
        out <- forward_filter_dlm_monitor_bilateral(
          y = y,
          F = model[["FF"]],
          G = model[["GG"]],
          D = model[["D"]],
          a = a,
          R = R,
          n = n,
          s = s,
          df_variance = model[["df_variance"]],
          bf_threshold = model[["monitor"]][["bf_threshold"]],
          location_shift = model[["monitor"]][["location_shift"]],
          scale_shift = model[["monitor"]][["scale_shift"]],
          exception_D = model[["monitor"]][["D"]],
          verbose = model[["monitor"]][["verbose"]],
          monitor_start = model[["monitor"]][["start_time"]])
      } else {
        out <- forward_filter_dlm_monitor(
          y = y,
          F = model[["FF"]],
          G = model[["GG"]],
          D = model[["D"]],
          a = a,
          R = R,
          n = n,
          s = s,
          df_variance = model[["df_variance"]],
          bf_threshold = model[["monitor"]][["bf_threshold"]],
          location_shift = model[["monitor"]][["location_shift"]],
          scale_shift = model[["monitor"]][["scale_shift"]],
          exception_D = model[["monitor"]][["D"]],
          verbose = model[["monitor"]][["verbose"]],
          monitor_start = model[["monitor"]][["start_time"]])
      }
    }
  } else {

    if (is.null(model[["monitor"]])) {
      out <- forward_filter_dlm_X(y = y, F = model[["FF"]], X = t(model[["X"]]),
                                  G = model[["GG"]], D = model[["D"]],
                                  a = a, R = R, n = n, s = s,
                                  df_variance = model[["df_variance"]])

    } else {
      # Check if monitor is bilateral
      if (model[["monitor"]][["bilateral"]]) {
        out <- forward_filter_dlm_monitor_bilateral_X(
          y = y,
          F = model[["FF"]],
          X = t(model[["X"]]),
          G = model[["GG"]],
          D = model[["D"]],
          a = a,
          R = R,
          n = n,
          s = s,
          df_variance = model[["df_variance"]],
          bf_threshold = model[["monitor"]][["bf_threshold"]],
          location_shift = model[["monitor"]][["location_shift"]],
          scale_shift = model[["monitor"]][["scale_shift"]],
          exception_D = model[["monitor"]][["D"]],
          verbose = model[["monitor"]][["verbose"]],
          monitor_start = model[["monitor"]][["start_time"]])
      }
      else {
        out <- forward_filter_dlm_monitor_X(
          y = y,
          F = model[["FF"]],
          X = t(model[["X"]]),
          G = model[["GG"]],
          D = model[["D"]],
          a = a,
          R = R,
          n = n,
          s = s,
          df_variance = model[["df_variance"]],
          bf_threshold = model[["monitor"]][["bf_threshold"]],
          location_shift = model[["monitor"]][["location_shift"]],
          scale_shift = model[["monitor"]][["scale_shift"]],
          exception_D = model[["monitor"]][["D"]],
          verbose = model[["monitor"]][["verbose"]],
          monitor_start = model[["monitor"]][["start_time"]])
      }

    }
  }

  structure(out, class = "dlm.forward_filter")
}

#' @rdname fit.dlm
#' @export
backward_smoother.dlm <- function(model, filtering_parameters, ...) {

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
