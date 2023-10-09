#' @title Forward filtering and Backward Smoothing for Bayesian Dynamic Linear Models
#'
#' @description Methods \code{forward_filter} and \code{backward_smoother}
#' for objects of class \code{dlm}. The implementation is written in
#' \code{C++}.
#'
#' @param model a model object of class \code{dlm}.
#' @param y vector. Observed value of time series.
#' @param m0,C0 vector and matrix. posterior mean vector and covariance matrix at time t-1.
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
fit.dlm <- function(model, y, m0, C0, n = 1, s = 1, smooth = TRUE, ...) {

  filtered <- forward_filter(model = model, y = y,
                             m0 = m0, C0 = C0, n = n, s = s)
  # Perform the smoothing
  smoothed <- NULL
  if (smooth) {
    smoothed <- backward_smoother(model = model,
                                  filtering_parameters = filtered)
  }

  # Update the model object
  t_end <- ncol(filtered[["m"]])
  model[["prior"]] <- list(a = filtered[["a"]][, t_end],
                           R = filtered[["R"]][, , t_end],
                           n = filtered[["n"]][t_end],
                           s = filtered[["s"]][t_end])
  model[["posterior"]] <- list(m = filtered[["m"]][, t_end],
                               C = filtered[["C"]][, , t_end],
                               n = filtered[["n"]][t_end],
                               s = filtered[["s"]][t_end])
  model[["m0"]]
  # Predictive log-likelihood
  s_q <- sqrt(filtered[["q"]][, 1L])
  e <- y - filtered[["f"]][, 1L]
  dof <- c(model[["df_variance"]] * n, filtered[["n"]][-length(y)])
  model["loglik"] <- sum(log(dt(x = e / s_q, df = dof) / s_q))
  model["time"] <- t_end

  # Return an object of class fit.dlm
  structure(list(model = model, y = y,
                 filtered = filtered, smoothed = smoothed,
                 n0 = n,
                 call = match.call()),
            class = "dlm.fit")

}


#' @rdname fit.dlm
#' @export
forward_filter.dlm <- function(model, y, m0, C0, n = 1, s = 1, ...) {

  # Expand the regressor vector to transform into a matrix over time
  FF <- matrix(rep(model[["FF"]], length(y)), ncol = length(y))
  if (!is.null(model[["xreg"]]))
    FF[model$i_regressor, ] <- t(model[["xreg"]])

  # Check if there is transfer function component, if not create an empty zero vector
  model[["xreg_tf"]] <- if (is.null(model[["xreg_tf"]])) rep(0, length(y)) else model[["xreg_tf"]]

  out <- forward_filter_dlm(y = y,
                            F = FF,
                            G = model[["GG"]],
                            D = model[["D"]],
                            m = m0,
                            C = C0,
                            n = n,
                            s = s,
                            df_variance = model[["df_variance"]],
                            n_parms = model[["n_parms"]],
                            ar_order = model[["ar_order"]],
                            tf_order = model[["tf_order"]],
                            i_ar = model[["i_autoregressive"]] - 1L,
                            i_tf = model[["i_transfer_function"]] - 1L,
                            xreg_tf = model[["xreg_tf"]]
                            )

  structure(out, class = "dlm.forward_filter")
}

#' @rdname fit.dlm
#' @export
backward_smoother.dlm <- function(model, filtering_parameters, ...) {

  # Expand the regressor vector to transform into a matrix over time
  FF <- matrix(rep(model[["FF"]], ncol(filtering_parameters[["a"]])),
               ncol = ncol(filtering_parameters[["a"]]))
  if (!is.null(model[["xreg"]]))
    FF[model[["i_regressor"]], ] <- t(model[["xreg"]])

  out <- backward_smoother_dlm(F = FF,
                               G_seq = filtering_parameters[["G"]],
                               m_seq = filtering_parameters[["m"]],
                               a_seq = filtering_parameters[["a"]],
                               C_seq = filtering_parameters[["C"]],
                               R_seq = filtering_parameters[["R"]])

  structure(out, class = "dlm.backward_smoother")
}
