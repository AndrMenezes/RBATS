#' @title Fitting Bayesian Dynamic Linear Model
#'
#' @description Method \code{fit} applied to objects of class \code{dlm}
#' to perform the filtering and smoothing.
#'
#' @param object a model object of class \code{dlm}.
#' @param y vector. Observed value of time series.
#' @param smooth logical. Should smoothing be performed? Default is \code{TRUE}.
#' @param a0,R0 vector and matrix. prior state mean vector and covariance matrix.
#' Default is \code{NULL}, then use the first \code{prior_length} observation to compute.
#' @param n,s numeric. Prior sample size and mean for the variance, respectively.
#' For both parameter the default is 1.
#' @param state_parameters list. Historical prior moments (\code{a} and \code{R}) and
#' historical posterior moments (\code{m} and \code{C})
#' @param prior_length integer. Number of observations to use to obtain the prior moments.
#' Default is 10.
#' @param interval logical. Should exact credible interval be returned? Default is \code{TRUE}.
#' @param level numeric. The probability level of the credible interval. Default is \code{0.10}.
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
fit <- function(object, y, ...) {
  UseMethod("fit", object)
}

#' @rdname fit.dlm
#' @export
forward_filter <- function(object, y, ...) {
  UseMethod("forward_filter", object)
}

#' @rdname fit.dlm
#' @export
backward_smoother <- function(object, state_parameters, ...) {
  UseMethod("forward_filter", object)
}


#' @rdname fit.dlm
#' @export
fit.dlm <- function(object, y, smooth = TRUE, a0 = NULL, R0 = NULL, n = 1, s = 1,
                    prior_length = 10L, interval = TRUE, level = 0.10, ...) {

  # Perform forward filter
  out <- forward_filter.dlm(
    object = object, y = y, a0 = a0, R0 = R0, n = n, s = s,
    prior_length = prior_length, interval = interval, level = level)

  # Perform backward smoother
  if (smooth) {
    out_smooth <- backward_smoother.dlm(
      object = object, state_parameters = out[["state_parameters"]],
      interval = interval, level = level)
    out[["smooth_parameters"]] <- out_smooth[["smooth_parameters"]]
    out[["data_posterior_smooth"]] <- out_smooth[["data_posterior"]]
    out[["data_predictive_smooth"]] <- out_smooth[["data_predictive"]]
  }

  out$call <- match.call()
  structure(out, class = "dlm.fit")
}

#' @rdname fit.dlm
#' @export
forward_filter.dlm <- function(object, y, a0 = NULL, R0 = NULL, n = 1, s = 1,
                               prior_length = 10L, interval = TRUE, level = 0.10, ...) {


  # Check if the prior was specified
  if (is.null(a0) | is.null(R0)) {
    object[["prior_length"]] <- prior_length
    tmp <- define_prior(object, y = y[seq_len(prior_length)])
    a0 <- tmp[["a"]]
    R0 <- tmp[["R"]]
  }
  colnames(R0) <- rownames(R0) <- names(a0) <- object[["parameters_names"]]
  list_prior <- list(list())
  list_prior[[1L]][["a"]] <- a0
  list_prior[[1L]][["R"]] <- R0
  list_prior[[1L]][["W"]] <- object[["D"]] * R0
  list_prior[[1L]][["n"]] <- n
  list_prior[[1L]][["s"]] <- s
  object[["prior"]] <- list_prior[[1L]]

  list_posterior <- list()
  list_predictive <- list()
  for (t in seq_along(y)) {

    # One-step ahead predictive parameters distribution
    list_predictive[[t]] <- forecast_marginal(
      object, xreg = object[["xreg"]][t, , drop = FALSE], horizon = 1)[["predictive"]]

    # Updating
    object <- update_moments(object, y = y[t], X = object[["xreg"]][t, , drop = FALSE])

    # Posterior parameters for t
    list_posterior[[t]] <- object[["posterior"]]

    # Prior parameters for t + 1
    list_prior[[t + 1]] <- object[["prior"]]
  }

  # Organize the results in data frames
  data_predictive <- .tidy_predictive(
    list_predictive = list_predictive, interval = interval, level = level)
  data_posterior <- .tidy_posterior(
    list_posterior = list_posterior, interval = interval, level = level)

  list(model = object,
       y = y,
       data_predictive = data_predictive,
       data_posterior = data_posterior,
       state_parameters = list(prior = list_prior, posterior = list_posterior))
}

#' @rdname fit.dlm
#' @export
backward_smoother.dlm <- function(object, state_parameters, interval = TRUE,
                                  level = 0.10, ...) {

  end_time <- length(state_parameters[["posterior"]])
  # Model components
  GG <- object[["GG"]]
  FF <- object[["FF"]]

  # Start smooth parameters
  a_T <- state_parameters[["posterior"]][[end_time]][["m"]]
  R_T <- state_parameters[["posterior"]][[end_time]][["C"]]
  s_T <- state_parameters[["posterior"]][[end_time]][["s"]]
  n_T <- state_parameters[["posterior"]][[end_time]][["n"]]

  # Start predictive smooth
  if (!is.null(object[["xreg"]]))
    FF <- rbind(FF, t(object[["xreg"]][end_time, , drop = FALSE]))
  f_T <- drop(crossprod(FF, a_T))
  q_T <- drop(crossprod(FF, R_T %*% FF))

  # To save results
  list_smooth <- list()
  list_smooth[[end_time]] <- list(a = a_T, R = R_T, f = f_T, q = q_T, s = s_T)
  list_smooth[["nT"]] <- n_T
  list_smooth[["sT"]] <- s_T

  # Perform smooth on parameters
  for (t in (end_time - 1):1) {

    FF <- object[["FF"]]
    if (!is.null(object[["xreg"]]))
      FF <- rbind(FF, t(object[["xreg"]][t, , drop = FALSE]))

    # Past smooth parameters: a_T(t - T + 1) and R_T(t - T + 1)
    a_T_t_1 <- list_smooth[[t + 1]][["a"]]
    R_T_t_1 <- list_smooth[[t + 1]][["R"]]

    # Posterior at t and prior at t + 1
    m_t <- state_parameters[["posterior"]][[t]][["m"]]
    C_t <- state_parameters[["posterior"]][[t]][["C"]]
    s_t <- state_parameters[["posterior"]][[t]][["s"]]
    a_t_1 <- state_parameters[["prior"]][[t + 1]][["a"]]
    R_t_1 <- state_parameters[["prior"]][[t + 1]][["R"]]
    
    # inv_R_t_1 <- solve(R_t_1)
    inv_R_t_1 <- chol2inv(chol(R_t_1))
    B_t <- tcrossprod(C_t, GG) %*% inv_R_t_1


    # Retrospective prior: a_T(t - T) and R_T(t - T)
    a_T <- drop(m_t - B_t %*% (a_t_1 - a_T_t_1))
    R_T <- C_t - tcrossprod(B_t %*% (R_t_1 - R_T_t_1), B_t)

    # f_t(-k) e q_t(-k)
    f_T <- drop(crossprod(FF, a_T))
    q_T <- drop(crossprod(FF, R_T %*% FF))

    # Saving results
    list_smooth[[t]][["a"]] <- a_T
    list_smooth[[t]][["R"]] <- R_T
    list_smooth[[t]][["f"]] <- f_T
    list_smooth[[t]][["q"]] <- q_T
    list_smooth[[t]][["s"]] <- s_t

  }

  # Tidy data in data frame
  data_posterior <- .tidy_smooth(
    list_smooth = list_smooth, what = "posterior", interval = interval, level = level)
  data_predictive <- .tidy_smooth(
    list_smooth = list_smooth, what = "predictive", interval = interval, level = level)

  list(
    smooth_parameters = list_smooth,
    data_posterior = data_posterior,
    data_predictive = data_predictive
  )

}


# DGEGM -----------------------------------------------------------------------------
#' @title Fitting Bayesian Dynamic Generalized Exponential Growth Model
#'
#' @description Method \code{fit} applied to objects of class \code{dgem}
#' to perform the filtering and smoothing.
#'
#' @param object a model object of class \code{dgem}.
#' @param y vector. Observed value of time series.
#' @param smooth logical. Should smoothing be performed? Default is \code{TRUE}.
#' @param a0,R0 vector and matrix. prior state mean vector and covariance matrix.
#' Default is \code{NULL}, then use the first \code{prior_length} observation to compute.
#' @param n,s numeric. Prior sample size and mean for the variance, respectively.
#' For both parameter the default is 1.
#' @param state_parameters list. Historical prior moments (\code{a} and \code{R}) and
#' historical posterior moments (\code{m} and \code{C})
#' @param prior_length integer. Number of observations to use to obtain the prior moments.
#' Default is 10.
#' @param interval logical. Should exact credible interval be returned? Default is \code{TRUE}.
#' @param level numeric. The probability level of the credible interval. Default is \code{0.10}.
#' @param ... currently not used.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' Migon, H. M. and Gamerman, D. (1993). Generalized exponential growth models: A Bayesian approach,  \emph{Journal of Forecasting}, \bold{12}, 573--584


#' @rdname fit.dgegm
#' @export
fit.dgegm <- function(object, y, smooth = FALSE, a0 = NULL, R0 = NULL, n = 1, s = 1,
                      prior_length = 10L, interval = TRUE, level = 0.10, ...) {

  # Perform forward filter
  out <- forward_filter(
    object = object, y = y, a0 = a0, R0 = R0, n = n, s = s,
    prior_length = prior_length, interval = interval, level = level)

  # Perform backward smoother
  if (smooth)
    warning("Smooth is not implemented yet.")

  out$call <- match.call()
  structure(out, class = "dgegm.fit")
}

#' @rdname fit.dgegm
#' @export
forward_filter.dgegm <- function(object, y, a0 = NULL, R0 = NULL, n = 1, s = 1,
                                 prior_length = 10L, interval = TRUE, level = 0.10, ...) {


  # Check if the prior was specified
  if (is.null(a0) | is.null(R0)) {
    tmp <- define_prior(object, y = y)
    a0 <- tmp[["a"]]
    R0 <- tmp[["R"]]
  }
  colnames(R0) <- rownames(R0) <- names(a0) <- object[["parameters_names"]]
  list_prior <- list(list())
  list_prior[[1L]][["a"]] <- a0
  list_prior[[1L]][["R"]] <- R0
  list_prior[[1L]][["W"]] <- object[["D"]] * R0
  list_prior[[1L]][["n"]] <- n
  list_prior[[1L]][["s"]] <- s
  object[["prior"]] <- list_prior[[1L]]

  list_posterior <- list()
  list_predictive <- list()
  for (t in seq_along(y)) {

    # One-step ahead predictive parameters distribution
    list_predictive[[t]] <- forecast_marginal(object, horizon = 1)[["predictive"]]

    # Updating
    object <- update_moments(object, y = y[t])

    # Posterior parameters for t
    list_posterior[[t]] <- object[["posterior"]]

    # Prior parameters for t + 1
    list_prior[[t + 1]] <- object[["prior"]]
  }


  # Organize the results in data frames
  data_predictive <- .tidy_predictive(
    list_predictive = list_predictive, interval = interval, level = level)
  data_posterior <- .tidy_posterior(
    list_posterior = list_posterior, interval = interval, level = level)

  list(model = object,
       y = y,
       data_predictive = data_predictive,
       data_posterior = data_posterior,
       state_parameters = list(prior = list_prior, posterior = list_posterior))
}
