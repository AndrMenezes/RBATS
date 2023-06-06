#' @title Forecast Bayesian Dynamic Linear Models.
#'
#' @description Method \code{forecast} applied to objects of class \code{dlm.fit}
#' to perform marginal forecast of predictive and/or state distributions.
#'
#' @param object fitted model object of class \code{dlm.fit}.
#' @param xreg matrix. Future values of possible regressors.
#' @param horizon integer. Forecast horizon.
#' @param interval logical. Should exact credible interval be returned?
#' @param prob_interval numeric. The probability level of the credible interval.
#' @param state_parameters logical. Should returned the future prior moments of state
#' parameters? Default is \code{FALSE}.
#' @param ... currently not used.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' Migon, H. M. and Gamerman, D. (1993). Generalized exponential growth models: A Bayesian approach,  \emph{Journal of Forecasting}, \bold{12}, 573--584
#'
#' @importFrom stats qt

# General methods -------------------------------------------------------------------

#' @rdname forecast.dlm
#' @export
forecast <- function(object, horizon, ...){
  UseMethod("forecast", object)
}

#' @rdname forecast.dlm
#' @export
forecast_marginal <- function(object, horizon, ...){
  UseMethod("forecast_marginal", object)
}

#' @rdname forecast.dlm
#' @export
forecast_aR <- function(object, horizon, ...){
  UseMethod("forecast_aR", object)
}


# DLM -------------------------------------------------------------------------------

#' @rdname forecast.dlm
#' @export
forecast.dlm <- function(object, horizon, xreg = NULL, interval = TRUE,
                         prob_interval = 0.10, state_parameters = FALSE, ...) {

  max_time_index <- object[["time"]]
  time_index_future <- seq.int(max_time_index + 1, by = 1,
                               length.out = horizon)

  if (!is.null(xreg)) {
    if (is.null(colnames(xreg)))
      colnames(xreg) <- colnames(object[["xreg"]])
  }

  list_prior <- list()
  list_predictive <- list()
  for (h in seq_len(horizon)) {
    tmp <- forecast_marginal.dlm(
      object = object, xreg = xreg[h, ,drop = FALSE], horizon = h)
    list_prior[[h]] <- tmp[["prior"]]
    list_predictive[[h]] <- data.frame(
      t = time_index_future[h], horizon = h, mean = tmp[["predictive"]][["f"]],
      variance = tmp[["predictive"]][["q"]], df = tmp[["predictive"]][["n"]])
  }
  data_predictive <- do.call(rbind, list_predictive)

  if (state_parameters) {
    parms_names <- object$parameters_names
    lt_parms <- lapply(parms_names, function(parm) {
      tmp <- lapply(list_prior, function(z) {
        diag_R <- diag(z[["R"]])
        c("mean" = unname(z[["a"]][which(names(z[["a"]]) == parm)]),
          "variance" = unname(diag_R[names(diag_R) == parm]))
      })
      values <- do.call(rbind, tmp)
      out <- data.frame(t = time_index_future, horizon = seq_len(horizon))
      out$parameter <- parm
      out$mean <- values[, "mean"]
      out$variance <- values[, "variance"]
      out$df <- data_predictive$df
      out
    })
    data_state_parameters <- do.call(rbind, lt_parms)
  }

  if (interval) {
    data_predictive$ci_lower <- data_predictive$mean + (
      qt(prob_interval / 2, df = data_predictive$df) * sqrt(data_predictive$variance))
    data_predictive$ci_upper <- data_predictive$mean + (
      qt(1 - prob_interval / 2, df = data_predictive$df) * sqrt(data_predictive$variance))
    if (state_parameters) {
      data_state_parameters$ci_lower <- data_state_parameters$mean + (
        qt(prob_interval / 2, df = data_state_parameters$df) * sqrt(data_state_parameters$variance))
      data_state_parameters$ci_upper <- data_state_parameters$mean + (
        qt(1 - prob_interval / 2, df = data_state_parameters$df) * sqrt(data_state_parameters$variance))
    }
  }

  structure(
    list(
      predictive = data_predictive,
      state_parameters = if (state_parameters) data_state_parameters else NULL,
      call = match.call()),
    class = "dlm.forecast"
  )

}

#' @rdname forecast.dlm
#' @export
forecast_marginal.dlm <- function(object, horizon, xreg = NULL, ...) {

  FF <- object[["FF"]]
  if (!is.null(xreg))
    FF <- rbind(FF, t(xreg))

  # Observational variance
  st <- object[["prior"]][["s"]]
  nt <- object[["prior"]][["n"]]

  # Forecast prior mean (a) and variance (R)
  prior_parms <- forecast_aR(object = object, horizon = horizon)
  a_h <- prior_parms[["a_h"]]
  R_h <- prior_parms[["R_h"]]

  # Forecast marginal predictive distribution
  ft_h <- drop(crossprod(FF, a_h))
  k <- .variance_law(type = object[["variance_law"]][["type"]], mu = ft_h,
                     p = object[["variance_law"]][["power"]])
  qt_h <- drop(crossprod(FF, R_h %*% FF)) + st * k

  list(
    prior = list(a = a_h, R = R_h),
    predictive = list(f = ft_h, q = qt_h, n = nt)
  )
}

#' @rdname forecast.dlm
#' @export
forecast_aR.dlm <- function(object, horizon, ...) {

  a <- if (is.null(object[["posterior"]])) object[["prior"]][["a"]] else object[["posterior"]][["m"]]
  R <- if (is.null(object[["posterior"]])) object[["prior"]][["R"]] else object[["posterior"]][["C"]]
  G <- object[["GG"]]
  W <- (1 - object[["D"]]) * R
  G_h <- .matrix_power(G, horizon)
  a_h <- drop(G_h %*% a)
  R_h <- tcrossprod(G_h %*% R, G_h) + W
  list(a_h = a_h, R_h = R_h)
}


# DGEGM -----------------------------------------------------------------------------
#' @rdname forecast.dgem
#' @export
forecast.dgegm <- function(object, horizon, interval = TRUE, level = 0.10,
                           state_parameters = FALSE, ...) {

  max_time_index <- object[["time"]]
  time_index_future <- seq.int(max_time_index + 1, by = 1, length.out = horizon)

  list_prior <- list()
  list_predictive <- list()
  for (h in seq_len(horizon)) {
    tmp <- forecast_marginal.dgegm(object = object, horizon = h)
    list_prior[[h]] <- tmp[["prior"]]
    list_predictive[[h]] <- data.frame(
      t = time_index_future[h], horizon = h,
      mean = tmp[["predictive"]][["f"]],
      variance = tmp[["predictive"]][["q"]],
      df = tmp[["predictive"]][["n"]])
  }
  data_predictive <- do.call(rbind, list_predictive)

  if (state_parameters) {
    parms_names <- object[["parameters_names"]]
    lt_parms <- lapply(parms_names, function(parm) {
      tmp <- lapply(list_prior, function(z) {
        diag_R <- diag(z[["R"]])
        c("mean" = unname(z[["a"]][which(names(z[["a"]]) == parm)]),
          "variance" = unname(diag_R[names(diag_R) == parm]))
      })
      values <- do.call(rbind, tmp)
      out <- data.frame(t = time_index_future, horizon = seq_len(horizon))
      out$parameter <- parm
      out$mean <- values[, "mean"]
      out$variance <- values[, "variance"]
      out$df <- data_predictive$df
      out
    })
    data_state_parameters <- do.call(rbind, lt_parms)
  }

  if (interval) {
    data_predictive$ci_lower <- data_predictive$mean + (
      qt(level / 2, df = data_predictive$df) * sqrt(data_predictive$variance))
    data_predictive$ci_upper <- data_predictive$mean + (
      qt(1 - level / 2, df = data_predictive$df) * sqrt(data_predictive$variance))
    if (state_parameters) {
      data_state_parameters$ci_lower <- data_state_parameters$mean + (
        qt(level / 2, df = data_state_parameters$df) * sqrt(data_state_parameters$variance))
      data_state_parameters$ci_upper <- data_state_parameters$mean + (
        qt(1 - level / 2, df = data_state_parameters$df) * sqrt(data_state_parameters$variance))
    }
  }

  structure(
    list(
      predictive = data_predictive,
      state_parameters = if (state_parameters) data_state_parameters else NULL,
      call = match.call()),
    class = "dgegm.forecast"
  )

}


#' @rdname forecast.dgegm
#' @export
forecast_marginal.dgegm <- function(object, horizon, ...) {

  FF <- object[["FF"]]
  FF_prime <- object[["FF_prime"]]
  g <- object[["g"]]
  GG <- object[["GG"]]
  lambda <- object[["lambda"]]

  # Observational variance
  st <- object[["prior"]][["s"]]
  nt <- object[["prior"]][["n"]]

  # Forecast prior mean (a) and variance (R)
  prior_parms <- forecast_aR.dgegm(object = object, horizon = horizon)
  a_h <- prior_parms[["a_h"]]
  R_h <- prior_parms[["R_h"]]

  # Forecast marginal predictive distribution
  ft_h <- FF(theta = a_h, lambda = lambda)
  FF_prime_a <- FF_prime(theta = a_h, lambda = lambda)
  k <- .variance_law(mu = ft_h, type = object[["variance_law"]][["type"]],
                     p = object[["variance_law"]][["p"]])

  qt_h <- drop(crossprod(FF_prime_a, R_h %*% FF_prime_a)) + st * k

  list(
    prior = list(a = a_h, R = R_h),
    predictive = list(f = ft_h, q = qt_h, n = nt)
  )
}

#' @rdname forecast.dgegm
#' @export
forecast_aR.dgegm <- function(object, horizon, ...) {

  g <- object[["g"]]
  GG <- object[["GG"]]

  a_h <- if (is.null(object[["posterior"]])) object[["prior"]][["a"]] else object[["posterior"]][["m"]]
  R_h <- if (is.null(object[["posterior"]])) object[["prior"]][["R"]] else object[["posterior"]][["C"]]
  W <- object[["prior"]][["W"]]

  i <- 1L
  while (i <= horizon) {
    GG_h <- GG(a_h)
    a_h <- g(a_h)
    R_h <- tcrossprod(GG_h %*% R_h, GG_h) + W
    i <- i + 1L
  }
  list(a_h = a_h, R_h = R_h)
}




