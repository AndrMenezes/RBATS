#' @title Forecast Bayesian Dynamic Linear Models.
#'
#' @description Method \code{forecast} applied to objects of class \code{dlm.fit}
#' to perform marginal forecast of predictive and/or state distributions.
#'
#' @param object fitted model object of class \code{dlm.fit}.
#' @param xreg matrix. Future values of possible regressors.
#' @param horizon integer. Forecast horizon.
#' @param interval logical. Should exact credible interval be returned?
#' @param level numeric. The probability level of the credible interval.
#' @param state_parameters logical. Should returned the future prior moments of state
#' parameters? Default is \code{FALSE}.
#' @param ... currently not used.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' @importFrom stats qt

#' @rdname forecast.dlm.fit
#' @export
forecast <- function(object, ...){
  UseMethod("forecast", object)
}

#' @rdname forecast.dlm.fit
#' @export
forecast.dlm.fit <- function(object, xreg = NULL, horizon, interval = TRUE,
                             level = 0.10, state_parameters = FALSE, ...) {

  max_time_index <- max(object[["time_index"]])
  time_index_future <- seq.int(max_time_index + 1, by = 1, length.out = horizon)

  if (!is.null(xreg)) {
    if (is.null(colnames(xreg)))
      colnames(xreg) <- colnames(object[["xreg"]])
  }

  lt_prior <- list()
  lt_predictive <- list()
  for (h in seq_len(horizon)) {
    tmp <- forecast_marginal.dlm.fit(
      object = object, xreg = xreg[h, ,drop = FALSE], horizon = h)
    lt_prior[[h]] <- tmp[["prior"]]
    lt_predictive[[h]] <- data.frame(
      t = time_index_future[h], horizon = h, mean = tmp[["predictive"]][["f"]],
      variance = tmp[["predictive"]][["q"]], df = tmp[["predictive"]][["n"]])
  }
  data_predictive <- do.call(rbind, lt_predictive)

  if (state_parameters) {
    parms_names <- rownames(object[["FF"]])
    lt_parms <- lapply(parms_names, function(parm) {
      tmp <- lapply(lt_prior, function(z) {
        diag_R <- diag(z[["R"]])
        c("mean" = unname(z[["a"]][which(rownames(z[["a"]]) == parm)]),
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
    class = "dlm.forecast"
  )

}

#' @rdname forecast.dlm.fit
#' @export
forecast_marginal.dlm.fit <- function(object, xreg = NULL, horizon) {

  FF <- object[["FF"]]
  if (!is.null(xreg))
    FF <- rbind(FF, t(xreg))

  # Observational variance
  last_period <- length(object[["predictive"]])
  st <- object[["predictive"]][[last_period]][["s"]]
  nt <- object[["predictive"]][[last_period]][["n"]]

  # Forecast prior mean (a) and variance (R)
  prior_parms <- forecast_aR.dlm.fit(object = object, horizon = horizon)
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

#' @rdname forecast.dlm.fit
#' @export
forecast_aR.dlm.fit <- function(object, horizon) {
  from <- length(object[["prior"]])
  G <- object[["GG"]]
  a <- object[["prior"]][[from]][["a"]]
  R <- object[["prior"]][[from]][["R"]]
  W <- object[["prior"]][[from]][["W"]]
  G_h <- .matrix_power(G, horizon)
  a_h <- G_h %*% a
  R_h <- tcrossprod(G_h %*% R, G_h) + W
  list(a_h = a_h, R_h = R_h)
}


