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
  stop("Not implemented yet")
}

#' @rdname forecast.dlm
#' @export
forecast_marginal.dlm <- function(object, horizon, xreg = NULL, ...) {
  stop("Not implemented yet")
}

#' @rdname forecast.dlm
#' @export
forecast_aR.dlm <- function(object, horizon, ...) {
  stop("Not implemented yet")
}
