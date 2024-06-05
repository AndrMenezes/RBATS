#' @name extract.dlm.fit
#'
#' @title Extract information of filtering and smoothing distribution from
#' \code{dlm.fit} object
#'
#' @description Tidy the results of the forward filter or backward smoother into
#' a \code{data.frame}.
#'
#' @return A \code{data.frame} with the results for each time.
#'
#' @param x an object of class \code{dlm.fit}.
#' @param horizon integer. The forecast horizon.
#' @param xreg,xreg_tf matrix. The corresponding matrix with the predictors for the regression and transfer function components.
#' @param prob_interval vector with the probability of the credibility interval.
#' @param state_parameters logical. Should returned the future prior moments of state
#' parameters? Default is \code{FALSE}.
#' @param discount logical. Whether to use discounting when forecasting. Default is \code{FALSE}.
#' @param ... currently not used.
#' @return If \code{state_parameters=FALSE} then return a \code{data.frame} with the response forecasts. Otherwise return a \code{list} with two components: \code{response} a data.frame with the response forecasts and \code{state} a data.frame with the state parameters forecasts.

#' @importFrom stats qt

#' @rdname forecast.dlm.fit
#' @export
forecast <- function(x, ...) {
  UseMethod("forecast", x)
}

#' @rdname forecast.dlm.fit
#' @export

forecast <- function(x, horizon, xreg = NULL, xreg_tf = NULL,
                     prob_interval = c(0.05, 0.20),
                     state_parameters = FALSE,
                     discount = FALSE) {

  model <- x$model

  last_time <- model$time
  time_index_future <- seq.int(last_time + 1L, by = 1L, length.out = horizon)

  # Discount forecast
  # set all elements equal to 1, for not discount!
  D <- matrix(1.0, nrow = nrow(model$D), ncol = ncol(model$D))
  if (discount) D <- model$D

  # Expand the regressor vector to transform into a matrix over time
  FF <- matrix(rep(model$FF, horizon), ncol = horizon)
  if (!is.null(xreg))
    FF[model$i_regressor, ] <- t(xreg)

  xreg_tf <- if (is.null(xreg_tf)) rep(0, horizon) else xreg_tf
  tmp <- forecast_dlm(horizon = horizon,
                      F = FF,
                      G = model$GG,
                      D = D,
                      m = as.matrix(model$posterior$m),
                      C = as.matrix(model$posterior$C),
                      s = model$posterior$s,
                      n_parms = model$n_parms,
                      ar_order = model$ar_order,
                      tf_order = model$tf_order,
                      fixed_tf_parm = model$fixed_tf_parm ,
                      i_ar = model$i_autoregressive - 1L,
                      i_tf = model$i_transfer_function - 1L,
                      xreg_tf = xreg_tf)

  # Last degree of freedom
  dof <- x$filtered$n[last_time, 1L]

  # Forecast the response
  data_response <- data.frame(t = time_index_future,
                              horizon = seq_len(horizon),
                              mean = tmp$f[, 1L],
                              variance = tmp$q[, 1L])

  for (p in prob_interval) {
    label <- as.character(100 * (1 - p))
    data_response[[paste0("ci_lower__", label)]] <- (
      data_response$mean + qt(p / 2, df = dof) * sqrt(data_response$variance))
    data_response[[paste0("ci_upper__", label)]] <- (
      data_response$mean + qt(
        1 - p / 2, df = dof) * sqrt(data_response$variance))
  }

  if (!state_parameters) {
    return(data_response)
  }

  # Forecast the state
  data_state <- data.frame(t = time_index_future)
  parms_names <- model$parameters_names
  list_to_append <- list()
  for (j in seq_along(parms_names)) {
    data_state$parameter <- parms_names[j]
    data_state$mean <- tmp$a[j, ]
    data_state$variance <- tmp$R[j, j, ]
    # Appending
    list_to_append[[j]] <- data_state
  }
  # Compute the sum of all seasonal effects
  if (!is.null(model$seasonal)) {
    i_seas <- model$i_seasonal
    F_seas <- FF[i_seas, ,drop = FALSE]
    m_k <- tmp$a[i_seas, ]
    C_k <- tmp$R[i_seas, i_seas, ]
    data_state$parameter <- "sum_seasonality"
    data_state$mean <- crossprod(F_seas, m_k)[1L, ]
    data_state$variance <- apply(C_k, 3L, function(Ct) crossprod(F_seas, Ct %*% F_seas)[1L, 1L])
    list_to_append[[j + 1L]] <- data_state
  }
  data_state <- do.call("rbind", list_to_append)

  # Compute the credibility intervals
  for (p in prob_interval) {
    label <- as.character(100 * (1 - p))
    data_state[[paste0("ci_lower__", label)]] <- (
      data_state$mean + qt(p / 2, df = dof) * sqrt(data_state$variance))
    data_state[[paste0("ci_upper__", label)]] <- (
      data_state$mean + qt(
        1 - p / 2, df = dof) * sqrt(data_state$variance))
  }
  list(response = data_response, state = data_state)

}
