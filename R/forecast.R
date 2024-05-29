

forecast <- function(object, horizon, xreg = NULL, xreg_tf = NULL,
                     prob_interval = c(0.05, 0.20), state_parameters = FALSE) {

  model <- object$model

  last_time <- model$time
  time_index_future <- seq.int(last_time + 1L, by = 1L, length.out = horizon)

  # Discount factor
  D <- model$D
  # set all elements equal to 1, for not discount!
  # D <- matrix(1.0, nrow = nrow(D), ncol = ncol(D))

  # Expand the regressor vector to transform into a matrix over time
  FF <- matrix(rep(model$FF, horizon), ncol = horizon)
  if (!is.null(model[["xreg"]]))
    FF[model$i_regressor, ] <- t(model[["xreg"]])

  xreg_tf <- if (is.null(xreg_tf)) rep(0, horizon) else rep(0, h)
  tmp <- forecast_dlm(horizon = horizon,
                      F = FF,
                      G = model$GG,
                      D = D,
                      m = model$posterior$m,
                      C = model$posterior$C,
                      s = model$posterior$s,
                      n_parms = model$n_parms,
                      ar_order = model$ar_order,
                      tf_order = model$tf_order,
                      fixed_tf_parm = model$fixed_tf_parm ,
                      i_ar = model$i_autoregressive - 1L,
                      i_tf = model$i_transfer_function - 1L,
                      xreg_tf = xreg_tf)

  # Last degree of freedom
  dof <- object$filtered$n[last_time, 1L]

  # Forecast the response
  data_response <- data.frame(t = time_index_future,
                              horizon = seq_len(horizon),
                              mean = tmp$f[, 1L],
                              variance = tmp$q[, 1L])
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
    data_response[[paste0("ci_lower__", label)]] <- (
      data_response$mean + qt(p / 2, df = dof) * sqrt(data_response$variance))
    data_response[[paste0("ci_upper__", label)]] <- (
      data_response$mean + qt(
        1 - p / 2, df = dof) * sqrt(data_response$variance))
    # State parameters
    data_state[[paste0("ci_lower__", label)]] <- (
      data_state$mean + qt(p / 2, df = dof) * sqrt(data_state$variance))
    data_state[[paste0("ci_upper__", label)]] <- (
      data_state$mean + qt(
        1 - p / 2, df = dof) * sqrt(data_state$variance))
  }
  list(response = data_response, state = data_state)

}
