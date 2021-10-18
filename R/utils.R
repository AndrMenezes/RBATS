.variance_law <- function(type, mu, p) {
  switch (type,
          "identity" = 1,
          "poisson" = mu,
          "binomial" = mu * (1 - mu),
          "power" = mu ^ p
  )
}

.diag_one <- function(x) {
  m <- diag(x, nrow = length(x), ncol = length(x))
  m[which(m != diag(m))] <- 1
  m
}

.matrix_power <- function(x, k) {
  x_k <- x
  i <- 2L
  while (i <= k) {
    x_k <- x_k %*% x
    i <- i + 1L
  }
  x_k
}

.inverse_link <- function(lambda) {
  if (lambda == 0) return(function(mu, lambda) exp(mu))
  else return(function(mu, lambda) mu^(1/lambda))
}
.inverse_link_diff <- function(lambda) {
  if (lambda == 0) return(function(mu, lambda) exp(mu))
  else return(function(mu, lambda) mu^(1/lambda - 1) / lambda)
}


# Organize the predictive distribution parameters in data.frame
#' @importFrom stats qt
.get_predictive <- function(x, type = c("filter", "smooth"), interval = TRUE,
                            level = 0.10) {
  type <- match.arg(type)
  y <- x[["y"]]
  n <- length(y)
  n_prior <- if (is.null(x[["prior_length"]])) 0 else x[["prior_length"]]
  data_output <- data.frame(
    t = x[["time_index"]], y = y,
    prior = rep(c(TRUE, FALSE), times = c(n_prior, n - n_prior)))
  if (type == "filter") {
    data_output$mean <- sapply(x[["predictive"]], function(z) unname(z[["f"]]))
    data_output$variance <- sapply(x[["predictive"]], function(z) unname(z[["q"]]))
    data_output$df <- sapply(x[["predictive"]], function(z) unname(z[["n"]]))
  }
  if (type == "smooth") {
    data_output$mean <- sapply(x[["smooth"]], function(z) unname(z[["f_T"]]))
    data_output$variance <- sapply(x[["smooth"]], function(z) unname(z[["q_T"]]))
    data_output$df <- x[["nT"]]
  }

  if (interval) {
    data_output$ci_lower <- data_output$mean + qt(level / 2, df = data_output$df) * sqrt(data_output$variance)
    data_output$ci_upper <- data_output$mean + qt(1 - level / 2, df = data_output$df) * sqrt(data_output$variance)
  }
  data_output$type <- type
  data_output
}

# Organize the posterior distribution parameters in data.frame
#' @importFrom stats qt
.get_posterior <- function(x, type = c("filter", "smooth"), interval = TRUE, level = 0.10) {

  n <- length(x[["y"]])
  n_prior <- if (is.null(x[["prior_length"]])) 0 else x[["prior_length"]]
  parms_names <- x[["parameters_names"]]

  # Degrees of freedom
  df <- sapply(x[["posterior"]], function(z) unname(z[["n"]]))

  if (type == "filter") {
    lt_parms <- lapply(parms_names, function(parm) {
      tmp <- lapply(x[["posterior"]], function(z) {
        diag_C <- diag(z[["C"]])
        c("mean" = unname(z[["m"]][which(names(z[["m"]]) == parm)]),
          "variance" = unname(diag_C[names(diag_C) == parm]))
      })
      values <- do.call(rbind, tmp)
      out <- data.frame(t = x[["time_index"]],
                        prior = rep(c(TRUE, FALSE), times = c(n_prior, n - n_prior)))
      out$parameter <- parm
      out$mean <- values[, "mean"]
      out$variance <- values[, "variance"]
      out$df <- df
      out
    })
  }
  if (type == "smooth") {
    lt_parms <- lapply(parms_names, function(parm) {
      tmp <- lapply(x[["smooth"]], function(z) {
        diag_R <- diag(z[["R_T"]])
        c("mean" = unname(z[["a_T"]][which(names(z[["a_T"]]) == parm)]),
          "variance" = unname(diag_R[names(diag_R) == parm]))
      })
      out <- data.frame(t = x[["time_index"]],
                        prior = rep(c(TRUE, FALSE), times = c(n_prior, n - n_prior)))
      out$parameter <- parm
      values <- do.call(rbind, tmp)
      out$mean <- values[, "mean"]
      out$variance <- values[, "variance"]
      out$df <- x[["nT"]]
      out
    })
  }
  data_output <- do.call(rbind, lt_parms)

  if (interval) {
    data_output$ci_lower <- data_output$mean + qt(level / 2, df = data_output$df) * sqrt(data_output$variance)
    data_output$ci_upper <- data_output$mean + qt(1 - level / 2, df = data_output$df) * sqrt(data_output$variance)
  }
  data_output$type <- type
  data_output
}

