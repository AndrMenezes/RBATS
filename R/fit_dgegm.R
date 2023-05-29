#' @title Fitting Bayesian Dynamic Generalized Exponential Growth Model
#'
#' @description Method \code{fit} applied to models of class \code{dgegm}
#' to perform the filtering and smoothing.
#'
#' @param model a object of class \code{dgem}.
#' @param y vector. Observed value of time series.
#' @param smooth logical. Should smoothing be performed? Default is \code{TRUE}.
#' @param a0,R0 vector and matrix. prior state mean vector and covariance matrix.
#' Default is \code{NULL}, then use the first \code{prior_length} observation to compute.
#' @param n,s numeric. Prior sample size and mean for the variance, respectively.
#' For both parameter the default is 1.
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


#' @importFrom stats update
#' @rdname fit.dgegm
#' @export
fit.dgegm <- function(model, y, smooth = FALSE, a0 = NULL, R0 = NULL, n = 1, s = 1,
                      prior_length = 10L, interval = TRUE, level = 0.10, ...) {

  # Perform forward filter
  out <- forward_filter(
    model = model, y = y, a0 = a0, R0 = R0, n = n, s = s,
    prior_length = prior_length, interval = interval, level = level)

  # Perform backward smoother
  if (smooth)
    warning("Smooth is not implemented yet.")

  out$call <- match.call()
  structure(out, class = "dgegm.fit")
}

#' @rdname fit.dgegm
#' @export
forward_filter.dgegm <- function(model, y, a0 = NULL, R0 = NULL, n = 1, s = 1,
                                 prior_length = 10L, interval = TRUE, level = 0.10, ...) {


  # Check if the prior was specified
  if (is.null(a0) | is.null(R0)) {
    tmp <- define_prior(model, y = y)
    a0 <- tmp[["a"]]
    R0 <- tmp[["R"]]
  }
  colnames(R0) <- rownames(R0) <- names(a0) <- model[["parameters_names"]]
  list_prior <- list(list())
  list_prior[[1L]][["a"]] <- a0
  list_prior[[1L]][["R"]] <- R0
  list_prior[[1L]][["W"]] <- model[["D"]] * R0
  list_prior[[1L]][["n"]] <- n
  list_prior[[1L]][["s"]] <- s
  model[["prior"]] <- list_prior[[1L]]

  list_posterior <- list()
  list_predictive <- list()
  for (t in seq_along(y)) {

    # One-step ahead predictive parameters distribution
    list_predictive[[t]] <- forecast_marginal(model, horizon = 1)[["predictive"]]

    # Updating
    model <- update(model, y = y[t])

    # Posterior parameters for t
    list_posterior[[t]] <- model[["posterior"]]

    # Prior parameters for t + 1
    list_prior[[t + 1]] <- model[["prior"]]
  }


  # Organize the results in data frames
  data_predictive <- .tidy_predictive(
    list_predictive = list_predictive, interval = interval, level = level)
  data_posterior <- .tidy_posterior(
    list_posterior = list_posterior, interval = interval, level = level)

  list(model = model,
       y = y,
       data_predictive = data_predictive,
       data_posterior = data_posterior,
       state_parameters = list(prior = list_prior, posterior = list_posterior))
}
