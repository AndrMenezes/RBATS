#' @title Fitting Bayesian DLM and DGEGM
#'
#' @description Method \code{fit} applied to objects of class \code{dlm} or \code{dgem}
#' to perform the filtering and smoothing.
#'
#' @param object a model object of class \code{dlm} or \code{dgem}.
#' @param y vector. Observed value of time series.
#' @param smooth logical. Should smoothing be performed? Default is \code{TRUE}.
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
#'
#' Migon, H. M. and Gamerman, D. (1993). Generalized exponential growth models: A Bayesian approach,  \emph{Journal of Forecasting}, \bold{12}, 573--584

#' @rdname fit
#' @export
fit <- function(object, y, ...) {
  UseMethod("fit", object)
}

#' @rdname fit
#' @export
fit.dlm <- function(object, y, smooth = TRUE, prior_length = 10L, interval = TRUE,
                    level = 0.10, ...) {

  object <- .forward_filter_dlm(object = object, y = y, prior_length = prior_length,
                                interval = interval, level = level)

  if (smooth) object <- .backward_smoother_dlm(object, y = y, interval = interval,
                                               level = level)
  object$call <- match.call()
  structure(object, class = "dlm.fit")
}


.forward_filter_dlm <- function(object, y, prior_length = 10L, interval = TRUE,
                            level = 0.10) {

  object[["y"]] <- y
  object[["n"]] <- length(y)
  object[["time_index"]] <- seq_len(length(y))
  object[["prior_length"]] <- prior_length

  # Check if the prior was specified
  if (is.null(object[["prior"]][["a"]])) {
    tmp <- define_prior(object)
    object[["prior"]][[1L]][["a"]] <- tmp[["a"]]
    object[["prior"]][[1L]][["R"]] <- tmp[["R"]]
  }

  object[["posterior"]] <- list()
  object[["predictive"]] <- list()
  for (t in seq_along(y)) {
    # Posterior at t
    lt_tmp <- update_moments(object, time = t)

    # Posterior parameters for t
    object[["posterior"]][[t]] <- lt_tmp[["posterior"]]

    # Parameters for the one-step ahead predictive distribution
    object[["predictive"]][[t]] <- lt_tmp[["predictive"]]

    # Prior parameters for t
    object[["prior"]][[t + 1]] <- lt_tmp[["prior"]]
  }
  object[["filtered_parameters"]] <- TRUE

  # data.frame with the results
  object[["data_posterior_filter"]] <- .get_posterior(
    x = object, type = "filter", interval = interval, level = level)
  object[["data_predictive_filter"]] <- .get_predictive(
    x = object, type = "filter", interval = interval, level = level)

  object
}

.backward_smoother_dlm <- function(object, y, interval = TRUE, level = 0.10) {

  end_time <- length(y)
  # Start smooth parameters
  object[["smooth"]][[end_time]] <- list(
    a_T = object[["posterior"]][[end_time]][["m"]],
    R_T = object[["posterior"]][[end_time]][["C"]],
    f_T = object[["predictive"]][[end_time]][["f"]],
    q_T = object[["predictive"]][[end_time]][["q"]])
  object[["sT"]] <- object[["predictive"]][[end_time]][["s"]]
  object[["nT"]] <- object[["predictive"]][[end_time]][["n"]]

  # Get the smooth distribution parameters
  for (t in (end_time - 1):1) {
    # Compute smooth parameters
    lt_tmp <- retrospective_moments(object = object,  time = t)
    # Saving results
    object[["smooth"]][[t]] <- lt_tmp
  }
  object[["smooth_parameters"]] <- TRUE
  object[["data_posterior_smooth"]] <- .get_posterior(
    x = object, type = "smooth", interval = interval, level = level)
  object[["data_predictive_smooth"]] <- .get_predictive(
    x = object, type = "smooth", interval = interval, level = level)

  object
}


#' @rdname fit
#' @export
fit.dgegm <- function(object, y, interval = TRUE, level = 0.10, ...) {

  object[["y"]] <- y
  object[["n"]] <- length(y)
  object[["time_index"]] <- seq_len(length(y))

  # Check if the prior was specified
  if (is.null(object[["prior"]][[1L]][["a"]])) {
    tmp <- define_prior(object)
    object[["prior"]][[1L]][["a"]] <- tmp[["a"]]
    object[["prior"]][[1L]][["R"]] <- tmp[["R"]]
  }

  object[["posterior"]] <- list()
  object[["predictive"]] <- list()
  for (t in seq_along(y)) {
    lt_tmp <- update_moments(object = object, time = t)
    object[["posterior"]][[t]] <- lt_tmp[["posterior"]]
    object[["predictive"]][[t]] <- lt_tmp[["predictive"]]
    object[["prior"]][[t + 1]] <- lt_tmp[["prior"]]
  }

  # data.frame with the results
  object[["data_posterior_filter"]] <- .get_posterior(
    x = object, type = "filter", interval = interval, level = level)
  object[["data_predictive_filter"]] <- .get_predictive(
    x = object, type = "filter", interval = interval, level = level)

  structure(object, class = "dgegm.fit")
}

