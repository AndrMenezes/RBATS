#' @title Define prior moments for the state parameters.
#'
#' @description Helper functions applied on objects of class \code{dlm} and \code{dgegm},
#' to define a Dynamic Linear Models (DLM) or Dynamic Generalized Exponential Growth Models
#' (DGEGM), respectively. This function is especially useful if you do not know how to
#' specify a prior mean, \eqn{\mathbf{a}_0}, and covariance matrix,
#' \eqn{\mathbf{R}_0}, for the state vector.
#'
#' @param object an object of class \code{dlm} or \code{dgegm}.
#'
#' @details For \code{dlm} objects use linear model to define the prior
#' moments, the mean and covariance matrix, of model components.
#'
#' For \code{dgegm} objects the prior moments are:
#' for \eqn{\theta_1} is \eqn{g(y_1)}, for \eqn{\theta_2} is 0.001, and for \eqn{\theta_3}
#' is 0.95, where \eqn{g(\cdot)} is the Box-Cox link function and \eqn{y_1} is the first
#' observation of time series. The covariance
#'
#' @author André F. B. Menezes
#'
#' @importFrom stats lm poly ts vcov

#' @rdname define_prior
#' @export
define_prior <- function(object){
  UseMethod("define_prior", object)
}

#' @rdname define_prior
#' @export
define_prior.dlm <- function(object) {

  if (is.null(object[["prior_length"]]))
    object[["prior_length"]] <- 10L

  FF <- object[["FF"]]
  n <- object[["prior_length"]]
  y <- object[["y"]][seq_len(n)]
  t <- seq_len(n)
  X <- matrix(1, nrow = n)
  if (object[["polynomial_order"]] > 1L)
    X <- cbind(X, poly(t, degree = object[["polynomial_order"]] - 1))
  colnames(X) <- rownames(FF)[seq_len(object[["polynomial_order"]])]
  if (!is.null(object[["xreg"]]))
    X <- cbind(X, object[["xreg"]][t, , drop = FALSE])

  # Fit linear model and using coef and vcov as posterior at time t - 1
  fit_lm <- lm(y ~ X - 1)
  dlm_mean <- fit_lm[["coefficients"]]
  g <- max(2, floor(n / 2))
  dlm_vcov <- (g / (1 + g)) * vcov(fit_lm)
  dlm_vcov[-which(dlm_vcov == diag(dlm_vcov))] <- 0
  names(dlm_mean) <- colnames(X)
  colnames(dlm_vcov) <- rownames(dlm_vcov) <- colnames(X)

  if (object[["seasonal"]][["type"]] != "none") {
    FF_seas <- FF[-seq_len(object[["polynomial_order"]]), , drop = FALSE]
    dlm_mean_seas <- rep(0, nrow(FF_seas))
    dlm_vcov_seas <- diag(100, nrow = nrow(FF_seas), ncol = nrow(FF_seas))
    names(dlm_mean_seas) <- rownames(FF_seas)

    dlm_mean <- c(dlm_mean, dlm_mean_seas)
    dlm_vcov <- .bdiag(dlm_vcov, dlm_vcov_seas)
    colnames(dlm_vcov) <- rownames(dlm_vcov) <- c(colnames(X), rownames(FF_seas))
  }

  list(a = dlm_mean, R = dlm_vcov)
}

#' @rdname define_prior
#' @export
define_prior.dgegm <- function(object) {
  parms_names <- object[["parameters_names"]]
  y1 <- object[["y"]][1L]
  lambda <- object[["lambda"]]
  a11 <- if (lambda == 0) log(y1) else y1 ^ lambda
  a0 <- c(a11, 0.001, 0.95)
  R0 <- diag(100, ncol = 3, nrow = 3)
  colnames(R0) <- rownames(R0) <- parms_names
  names(a0) <- parms_names
  list(a = a0, R = R0)
}
