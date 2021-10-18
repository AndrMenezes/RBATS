#' @title Update the moments of state parameters
#'
#' @description Methods to update the moments of state parameters and compute the
#' one-step ahead predictive distributions for objects of class  \code{dlm} and
#' \code{dgegm}.
#'
#' @param object an object of \code{dlm} or \code{dgegm}.
#' @param time integer. The time index to update the moments.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' Migon, H. M. and Gamerman, D. (1993). Generalized exponential growth models: A Bayesian approach,  \emph{Journal of Forecasting}, \bold{12}, 573--584

#' @rdname update_moments
#' @export
update_moments <- function(object, time){
  UseMethod("update_moments", object)
}

#' @rdname update_moments
#' @export
update_moments.dlm <- function(object, time) {

  # Model components
  FF <- object[["FF"]]
  GG <- object[["GG"]]
  D <- object[["D"]]
  df_variance <- object[["df_variance"]]

  # Get the observation time series and possible covariates
  y <- object[["y"]][time]
  if (!is.null(object[["xreg"]]))
    FF <- rbind(FF, t(object[["xreg"]][time, , drop = FALSE]))

  # Prior at t
  a0 <- object[["prior"]][[time]][["a"]]
  R0 <- object[["prior"]][[time]][["R"]]
  # Prior for the observational variance
  n0 <- object[["prior"]][[time]][["n"]]
  s0 <- object[["prior"]][[time]][["s"]]

  # One-step ahead predictive distribution
  f <- drop(crossprod(FF, a0))
  k <- .variance_law(type = object[["variance_law"]][["type"]], mu = f,
                     p = object[["variance_law"]][["power"]])
  q <- s0 * k + drop(crossprod(FF, R0 %*% FF))

  # If data is missing then skip the Kalman updating, posterior = prior
  if (is.na(y)) {
    m <- a0
    C <- R0
    n <- n0
    s <- s0
  } else {
    # Update the posterior parameters:
    e <- y - f

    # Adaptive coefficient vector
    A <- (R0 %*% FF) / q

    # Volatility estimate ratio
    r <- (n0 + e^2 / q) / (n0 + 1)

    # Kalman filter update
    n <- n0 + 1
    s <- r * s0
    m <- drop(a0 + A %*% e)
    C <- r * (R0 - q * tcrossprod(A, A))
  }
  # Get priors a, R for time t + 1 from the posteriors m, C
  P <- tcrossprod(GG %*% C, GG)
  a <- drop(GG %*% m)

  # Discount information
  R <- D * P
  n <- df_variance * n

  colnames(R) <- rownames(R) <- colnames(C)
  names(a) <- names(m)

  list(posterior = list(m = m, C = C, n = n, s = s),
       predictive = list(f = f, q = q, n = n0, s = s0),
       prior = list(a = a, R = R, W = R - P, n = n, s = s))

}

#' @rdname update_moments
#' @export
update_moments.dgegm <- function(object, time) {

  parms_names <- object[["parameters_names"]]
  # Model components
  FF <- object[["FF"]]
  FF_prime <- object[["FF_prime"]]
  g <- object[["g"]]
  GG <- object[["GG"]]
  D <- object[["D"]]
  df_variance <- object[["df_variance"]]
  variance_law <- object[["variance_law"]]
  lambda <- object[["lambda"]]
  # Observed value at time t
  y <- object[["y"]][time]

  # Prior at t
  a0 <- object[["prior"]][[time]][["a"]]
  R0 <- object[["prior"]][[time]][["R"]]
  n0 <- object[["prior"]][[time]][["n"]]
  s0 <- object[["prior"]][[time]][["s"]]

  # Predictive t
  f <- FF(theta = a0, lambda = lambda)
  FF_prime_a <- FF_prime(theta = a0, lambda = lambda)
  q <- drop(crossprod(FF_prime_a, R0 %*% FF_prime_a))
  v <- .variance_law(mu = f, type = variance_law[["type"]],
                     p = variance_law[["p"]])
  Q <- q + v * s0

  # Posterior at t
  e <- y - f

  # Adaptive coefficient vector
  A <- (R0 %*% FF_prime_a) / Q

  # Volatility estimate ratio
  r <- (n0 + e^2 / Q) / (n0 + 1)

  # Kalman filter update
  n <- n0 + 1
  s <- r * s0
  m <- drop(a0 + A %*% e)
  C <- r * (R0 - Q * tcrossprod(A, A))

  # Prior t + 1
  a <- g(m)
  GG_m <- GG(m)
  P <- tcrossprod(GG_m %*% C, GG_m)

  # Discount information
  R <- D * P
  n <- df_variance * n

  colnames(R) <- rownames(R) <- colnames(C) <- parms_names
  names(a) <- names(m) <- parms_names

  list(posterior = list(m = m, C = C, n = n, s = s),
       predictive = list(f = f, q = q, n = n0, s = s0),
       prior = list(a = a, R = R, W = R - P, n = n, s = s))
}

