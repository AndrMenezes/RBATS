#' @title Update the moments of state parameters
#'
#' @description Methods to update the moments of state parameters and compute the
#' one-step ahead predictive distributions for objects of class  \code{dlm} and
#' \code{dgegm}.
#'
#' @param object an object of \code{dgegm} class.
#' @param y vector. Observed value of time series.
#' @param ... currently not used.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' Migon, H. M. and Gamerman, D. (1993). Generalized exponential growth models: A Bayesian approach,  \emph{Journal of Forecasting}, \bold{12}, 573--584

#' @rdname update.dgegm
#' @export
update.dgegm <- function(object, y, ...) {

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

  # Prior at t
  a0 <- object[["prior"]][["a"]]
  R0 <- object[["prior"]][["R"]]
  n0 <- object[["prior"]][["n"]]
  s0 <- object[["prior"]][["s"]]

  # If data is missing then skip the Kalman updating, posterior = prior
  if (is.na(y)) {
    m <- a0
    C <- R0
    n <- n0
    s <- s0
    lpl <- 0
  } else {
    # Predictive t
    f <- FF(theta = a0, lambda = lambda)
    FF_prime_a <- FF_prime(theta = a0, lambda = lambda)
    k <- .variance_law(mu = f, type = variance_law[["type"]],
                       p = variance_law[["p"]])
    q <- s0 * k + drop(crossprod(FF_prime_a, R0 %*% FF_prime_a))

    # Posterior at t
    e <- y - f

    # Adaptive coefficient vector
    A <- (R0 %*% FF_prime_a) / q

    # Volatility estimate ratio
    r <- (n0 + e^2 / q) / (n0 + 1)

    # Kalman filter update
    n <- n0 + 1
    s <- r * s0
    m <- drop(a0 + A %*% e)
    C <- r * (R0 - q * tcrossprod(A, A))

    # Log-predictive likelihood
    lpl <- log(1/sqrt(q) * dt((y - f)/sqrt(q), n0))
  }
  # Prior t + 1
  a <- g(m)
  GG_m <- GG(m)
  P <- tcrossprod(GG_m %*% C, GG_m)

  # Discount information
  R <- D * P
  n <- df_variance * n

  # Set names
  colnames(R) <- rownames(R) <- colnames(C) <- parms_names
  names(a) <- names(m) <- parms_names

  # Save posterior
  object[["posterior"]] <- list(m = m, C = C, n = n, s = s)

  # Update the prior to t + 1
  object[["prior"]] <- list(a = a, R = R, W = R - P, n = n, s = s)

  # Update the log-predictive likelihood
  object[["loglik"]] <- object[["loglik"]] + lpl

  object[["predictive"]] <- list(f = f, q = q)

  object
}

