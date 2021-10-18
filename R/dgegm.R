#' @title Bayesian Dynamic Generalized Exponential Growth Model
#'
#' @description Create object of class \code{dgegm}.
#'
#' @param lambda numeric. Box-Cox power parameter \eqn{\lambda}. Special cases are:
#' \eqn{\lambda=1} the Modified Exponential curve, \eqn{\lambda=-1} the Logistic curve,
#' and \eqn{\lambda=0} the Gompertz curve.
#' @param discount_factors list. Discount factors for each model parameters,
#' \eqn{\boldsymbol{\theta}_t = (\theta_{1t}, \theta_{2t}, \theta_{3t})}.
#' @param df_variance Discount factor for observation variance. Use a beta-gamma random walk.
#' @param variance_law list. Variance law \code{type} and \code{power} parameter.
#' The variance law \code{type} are \code{identity}, \code{poisson}, \code{binomial},
#' and \code{power}. The variance law \code{power} should be numeric \eqn{p \geq 1}.
#' @param prior list. If specified it should contains prior for state parameters.
#' @param n,s numeric. Prior sample size and mean for the variance, respectively.
#' For both parameter the default is 1.
#'
#' @author André F. B. Menezes
#'
#' @references
#' West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
#'
#' Migon, H. M. and Gamerman, D. (1993). Generalized exponential growth models: A Bayesian approach,  \emph{Journal of Forecasting}, \bold{12}, 573--584
#'
#' Gamerman, D. and Migon, H. S. (1991). Forecasting the number of AIDS cases in Brasil,  \emph{The Statistician}, \bold{40}, 427--442.
#'


#' @rdname dgegm
#' @export
dgegm <- function(lambda = 1, discount_factors, df_variance = 1,
                  variance_law = list(type = "identity", p = 2),
                  prior = NULL, n = 1, s = 1) {

  # Model specification
  FF <- function(theta, lambda) {
    h <- .inverse_link(lambda = lambda)
    h(mu = theta[1L], lambda = lambda)
  }
  FF_prime <- function(theta, lambda) {
    h_prime <- .inverse_link_diff(lambda = lambda)
    a11 <- h_prime(mu = theta[1L], lambda = lambda)
    matrix(c(a11, 0, 0), ncol = 1)
  }
  g <- function(theta) {
    c(theta[1L] + theta[2L], theta[2L]*theta[3L], theta[3L])
  }
  GG <- function(theta) {
    matrix(c(1, 1, 0,
             0, theta[3L], theta[2L],
             0, 0, 1), ncol = 3, nrow = 3, byrow = TRUE)
  }
  # Matrix of discount factors
  D <- diag(x = 1/discount_factors, nrow = 3, ncol = 3)
  D[-which(D==diag(D))] <- 1/0.98

  # Prior
  if (is.null(prior)) {
    prior[[1L]][["a"]] <- NULL
    prior[[1L]][["R"]] <- NULL
    prior[[1L]][["n"]] <- n
    prior[[1L]][["s"]] <- n
  }

  structure(
    list(call = match.call(), lambda = lambda, FF = FF, FF_prime = FF_prime, g = g,
         GG = GG, D = D, prior = prior, df_variance = df_variance,
         variance_law = variance_law, parameters_names = paste0("theta_", 1:3)),
    class = "dgegm")
}
