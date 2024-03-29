% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlm.R
\name{dlm}
\alias{dlm}
\title{Bayesian Dynamic Linear Model}
\usage{
dlm(
  polynomial_order = 1L,
  seasonal = list(type = c("none", "free", "fourier"), period = NULL, harmonics = NULL),
  xreg = NULL,
  discount_factors = list(polynomial = 0.95, seasonal = 0.98, regressors = 0.99),
  df_variance = 1,
  variance_law = list(type = "identity", power = 2)
)
}
\arguments{
\item{polynomial_order}{integer. Order of polynomial model.}

\item{seasonal}{list. Components to specify seasonal components.}

\item{xreg}{matrix. Regressors.}

\item{discount_factors}{list. Discount factors for each model components.}

\item{df_variance}{numeric. Discount factor for observation variance. Use a beta-gamma random walk.}

\item{variance_law}{list. Variance law \code{type} and \code{power} parameter.
The variance law \code{type} are \code{identity}, \code{poisson}, \code{binomial},
and \code{power}. The variance law \code{power} should be numeric \eqn{p \geq 1}.}
}
\description{
Create object of class \code{dlm}.
}
\references{
West, M.; Harrison, J. Bayesian Forecasting and Dynamic Models. Springer, 1997.
}
\author{
André F. B. Menezes
}
