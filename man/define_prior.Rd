% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_prior.R
\name{define_prior}
\alias{define_prior}
\alias{define_prior.dlm}
\alias{define_prior.dgegm}
\title{Define prior moments for the state parameters.}
\usage{
define_prior(object, ...)

\method{define_prior}{dlm}(object, y, ...)

\method{define_prior}{dgegm}(object, y, ...)
}
\arguments{
\item{object}{an object of class \code{dlm} or \code{dgegm}.}

\item{...}{currently not used.}

\item{y}{numeric. prior first observations use to compute the prior moments.}
}
\description{
Helper functions applied on objects of class \code{dlm} and \code{dgegm},
to define a Dynamic Linear Models (DLM) or Dynamic Generalized Exponential Growth Models
(DGEGM), respectively. This function is especially useful if you do not know how to
specify a prior mean, \eqn{\mathbf{a}_0}, and covariance matrix,
\eqn{\mathbf{R}_0}, for the state vector.
}
\details{
For \code{dlm} objects use linear model to define the prior
moments, the mean and covariance matrix, of model components.

For \code{dgegm} objects the prior moments are:
for \eqn{\theta_1} is \eqn{g(y_1)}, for \eqn{\theta_2} is 0.001, and for \eqn{\theta_3}
is 0.95, where \eqn{g(\cdot)} is the Box-Cox link function and \eqn{y_1} is the first
observation of time series. The covariance
}
\author{
André F. B. Menezes
}
