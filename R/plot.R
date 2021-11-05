#' @title Plot method for Bayesian Dynamic Linear Models.
#'
#' @description Method \code{plot} applied to objects of class \code{dlm.fit}.
#' Plot historical data with forecasts and credible interval. Also, plot the filter and
#' smooth state parameters.
#'
#' @param x fitted model object of class \code{dlm.fit}.
#' @param what character. Should plot the \code{"predictive"} or \code{"posterior"}
#' state moments?. Default is \code{"predictive"}.
#' @param type character. Should plot the \code{"filter"}, \code{"smooth"}, or
#' \code{"both"}?. Default is \code{"filter"}.
#' @param interval logical. Should exact credible interval be plotted? Default is
#' \code{TRUE} for filter distribution and FALSE for \code{"smooth"} and \code{"both"}.
#' @param parm vector. a specification of which parameters are to be plotted, a vector of names. If missing, all parameters are considered.
#' @param ... additional plotting arguments that affect the plot.
#'
#' @author André F. B. Menezes

#' @rdname plot.dlm
#' @export
plot.dlm.fit <- function(x, what = c("predictive", "posterior"),
                         type = c("filter", "smooth", "both"),
                         interval = type != "both", parm = NULL, ...) {

  what <- match.arg(what)
  type <- match.arg(type)

  if (what == "predictive") {
    return(.plot_predictive(x, type = type, interval = interval, ...))
  }
  if (what == "posterior") {
    return(.plot_posterior(x, type = type, interval = interval, parm = parm, ...))
  }
}
