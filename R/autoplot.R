#' @title Fit plot method for Bayesian Dynamic Linear Models.
#'
#' @description Method \code{autoplot} applied to objects of class \code{dlm.fit}.
#' Plot historical data with forecasts and credible interval. Also, plot the filter and
#' smooth state parameters.
#'
#' @param object fitted model object of class \code{dlm.fit}.
#' @param what character. Should plot the \code{"predictive"} or \code{"posterior"}
#' state moments?. Default is \code{"predictive"}.
#' @param type character. Should plot the \code{"filter"}, \code{"smooth"}, or
#' \code{"both"}?. Default is \code{"filter"}.
#' @param interval logical. Should exact credible interval be plotted? Default is
#' \code{TRUE} for filter distribution and FALSE for \code{"smooth"} and \code{"both"}.
#' @param ... additional plotting arguments that affect the plot.
#'
#' @author André F. B. Menezes

#' @rdname autoplot.dlm
#' @export
autoplot.dlm.fit <- function(object, what = c("predictive", "posterior"),
                             type = c("filter", "smooth", "both"), interval = type == "filter", ...) {

  what <- match.arg(what)
  type <- match.arg(type)

  if (what == "predictive") {
    return(.ggplot_predictive(object, type = type, interval = interval, ...))
  }
  if (what == "posterior") {
    return(.ggplot_posterior(object, type = type, interval = interval, ...))
  }
}
