#' @name us_retail_employment
#' @title US monthly retail employment
#' @description Data extract from \code{us_employment} of \pkg{fpp3} package.
#' It contains US retail employment data from January 1990 to June 2019.
#' @format A \code{\link{data.frame}} with 357 observations and 2 columns.
#' @author André Felipe Menezes
#' @usage data(us_retail_employment, package = "RBATS")
#' @source According to \code{fpp3} package: U.S. Bureau of Labor Statistics.
"us_retail_employment"

#' @name vic_electricity_daily
#' @title Daily electricity demand of Victoria, Australia
#' @description Data extract from \code{vic_elec} of \pkg{tsibbledata} package.
#' It contains the daily electricity demand per 1000 and maximum temperature in Victoria, Australia.
#' for 2014.
#' @format A \code{\link{data.frame}} with 367 observations and 5 columns:
#' \itemize{
#' \item \code{time}: Date. The correspond day.
#' \item \code{demand}: Numeric. Total electricity demand per 1000.
#' \item \code{max_temperature}: Numeric. The maximum temperature in that day.
#' \item \code{holiday}: Logical. Indicator for if that day is a public holiday.
#' \item \code{day_type}: Factor. The type of day, it could be \code{holiday}, \code{weekday}, and \code{weekend}.
#' }
#' @author André Felipe Menezes
#' @usage data(us_retail_employment, package = "RBATS")
#' @source According to \code{tsibbledata} package: Australian Energy Market Operator.
"vic_electricity_daily"

#' @name aids_brasil
#' @title Monthly notified cases of AIDS in Brazil
#' @description Data from Gamerman and Migon (1991).
#' It contains the cumulative number of monthly notifications of AID cases in Brazil
#' from September 1985 to December 1988.
#' @format A \code{\link{data.frame}} with 40 observations and 2 columns:
#' \itemize{
#' \item \code{time}: Date. The correspond month.
#' \item \code{cases}: Numeric. Cumulative monthly notifications of AID cases.
#' }
#' @author André Felipe Menezes
#' @usage data(aids_brasil, package = "RBATS")
#' @source Gamerman, D. and Migon, H. S. (1991). Forecasting the number of AIDS cases in Brasil,  \emph{The Statistician}, \bold{40}, 427--442.
"aids_brasil"
