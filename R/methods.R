#' @name s3methods
#' @title S3 methods for \code{dlm} objects
#'
#' @description Usual S3 methods for extracting information from objects of class
#' \code{\link{dlm}}.
#'
#' @param x,object an object of class \code{\link{dlm}}.
#' @param ... currently not used.
#'
#' @importFrom utils head
#' @importFrom stats dt logLik

#' @rdname s3methods
#' @export
print.dlm <- function(x, ...) {

  cat("\nBayesian Dynamic Linear Model \n\n", sep = "")

  cat(paste0("Call:\n", paste(deparse(x$call, width.cutoff = 80L), collapse = "\n")),
      "\n\n")

  cat("Model specfication:\n")
  discount_factors <- 1/diag(x$D)
  comp_names <- x[["parameters_names"]]
  tab <- data.frame(parameter = comp_names,
                    discount_factors = discount_factors, row.names = NULL)
  print(tab)
  cat("\n")
  var_law <- if (x$variance_law$type == "power") paste0("power with p = ", x$variance_law$power) else x$variance_law$type
  cat("Variance law:", var_law, "\n")
}

#' @rdname s3methods
#' @export
print.dgegm <- function(x, ...) {

  cat("\nBayesian Dynamic Generalized Exponential Growth with lambda = ", x$lambda,
      "\n\n", sep = "")

  cat(paste0("Call:\n", paste(deparse(x$call, width.cutoff = 80L), collapse = "\n")),
      "\n\n")

  cat("Model specfication:\n")
  discount_factors <- 1/diag(x$D)
  comp_names <- x[["parameters_names"]]
  tab <- data.frame(parameter = comp_names,
                    discount_factors = discount_factors, row.names = NULL)
  print(tab)
  cat("\n")
  var_law <- if (x$variance_law$type == "power") paste0("power with p = ", x$variance_law$power) else x$variance_law$type
  cat("Variance law:", var_law, "\n")
}

#' @rdname s3methods
#' @export
print.dlm.fit <- function(x, ...) {

  cat("\nFit from Bayesian Dynamic Linear Model \n\n", sep = "")

  cat(paste0("Call:\n", paste(deparse(x$call, width.cutoff = 80L), collapse = "\n")),
      "\n\n")

  cat("The model used the first", x[["prior_length"]],
      "observations to construct the prior parameters\n\n")

  last_index <- length(x[["y"]])
  cat("Filtering parameters at time ", last_index, ":\n", sep = "")
  last_posterior <- x[["posterior"]][[last_index]]
  tab_filter <- data.frame(parameter = names(last_posterior[["m"]]),
                           mean = last_posterior[["m"]],
                           variance = diag(last_posterior[["C"]]), row.names = NULL)
  print(tab_filter)
  cat("\n")
  cat("Predictive log-likelihod: ", logLik(x)[1L], "\n", sep = "")
}

#' @rdname s3methods
#' @export
print.dgegm.fit <- function(x, ...) {

  cat("\nFit from Bayesian Dynamic Generalized Exponential Growth Model \n\n", sep = "")

  cat(paste0("Call:\n", paste(deparse(x$call, width.cutoff = 80L), collapse = "\n")),
      "\n\n")

  last_index <- length(x[["y"]])
  cat("Filtering parameters at time ", last_index, ":\n", sep = "")
  last_posterior <- x[["posterior"]][[last_index]]
  tab_filter <- data.frame(parameter = names(last_posterior[["m"]]),
                           mean = last_posterior[["m"]],
                           variance = diag(last_posterior[["C"]]), row.names = NULL)
  print(tab_filter)
  cat("\n")
  cat("Predictive log-likelihod: ", logLik(x)[1L], "\n", sep = "")
}

#' @rdname s3methods
#' @export
print.dlm.forecast <- function(x, ...) {
  cat("\nForecast from Bayesian Dynamic Linear Model \n\n", sep = "")
  cat(paste0("Call:\n", paste(deparse(x$call, width.cutoff = 80L), collapse = "\n")),
      "\n\n")
  h <- nrow(x$predictive)
  cat(h, "-step ahead forecast distribution parameters:\n", sep = "")
  print(head(x$predictive, 5))
  if (h > 5)
    message("# ... with ", h - 5, " more rows")
  cat("\n")
}

#' @rdname s3methods
#' @export
logLik.dlm.fit <- function(object, ...) {
  mu <- sapply(object[["predictive"]], function(z) z[["f"]])
  sigma <- sqrt(sapply(object[["predictive"]], function(z) z[["q"]]))
  df <- sapply(object[["predictive"]], function(z) z[["n"]])
  ll <- sum(log(1/sigma * dt((object[["y"]] - mu)/sigma, df)), na.rm = TRUE)
  attr(ll, "nobs") <- length(mu)
  class(ll) <- "logLik"
  ll
}

#' @rdname s3methods
#' @export
logLik.dgegm.fit <- function(object, ...) {
  mu <- sapply(object[["predictive"]], function(z) z[["f"]])
  sigma <- sqrt(sapply(object[["predictive"]], function(z) z[["q"]]))
  df <- sapply(object[["predictive"]], function(z) z[["n"]])
  ll <- sum(log(1/sigma * dt((object[["y"]] - mu)/sigma, df)), na.rm = TRUE)
  attr(ll, "nobs") <- length(mu)
  class(ll) <- "logLik"
  ll
}
