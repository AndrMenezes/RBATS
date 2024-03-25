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

  cat("Model specification:\n")
  discount_factors <- 1/diag(x$D)
  comp_names <- x[["parameters_names"]]
  tab <- data.frame(parameter = comp_names,
                    discount_factors = discount_factors, row.names = NULL)
  print(tab)
  cat("\n")
  var_law <- if (x$variance_law$type == "power") paste0("power with p = ", x$variance_law$power) else x$variance_law$type
  cat("Variance law:", var_law, "\n")

  if (!is.null(x[["posterior"]])) {
    cat("\nPosterior parameters at time ", x[["time"]], ":\n", sep = "")
    tab_filter <- data.frame(parameter = x[["parameters_names"]],
                             mean = x[["posterior"]][["m"]],
                             variance = diag(as.matrix(x[["posterior"]][["C"]])),
                             row.names = NULL)
    print(tab_filter)
    cat("\n")
    cat("Predictive log-likelihood: ", logLik(x)[1L], "\n", sep = "")
  }
  invisible()
}

#' @rdname s3methods
#' @export
print.dlm.fit <- function(x, ...) {
  print(x$model)
}

#' @rdname s3methods
#' @export
print.dgegm <- function(x, ...) {

  cat("\nBayesian Dynamic Generalized Exponential Growth with lambda = ", x$lambda,
      "\n\n", sep = "")

  cat(paste0("Call:\n", paste(deparse(x$call, width.cutoff = 80L), collapse = "\n")),
      "\n\n")

  cat("Model specification:\n")
  discount_factors <- 1/diag(x$D)
  comp_names <- x[["parameters_names"]]
  tab <- data.frame(parameter = comp_names,
                    discount_factors = discount_factors, row.names = NULL)
  print(tab)
  cat("\n")
  var_law <- if (x$variance_law$type == "power") paste0("power with p = ", x$variance_law$power) else x$variance_law$type
  cat("Variance law:", var_law, "\n")

  if (!is.null(x[["prior"]])) {
    cat("\nPosterior parameters at time ", x[["time"]], ":\n", sep = "")
    tab_filter <- data.frame(parameter = x[["parameters_names"]],
                             mean = x[["posterior"]][["m"]],
                             variance = diag(x[["posterior"]][["C"]]), row.names = NULL)
    print(tab_filter)
    cat("\n")
    cat("Predictive log-likelihod: ", logLik(x)[1L], "\n", sep = "")
  }
  invisible()
}


#' @rdname s3methods
#' @export
print.dgegm.fit <- function(x, ...) {
  print(x$model)
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
logLik.dlm <- function(object, ...) {
  if (!missing(...))
    warning("Extra arguments discarded")
  ll <- object$loglik
  attr(ll, "df") <- length(object$parameters_names)
  attr(ll, "nobs") <- object$time
  class(ll) <- "logLik"
  ll
}

#' @rdname s3methods
#' @export
logLik.dlm.fit <- function(object, ...) {
  logLik(object[["model"]])
}
