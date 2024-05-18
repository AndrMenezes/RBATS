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
