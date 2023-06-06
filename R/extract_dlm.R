#' @name extract.dlm.fit
#'
#' @title Extract information of filtering and smoothing distribution from
#' \code{dlm.fit} object
#'
#' @description Tidy the results of the forward filter or backward smoother into
#' a \code{data.frame}.
#'
#' @return A \code{data.frame} with the results for each time.
#'
#' @param x an object of class \code{dlm.fit}.
#' @param prob_interval vector with the probability of the credibility interval.
#' @param distribution character. The options are  \code{filter} or \code{smooth}. Default is \code{filter}
#' @param component character. The options are \code{state} or \code{response}.
#' @param ... currently not used.

#' @rdname extract.dlm.fit
#' @export
extract <- function(x, ...) {
  UseMethod("extract", x)
}

#' @rdname extract.dlm.fit
#' @export
extract.dlm.fit <- function(x, prob_interval = c(0.05, 0.20),
                            distribution = c("filter", "smooth"),
                            component = c("response", "state"), ...) {

  distribution <- match.arg(distribution)
  component <- match.arg(component)

  # Model info
  y <- x$y
  nobs <- length(y)
  df_variance <- x$model$df_variance
  seas_type <- x$model$seasonal$type

  data_out <- data.frame(t = seq_along(y), y = y, distribution = distribution)

  # Response (predictive or smoothed mean response)
  if (component == "response") {
    f <- x$filtered$f[, 1L]
    q <- x$filtered$q[, 1L]
    n <- df_variance * x$filtered$n[-(nobs + 1)]
    if (distribution == "smooth") {
      f <- x$smoothed$fk[, 1L]
      q <- x$smoothed$qk[, 1L]
    }
    data_out$mean <- f
    data_out$variance <- q
    data_out$degrees_freedom <- n
    # Compute the credibility intervals
    for (i in seq_along(prob_interval)) {
      p <- prob_interval[i]
      label <- as.character(100 * (1 - p))
      data_out[[paste0("ci_lower__", label)]] <- (
        data_out$mean + qt(
          p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
      data_out[[paste0("ci_upper__", label)]] <- (
        data_out$mean + qt(
          1 - p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
    }
  }

  # State parameters
  if (component == "state") {
    n <- x$filtered$n[-1L]
    parms_names <- x$model$parameters_names
    list_to_append <- list()

    for (j in seq_along(parms_names)) {
      m <- x$filtered$m[j, ]
      C <- x$filtered$C[j, j, ]
      if (distribution == "smooth") {
        m <- x$smoothed$ak[j, ]
        C <- x$smoothed$Rk[j, j, ]
      }
      data_out$parameter <- parms_names[j]
      data_out$mean <- m
      data_out$variance <- C
      data_out$degrees_freedom <- n
      for (i in seq_along(prob_interval)) {
        p <- prob_interval[i]
        label <- as.character(100 * (1 - p))
        data_out[[paste0("ci_lower__", label)]] <- (
          data_out$mean + qt(
            p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
        data_out[[paste0("ci_upper__", label)]] <- (
          data_out$mean + qt(
            1 - p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
      }
      # Appending
      list_to_append[[j]] <- data_out
    }
    # Compute the sum of all seasonal effects
    if (seas_type != "none") {
      i_seas <- x$model$i_seasonal
      F_seas <- x$model$FF[i_seas, ,drop = FALSE]
      m <- if (distribution == "filter") x$filtered$m[i_seas, ] else x$smoothed$ak[i_seas, ]
      C <- if (distribution == "filter") x$filtered$C[i_seas, i_seas, ] else x$smoothed$Rk[i_seas, i_seas, ]
      data_out$parameter <- "sum_seasonality"
      data_out$mean <- crossprod(F_seas, m)[1, ]
      data_out$variance <- apply(C, 3, function(Ct) crossprod(F_seas, Ct %*% F_seas)[1L,1L])
      data_out$degrees_freedom <- n
      for (i in seq_along(prob_interval)) {
        p <- prob_interval[i]
        label <- as.character(100 * (1 - p))
        data_out[[paste0("ci_lower__", label)]] <- (
          data_out$mean + qt(
            p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
        data_out[[paste0("ci_upper__", label)]] <- (
          data_out$mean + qt(
            1 - p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
      }
      list_to_append[[j + 1L]] <- data_out
    }

    data_out <- do.call("rbind", list_to_append)
  }

  data_out
}

