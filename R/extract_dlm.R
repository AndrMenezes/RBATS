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
  df_variance <- x$model$df_variance
  seas_type <- x$model$seasonal$type

  # Template
  data_out <- data.frame(t = seq_along(y), y = y, distribution = distribution)

  # Response (predictive or smoothed mean response)
  if (component == "response") {
    data_out$mean <- if (distribution == "filter") x$filtered$f[, 1L] else x$smoothed$fk[, 1L]
    data_out$variance <- if (distribution == "filter") x$filtered$q[, 1L] else x$smoothed$qk[, 1L]
    data_out$degrees_freedom <- df_variance * (x$filtered$n - 1)
  }

  # State parameters
  if (component == "state") {
    parms_names <- x$model$parameters_names
    list_to_append <- list()
    for (j in seq_along(parms_names)) {
      data_out$parameter <- parms_names[j]
      data_out$mean <- if (distribution == "filter") x$filtered$m[j, ] else x$smoothed$ak[j, ]
      data_out$variance <- if (distribution == "filter") x$filtered$C[j, j, ] else x$smoothed$Rk[j, j, ]
      data_out$degrees_freedom <- x$filtered$n
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
      data_out$mean <- crossprod(F_seas, m)[1L, ]
      data_out$variance <- apply(C, 3L, function(Ct) crossprod(F_seas, Ct %*% F_seas)[1L, 1L])
      data_out$degrees_freedom <- x$filtered$n
      list_to_append[[j + 1L]] <- data_out
    }
    data_out <- do.call("rbind", list_to_append)
  }

  # Compute the credibility intervals
  for (p in prob_interval) {
    label <- as.character(100 * (1 - p))
    data_out[[paste0("ci_lower__", label)]] <- (
      data_out$mean + qt(
        p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
    data_out[[paste0("ci_upper__", label)]] <- (
      data_out$mean + qt(
        1 - p / 2, df = data_out$degrees_freedom) * sqrt(data_out$variance))
  }

  data_out
}

