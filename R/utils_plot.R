# plot base -------------------------------------------------------------------------

#' @importFrom graphics par plot lines title grid legend polygon
#' @importFrom grDevices n2mfrow adjustcolor
.plot_predictive <- function(x, type = c("filter", "smooth", "both"),
                             interval = type != "both", exclude_prior = FALSE) {

  cols <- c("t", "mean", "ci_lower", "ci_upper")
  to_plot <- seq_along(x$y)
  if (exclude_prior) to_plot <- to_plot[-seq_len(x$model$prior_length)]
  y <- x$y[to_plot]
  if (type == "both") {
    yl <- range(x$data_predictive_smooth$mean[to_plot], x$data_predictive$mean[to_plot], y)
    plot(y, ylim = yl, type = "b", lty = 2)
    lines(x$data_predictive$mean[to_plot], col = "red", lwd = 2)
    lines(x$data_predictive_smooth$mean[to_plot], col = "blue", lwd = 2)
    legend("topright", legend = c("filter", "smooth"), lwd = 2, col = c('red', 'blue'),
           bty = "n")
    title("Filter and smooth predictive distribution")
    grid()
  }
  else {
    data <- if(type == "filter") x$data_predictive[to_plot, ] else x$data_predictive_smooth[to_plot, ]
    yl <- range(y, data$mean)
    if (interval) yl <- range(yl, data$ci_lower, data$ci_upper)
    plot(y, ylim = yl, type = "b", lty = 2)
    lines(data$mean, col = "blue", lwd = 2)
    if (interval) {
      polygon(x = c(data$t, rev(data$t)),
              y = c(data$ci_lower, rev(data$ci_upper)),
              border = FALSE, col = adjustcolor("#BFBFBF", alpha.f = 0.5))
    }
    title(paste0(type, " predictive distribution"))
    grid()
  }
  invisible()
}

.plot_posterior <- function(x, parm = NULL,
                            type = c("filter", "smooth", "both"),
                            interval = type != "both", exclude_prior = FALSE,
                            mfrow = NULL, mar = NULL) {
  if (is.null(parm))
    parm <- x$model$parameters_names

  to_plot <- seq_along(x$y)
  if (exclude_prior) to_plot <- to_plot[-seq_len(x$model$prior_length)]

  mfrow_orig <- par("mfrow")
  mar_orig <- par("mar")

  y <- x$y[to_plot]
  data <- if(type == "filter") {
    data <- x$data_posterior[x$data_posterior$t %in% to_plot, ]
  }
  if (type == "smooth")
    data <- x$data_posterior_smooth[x$data_posterior_smooth$t %in% to_plot, ]
  if (type == "both") {
    cols <- c("t", "parameter", "mean", "ci_lower", "ci_upper")
    data_filter <- x$data_posterior[x$data_posterior$t %in% to_plot, cols]
    data_filter$type <- "filter"
    data_smooth <- x$data_posterior_smooth[x$data_posterior_smooth$t %in% to_plot, cols]
    data_smooth$type <- "smooth"
    data <- rbind(data_filter, data_smooth)
  }

  # Filter the chosen parameters
  data <- data[data$parameter %in% parm, ]

  if (type == "both") {
    if (length(parm) == 1L) {
      if (parm == "level") {
        yl <- range(y, data$mean)
        plot(y, ylim = yl, type = "b", lty = 2)
        lines(data[data$type == "filter", ]$mean, col = "red", lwd = 2)
        lines(data[data$type == "smooth", ]$mean, col = "blue", lwd = 2)
        legend("topright", legend = c("filter", "smooth"), lwd = 2,
               col = c('red', 'blue'), bty = "n")
        title("Filter and smooth level")
        grid()
      }
      else {
        yl <- range(data$mean)
        plot(data[data$type == "filter", ]$mean, ylim = yl, col = "red", lwd = 2,
             type = "l", ylab = parm)
        lines(data[data$type == "smooth", ]$mean, col = "blue", lwd = 2)
        legend("topright", legend = c("filter", "smooth"), lwd = 2,
               col = c('red', 'blue'), bty = "n")
        title(paste0("Filter and smooth ", parm))
        grid()
      }
    }

    if (length(parm) > 1L) {

      if (is.null(mfrow))
        mfrow <- n2mfrow(length(parm))
      if (is.null(mar))
        mar <- c(4, 4, 1, 1)
      par(mfrow = mfrow, mar = mar)
      for (j in seq_along(parm)) {
        data_cur <- data[data$parameter == parm[j], ]
        yl <- range(data_cur$mean)
        if (parm[j] == "level") {
          yl <- range(yl, y)
          plot(y, ylim = yl, type = "b", lty = 2)
          lines(data_cur[data_cur$type == "filter", ]$mean, col = "red", lwd = 2)
          lines(data_cur[data_cur$type == "smooth", ]$mean, col = "blue", lwd = 2)
          legend("topright", legend = c("filter", "smooth"), lwd = 2,
                 col = c('red', 'blue'), bty = "n")
        } else {
          plot(seq_along(to_plot), data_cur[data_cur$type == "filter", ]$mean, ylim = yl,
               ylab = parm[j], col = "red", lwd = 2, type = "l")
          lines(seq_along(to_plot), data_cur[data_cur$type == "smooth", ]$mean,
                col = "blue", lwd = 2)
          legend("topright", legend = c("filter", "smooth"), lwd = 2,
                 col = c('red', 'blue'), bty = "n")
        }
        title(paste0("Filter and smooth ", parm[j]))
        grid()
      }
    }
  }

  if (type != "both") {
    if (length(parm) == 1L) {
      if (parm == "level") {
        yl <- range(y, data$mean)
        if (interval) yl <- range(yl, data$ci_lower, data$ci_upper)
        plot(y, ylim = yl, type = "b", lty = 2)
        lines(data$mean, col = "blue", lwd = 2)
        if (interval) {
          t_index <- seq_len(nrow(data))
          polygon(x = c(t_index, rev(t_index)),
                  y = c(data$ci_lower, rev(data$ci_upper)),
                  border = FALSE, col = adjustcolor("#BFBFBF", alpha.f = 0.5))
        }
        title(paste0(type, " level"))
        grid()
      }
      else {
        if (interval) yl <- range(data$ci_lower, data$ci_upper)
        plot(data$mean, ylim = yl, ylab = parm, col = "blue", lwd = 2, type = "l")
        if (interval) {
          t_index <- seq_len(nrow(data))
          polygon(x = c(t_index, rev(t_index)),
                  y = c(data$ci_lower, rev(data$ci_upper)),
                  border = FALSE, col = adjustcolor("#BFBFBF", alpha.f = 0.5))
        }
        title(paste0(type, " ", parm))
        grid()
      }
    }
    if (length(parm) > 1L) {
      if (is.null(mfrow))
        mfrow <- n2mfrow(length(parm))
      if (is.null(mar))
        mar <- c(4, 4, 1, 1)
      par(mfrow = mfrow, mar = mar)
      for (j in seq_along(parm)) {
        data_cur <- data[data$parameter == parm[j], ]
        yl <- range(data_cur$mean)
        if (interval) yl <- range(yl, data_cur$ci_lower, data_cur$ci_upper)
        if (parm[j] == "level") {
         yl <- range(yl, y)
         plot(y, ylim = yl, type = "b", lty = 2)
         lines(data_cur$mean, col = "blue", lwd = 2)
        } else {
          plot(data_cur$mean, ylab = parm[j], ylim = yl, type = "l", col = "blue",
               lwd = 2)
        }
        if (interval) {
          t_index <- seq_len(nrow(data_cur))
          polygon(x = c(t_index, rev(t_index)),
                  y = c(data_cur$ci_lower, rev(data_cur$ci_upper)),
                  border = FALSE, col = adjustcolor("#BFBFBF", alpha.f = 0.5))
        }
        title(paste0(type, " ", parm[j]))
        grid()
      }
    }
  }
  par(mfrow = mfrow_orig, mar = mar_orig)
  invisible()
}
