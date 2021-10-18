#' @importFrom ggplot2 autoplot ggplot geom_line geom_point geom_ribbon scale_color_manual facet_wrap aes
#'
.ggplot_predictive <- function(object, type = c("filter", "smooth", "both"),
                               interval = type == "filter", exclude_prior = FALSE,
                               ts_geom = "line", ts_colour = "black", ts_size = 1,
                               fitted_geom = "line",
                               fitted_colour = "#FF0000", fitted_size = 1,
                               fitted_shape = 21, fitted_fill = "black",
                               fitted_linetype = 1,
                               interval_geom = c("ribbon", "line"),
                               interval_colour = '#0000FF',
                               interval_fill = 'grey69',
                               interval_alpha = 0.6,
                               interval_linetype = "dashed") {
  type <- match.arg(type)

  data_to_plot <- switch (type,
                          "filter" = object$data_predictive_filter,
                          "smooth" = object$data_predictive_smooth,
                          "both" = rbind(object$data_predictive_smooth, object$data_predictive_filter)
  )

  if (exclude_prior)
    data_to_plot <- data_to_plot[!data_to_plot$prior, ]

  geom_ts <- switch (ts_geom,
                     "line" = ggplot2::geom_line,
                     "point" = ggplot2::geom_point
  )
  y <- data_to_plot$y
  g_out <- ggplot(data = data_to_plot, aes(x = t, y = y)) +
    geom_ts(col = ts_colour, size = ts_size)
  y_mean <- data_to_plot$mean
  if (type != "both") {
    geom_fitted <- switch (fitted_geom,
                           "line" = ggplot2::geom_line,
                           "point" = ggplot2::geom_point)
    g_out <- g_out +
      geom_fitted(aes(y = y_mean), col = fitted_colour,
                  size = fitted_size, linetype = fitted_linetype)
    if (interval) {
      ci_lower <- data_to_plot$ci_lower
      ci_upper <- data_to_plot$ci_upper
      interval_geom <- match.arg(interval_geom)
      if (interval_geom == "ribbon") {
        g_out <- g_out +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = interval_fill,
                      alpha = interval_alpha, col = interval_colour)
      }
      if (interval_geom == "line") {
        g_out <- g_out +
          geom_line(aes(y = ci_lower), col = interval_colour,
                    linetype = interval_linetype) +
          geom_line(aes(y = ci_upper), col = interval_colour,
                    linetype = interval_linetype)
      }
    }
  }
  if (type == "both") {
    g_out <- g_out +
      geom_line(data = data_to_plot, aes(y = y_mean, col = type), size = fitted_size,
                linetype = fitted_linetype) +
      scale_color_manual(values = c("#0000FF", "#FF0000"))
  }
  g_out
}

.ggplot_posterior <- function(object, type = c("filter", "smooth", "both"),
                              interval = type == "filter", exclude_prior = FALSE,
                              fitted_colour = "black", fitted_size = 1,
                              fitted_linetype = 1,
                              interval_geom = c("ribbon", "line"),
                              interval_colour = '#0000FF',
                              interval_fill = 'grey69',
                              interval_alpha = 0.6,
                              interval_linetype = "dashed") {
  type <- match.arg(type)

  data_to_plot <- switch (type,
                          "filter" = object$data_posterior_filter,
                          "smooth" = object$data_posterior_smooth,
                          "both" = rbind(object$data_posterior_filter,
                                         object$data_posterior_smooth)
  )

  if (exclude_prior)
    data_to_plot <- data_to_plot[!data_to_plot$prior, ]

  if (type != "both") {
    g_out <- ggplot(data = data_to_plot, aes(x = t, y = mean)) +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line(col = fitted_colour, size = fitted_size, linetype = fitted_linetype)
    if (interval) {
      ci_lower <- data_to_plot$ci_lower
      ci_upper <- data_to_plot$ci_upper
      interval_geom <- match.arg(interval_geom)
      if (interval_geom == "ribbon") {
        g_out <- g_out +
          geom_ribbon(data = data_to_plot, aes(ymin = ci_lower, ymax = ci_upper),
                      fill = interval_fill, alpha = interval_alpha, col = interval_colour)
      }
      if (interval_geom == "line") {
        g_out <- g_out +
          geom_line(data = data_to_plot, aes(y = ci_lower), col = interval_colour,
                    linetype = interval_linetype) +
          geom_line(data = data_to_plot, aes(y = ci_upper), col = interval_colour,
                    linetype = interval_linetype)
      }
    }
  }
  if (type == "both") {
    g_out <- ggplot(data_to_plot, aes(x = t, y = mean, col = type)) +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line(size = fitted_size, linetype = fitted_linetype) +
      scale_color_manual(values = c("#0000FF", "#FF0000"))
  }
  g_out
}
