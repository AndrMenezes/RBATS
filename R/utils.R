.variance_law <- function(type, mu, p) {
  switch(type,
        "identity" = 1,
        "poisson" = mu,
        "binomial" = mu * (1 - mu),
        "power" = mu ^ p
  )
}

.diag_one <- function(x) {
  m <- diag(x, nrow = length(x), ncol = length(x))
  m[which(m != diag(m))] <- 1
  m
}

.matrix_power <- function(x, k) {
  x_k <- x
  i <- 2L
  while (i <= k) {
    x_k <- x_k %*% x
    i <- i + 1L
  }
  x_k
}

.inverse_link <- function(lambda) {
  if (lambda == 0) return(function(mu, lambda) exp(mu))
  else return(function(mu, lambda) mu^(1/lambda))
}

.inverse_link_diff <- function(lambda) {
  if (lambda == 0) return(function(mu, lambda) exp(mu))
  else return(function(mu, lambda) mu^(1/lambda - 1) / lambda)
}


# Organize the predictive distribution parameters in data.frame
#' @importFrom stats qt
.tidy_predictive <- function(list_predictive, interval = TRUE, level = 0.10) {

  data_output <- data.frame(t = seq_along(list_predictive))
  data_output$mean <- sapply(list_predictive, function(z) unname(z[["f"]]))
  data_output$variance <- sapply(list_predictive, function(z) unname(z[["q"]]))
  data_output$df <- sapply(list_predictive, function(z) unname(z[["n"]]))

  if (interval) {
    data_output$ci_lower <- data_output$mean + qt(level / 2, df = data_output$df) * sqrt(data_output$variance)
    data_output$ci_upper <- data_output$mean + qt(1 - level / 2, df = data_output$df) * sqrt(data_output$variance)
  }
  data_output
}

# Organize the posterior distribution parameters in data.frame
#' @importFrom stats qt qnorm
.tidy_posterior <- function(list_posterior, interval = TRUE, level = 0.10) {

  # Parameters names
  parms_names <- names(list_posterior[[1L]][["m"]])

  # Time index
  time_index <- seq_along(list_posterior)

  # Degrees of freedom
  df <- sapply(list_posterior, function(z) unname(z[["n"]]))

  lt_parms <- lapply(parms_names, function(parm) {
    tmp <- lapply(list_posterior, function(z) {
      diag_C <- diag(z[["C"]])
      c("mean" = unname(z[["m"]][which(names(z[["m"]]) == parm)]),
        "variance" = unname(diag_C[names(diag_C) == parm]))
    })
    values <- do.call(rbind, tmp)
    out <- data.frame(t = time_index)
    out$parameter <- parm
    out$mean <- values[, "mean"]
    out$variance <- values[, "variance"]
    out$df <- df
    out
  })
  data_output <- do.call(rbind, lt_parms)

  if (interval) {
    data_output$ci_lower <- data_output$mean + qt(level / 2, df = data_output$df) * sqrt(data_output$variance)
    data_output$ci_upper <- data_output$mean + qt(1 - level / 2, df = data_output$df) * sqrt(data_output$variance)
  }
  data_output
}

# Organize the posterior distribution parameters in data.frame
#' @importFrom stats qt
.tidy_smooth <- function(list_smooth, what = c("posterior", "predictive"),
                         interval = TRUE, level = 0.10) {

  # Parameters names
  parms_names <- names(list_smooth[[1L]][["a"]])


  # Final degrees of freedom and variance
  nT <- list_smooth[["nT"]]
  sT <- list_smooth[["sT"]]
  list_smooth[["sT"]] <- NULL
  list_smooth[["nT"]] <- NULL

  # Time index
  time_index <- seq_along(list_smooth)

  if (what == "posterior") {
    lt_parms <- lapply(parms_names, function(parm) {
      tmp <- lapply(seq_along(list_smooth), function(j) {
        z <- list_smooth[[j]]
        diag_R <- diag(z[["R"]])
        c("mean" = unname(z[["a"]][which(names(z[["a"]]) == parm)]),
          "variance" = unname(diag_R[names(diag_R) == parm]),
          "s" = unname(z[["s"]]))
      })
      values <- do.call(rbind, tmp)
      out <- data.frame(t = time_index)
      out$parameter <- parm
      out$mean <- values[, "mean"]
      out$variance <- values[, "variance"]
      out$s <- values[, "s"]
      out
    })
    data_output <- do.call(rbind, lt_parms)
  }
  if (what == "predictive") {
    data_output <- data.frame(t = time_index)
    data_output$mean <- sapply(list_smooth, function(z) unname(z[["f"]]))
    data_output$variance <- sapply(list_smooth, function(z) unname(z[["q"]]))
    data_output$s <- sapply(list_smooth, function(z) unname(z[["s"]]))
  }
  if (interval) {
    qa <- qnorm(1 - level / 2)
    data_output$ci_lower <- data_output$mean - qa * sqrt(data_output$variance)
    data_output$ci_upper <- data_output$mean + qa * sqrt(data_output$variance)
  }
  data_output
}

