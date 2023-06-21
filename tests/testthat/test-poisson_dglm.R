test_that("local level model", {

  rm(list = ls())
  devtools::load_all()

  # Simulate some poisson data (local level)
  set.seed(10)
  n <- 100
  var_W <- 0.0015
  lambda <- numeric(n)
  lambda[1] <- log(10)
  omega <- rnorm(n, sd = sqrt(var_W))
  for (t in seq_len(n)[-1]) {
    lambda[t] <- lambda[t - 1] + omega[t]
  }
  y <- rpois(n, lambda = exp(lambda))
  plot(y)
  mean(y)

  model_poisson <- poisson_dglm(polynomial = list(order = 1,
                                                  discount_factor = 0.95))

  out <- forward_filter_poisson_dglm(
    y = as.integer(y),
    F = model_poisson[["FF"]],
    G = model_poisson[["GG"]],
    D = model_poisson[["D"]],
    a = matrix(log(y[1L]), ncol = 1L),
    R = diag(100, 1))

  p <- out$parm2[, 1L] / (1 + out$parm2[, 1L])
  range(p)

  f_median <- qnbinom(p = 0.5, size = out$parm1[, 1L], prob = p)
  f_mean <- out$parm1[, 1L]*(1 - p) / p
  f_mode <- ifelse(out$parm1[, 1L] <= 1, 0, trunc((out$parm1[, 1L] - 1) * (1 - p) / p))
  ci_upper <- qnbinom(p = 1 - 0.05/2, size = out$parm1[, 1L], prob = p)
  ci_lower <- qnbinom(p = 0.05/2, size = out$parm1[, 1L], prob = p)
  plot(y)
  lines(f_median, col = "blue")
  lines(f_mean, col = "red")
  lines(f_mode, col = "green")
  lines(ci_lower, col = "blue", lty = 2)
  lines(ci_upper, col = "blue", lty = 2)

  d <- data.frame(t = 1:length(y), f = out$f[, 1L], q = out$q[, 1L],
                  r = out$parm1[, 1L], beta = out$parm2[, 1L],
                  mean = f_mean, mode = f_mode)
  head(d)

})

test_that("linear growth model", {

  rm(list = ls())
  devtools::load_all()

  # Simulate some poisson data (local level)
  set.seed(6669)
  n <- 100
  # var_level <- 0.002
  # var_growth <- 0.00015
  # th_level <- th_growth <- numeric(n)
  # th_level <- 2#log(5)
  # th_growth <- 0.1
  # omega_level <- rnorm(n, sd = sqrt(var_level))
  # omega_growth <- rnorm(n, sd = sqrt(var_growth))
  # for (t in seq_len(n)[-1]) {
  #   th_growth[t] <- th_growth[t - 1] + omega_growth[t]
  #   th_level[t] <- th_level[t - 1] + th_growth[t] + omega_level[t]
  # }
  # y <- rpois(n, lambda = exp(th_level))
  y <- cumsum(rpois(n, lambda = 2))
  plot(y)
  mean(y)

  model_poisson <- poisson_dglm(polynomial = list(order = 2,
                                                  discount_factor = 0.95))

  out <- forward_filter_poisson_dglm(
    y = as.integer(y),
    F = model_poisson[["FF"]],
    G = model_poisson[["GG"]],
    D = model_poisson[["D"]],
    a = matrix(c(.1, .1), ncol = 2L),
    R = diag(100, 2))

  p <- out$parm2[, 1L] / (1 + out$parm2[, 1L])
  range(p)

  f_median <- qnbinom(p = 0.5, size = out$parm1[, 1L], prob = p)
  f_mean <- out$parm1[, 1L]*(1 - p) / p
  f_mode <- ifelse(out$parm1[, 1L] <= 1, 0, trunc((out$parm1[, 1L] - 1) * (1 - p) / p))
  ci_upper <- qnbinom(p = 1 - 0.10/2, size = out$parm1[, 1L], prob = p)
  ci_lower <- qnbinom(p = 0.10/2, size = out$parm1[, 1L], prob = p)
  plot(y)
  lines(f_median, col = "blue")
  lines(f_mean, col = "red")
  lines(f_mode, col = "green")
  lines(ci_lower, col = "blue", lty = 2)
  lines(ci_upper, col = "blue", lty = 2)

  d <- data.frame(t = 1:length(y), f = out$f[, 1L], q = out$q[, 1L],
                  r = out$parm1[, 1L], beta = out$parm2[, 1L],
                  mean = f_mean, mode = f_mode)
  head(d)

})
