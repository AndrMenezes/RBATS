test_that("multiplication works", {

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

  a <- matrix(log(y[1L]), ncol = 1L)
  R <- diag(100, 1)
  update_poisson_dglm(y = as.integer(y[1L]),
                      F = model_poisson[["FF"]],
                      G = model_poisson[["GG"]],
                      D = model_poisson[["D"]],
                      a = matrix(log(y[1L]), ncol = 1L),
                      R = diag(100, 1))
  forward_filter_poisson_dglm(
    y = as.integer(y),
    F = model_poisson[["FF"]],
    G = model_poisson[["GG"]],
    D = model_poisson[["D"]],
    a = matrix(log(y[1L]), ncol = 1L),
    R = diag(100, 1))
  print(a)
  print(R)






})
