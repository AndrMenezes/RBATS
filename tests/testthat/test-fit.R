test_that("filter and smooth functions for growth model for Nile data works", {
  y <- c(Nile)

  # Define the model
  model <- dlm(
    polynomial = list(order = 1, discount_factor = 0.90)
  )

  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1]), ncol = 1),
                       C0 = diag(100, 1))

  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
})

test_that("filter and smooth functions with missing value for Nile data works", {
  y <- c(Nile)
  y[c(10, 30)] <- NA_real_
  model <- dlm(
    polynomial = list(order = 1, discount_factor = 0.90)
  )

  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1]), ncol = 1),
                       C0 = diag(100, 1))

  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
})

test_that("level and ar for Nile data works", {
  y <- c(Nile)
  y[c(10, 30)] <- NA_real_
  model <- dlm(
    polynomial = list(order = 1, discount_factor = 0.90),
    autoregressive = list(order = 1, discount_factor = 0.998)
  )

  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], 0, 0), ncol = 1),
                       C0 = diag(c(100, 2, 1)))
  # plot(y)
  # lines(fitted_object$filtered$f[, 1L], col = "blue")
  # lines(fitted_object$smoothed$fk[, 1L], col = "red")
  # lines(fitted_object$filtered$m[1L, ], col = "green")
  # fitted_object$filtered$m[3L, ]

  model <- dlm(
    polynomial = list(order = 1, discount_factor = 0.90),
    autoregressive = list(order = 2, discount_factor = 0.998)
  )

  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], 0, 0, 0, 0), ncol = 1),
                       C0 = diag(c(100, 2, 2, 1, 1)))
  # plot(y)
  # lines(fitted_object$filtered$f[, 1L], col = "blue")
  # lines(fitted_object$smoothed$fk[, 1L], col = "red")
  # lines(fitted_object$filtered$m[1L, ], col = "green")

  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
})

test_that("trend seasonal model for AirPassengers data works", {
  y <- c(AirPassengers)

  # Seasonal fourier (sin and cos)
  model_object <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 12, harmonics = c(1, 2, 3),
                    discount_factor = 0.98)
  )

  fitted_object <- fit(model = model_object, y = y,
                       m0 = matrix(c(y[1], rep(0, 7)), ncol = 1),
                       C0 = diag(100, 8))
  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
  # plot(y)
  # lines(fitted_object$filtered$f[, 1L], col = "blue")
  # lines(fitted_object$smoothed$fk[, 1L], col = "red")
  # lines(fitted_object$filtered$m[1L, ], col = "green")

  # Seasonal free
  model_object <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "free", period = 12,
                    discount_factor = 0.99)
  )
  fitted_object <- fit(model = model_object, y = y,
                       m0 = matrix(c(y[1], rep(0, 12)), ncol = 1),
                       C0 = diag(100, 13))
  # plot(y)
  # lines(fitted_object$filtered$f[, 1L], col = "blue")
  # lines(fitted_object$smoothed$ak[1L, ], col = "red")

  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
})

test_that("trend seasonal model for us_retail_employment data works", {
  data(us_retail_employment)
  y <- us_retail_employment$employed

  model_object <- dlm(
    polynomial = list(order = 2, discount_factor = 0.90),
    seasonal = list(type = "fourier", period = 12, harmonics = c(1, 2, 3),
                    discount_factor = 0.98),
    autoregressive = list(order = 2, discount_factor = 0.998)
  )
  fitted_object <- fit(model = model_object, y = y,
                       m0 = matrix(c(y[1], rep(0, 7), rep(0, 4)), ncol = 1),
                       C0 = diag(c(rep(100, 8), 2, 2, 1, 1), 12))
  # plot(y)
  # lines(fitted_object$filtered$f[, 1L], col = "red")
  # lines(fitted_object$smoothed$ak[1L, ], col = "blue")

  fitted_object$filtered$m[11L, length(y)]
  fitted_object$filtered$m[12L, length(y)]

  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")


})

test_that("ar(2) model with simulated data", {

  nobs <- 200
  sd_y <- 0.2
  true_phi_1 <- 0.9
  true_phi_2 <- -0.5

  xi_1 <- numeric(nobs)
  xi_2 <- numeric(nobs)
  y <- numeric(nobs)
  xi_1[1L] <- 0

  # Random errors
  set.seed(121)
  nu <- rnorm(n = nobs, sd = sd_y)

  # First observation
  y[1L] <- xi_1[1L] + nu[1L]

  for (t in seq_len(nobs)[-1]) {
    xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
    xi_2[t] <- xi_1[t - 1]
    y[t] <- xi_1[t]
  }

  mod <- dlm(autoregressive = list(order = 2, discount_factor = 0.998))
  m0 <- matrix(c(0, 0, 1, 0), ncol = 1)
  C0 <- diag(x = c(2, 2, 1, 1), nrow = 4)
  out <- fit(model = mod, y = y, m0 = m0, C0 = C0)
  expect_equal(out$filtered$m[3L, nobs] - true_phi_1, -0.0792, tolerance = 1e-4)
  expect_equal(out$filtered$m[4L, nobs] - true_phi_2, 0.0597, tolerance = 1e-4)

  # plot(y)
  # lines(out$filtered$f[, 1L], col = "blue")
  # lines(out$smoothed$fk[, 1L], col = "red")
  # # Phi's
  # plot(out$filtered$m[3L, ]); abline(h = true_phi_1)
  # lines(out$smoothed$ak[3L, ], col = "blue")
  # plot(out$filtered$m[4L, ]); abline(h = true_phi_2)
  # lines(out$smoothed$ak[4L, ], col = "blue")


})

test_that("level + ar(2) with simulated data", {

  # Simulate the data
  nobs <- 100
  sd_y <- 1
  sd_mu <- 0.06
  true_phi_1 <- 0.8
  true_phi_2 <- -0.5

  mu <- numeric(nobs)
  xi_1 <- numeric(nobs)
  xi_2 <- numeric(nobs)
  y <- numeric(nobs)
  xi_1[1L] <- 0
  mu[1L] <- 5

  # Random errors
  set.seed(1111)
  nu <- rnorm(n = nobs, sd = sd_y)
  omega_mu <- rnorm(n = nobs, sd = sd_mu)

  # First observation
  y[1L] <- mu[1L] + xi_1[1L] + nu[1L]

  for (t in seq_len(nobs)[-1]) {
    xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
    xi_2[t] <- xi_1[t - 1]
    mu[t] <- mu[t - 1] + omega_mu[t]
    y[t] <- mu[t] + xi_1[t]
  }

  # Model definition
  mod <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    autoregressive = list(order = 2, discount_factor = 0.998)
  )
  out <- fit(model = mod, y = y,
             m0 = matrix(c(y[1L], rep(0, 4)), ncol = 1),
             C0 = diag(c(10, 2, 2, 1, 1), nrow = 5))

  expect_equal(out$filtered$m[4L, nobs] - true_phi_1, 0.0447, tolerance = 1e-3)
  expect_equal(out$filtered$m[5L, nobs] - true_phi_2, 0.0463, tolerance = 1e-3)

  # plot(out$filtered$m[4L, ]); abline(h = true_phi_1)
  # lines(out$smoothed$ak[4L, ], col = "blue")
  # plot(out$filtered$m[5L, ]); abline(h = true_phi_2)
  # lines(out$smoothed$ak[5L, ], col = "blue")
  #
  # plot(mu, ylim = c(4.5, 7))
  # lines(out$filtered$m[1L, ], col = "blue")
  # lines(out$smoothed$ak[1L, ], col = "red")
  #
  # plot(y)
  # lines(out$filtered$f[, 1L], col = "blue")
  # lines(out$smoothed$fk[, 1L], col = "red")


})

test_that("level + seasonal + ar(2) model with simulated data", {

  # Setting parameters
  nobs <- 500
  sd_y <- 1
  sd_mu <- 0.02
  sd_a <- 1e-6
  sd_b <- 1e-6
  true_phi_1 <- 0.8
  true_phi_2 <- -0.5

  y <- numeric(nobs)
  mu <- numeric(nobs)
  a <- numeric(nobs)
  b <- numeric(nobs)
  xi_1 <- numeric(nobs)
  xi_2 <- numeric(nobs)
  xi_1[1L] <- 0
  mu[1L] <- 5
  a[1L] <- 0.1
  b[1L] <- 0.1

  period <- 12L # annual
  r <- 1L # number of harmonics
  w <- 2 * pi * r / period

  # Random errors
  set.seed(1111)
  nu <- rnorm(n = nobs, sd = sd_y)
  omega_mu <- rnorm(n = nobs, sd = sd_mu)
  omega_a <- rnorm(n = nobs, sd = sd_a)
  omega_b <- rnorm(n = nobs, sd = sd_b)

  # First observation
  y[1L] <- mu[1L] + a[1L] + xi_1[1L] + nu[1L]

  for (t in seq_len(nobs)[-1]) {
    # AR2
    xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
    xi_2[t] <- xi_1[t - 1]
    # Level
    mu[t] <- mu[t - 1] + omega_mu[t]
    # Seasonality
    a[t] <- a[t - 1] * cos(w) + b[t - 1] * sin(w) + omega_a[t]
    b[t] <- b[t - 1] * cos(w) - a[t - 1] * sin(w) + omega_b[t]
    # Response
    y[t] <- mu[t] + a[t] + xi_1[t]
  }

  mod <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 12, harmonics = 1,
                    discount_factor = 0.99),
    autoregressive = list(order = 2, discount_factor = 0.998)
  )
  m0 <- matrix(c(y[1L], rep(0, 6)), ncol = 1)
  C0 <- diag(x = c(4, 2, 2, rep(1, 4)), nrow = 7L)
  out <- fit(model = mod, y = y, m0 = m0, C0 = C0)

  expect_equal(out$filtered$m[6L, nobs] - true_phi_1, -0.0149, tolerance = 1e-3)
  expect_equal(out$filtered$m[7L, nobs] - true_phi_2, 0.0334, tolerance = 1e-3)

  # plot(out$filtered$m[6L, ]); abline(h = true_phi_1)
  # lines(out$smoothed$ak[6L, ], col = "blue")
  # plot(out$filtered$m[7L, ]); abline(h = true_phi_2)
  # lines(out$smoothed$ak[7L, ], col = "blue")

})

test_that("level + regression + ar(2)", {


  # Simulate the data
  nobs <- 500
  sd_y <- 1
  sd_mu <- 0.06
  sd_beta_1 <- 0.003
  true_phi_1 <- 0.8
  true_phi_2 <- -0.5

  mu <- numeric(nobs)
  xi_1 <- numeric(nobs)
  xi_2 <- numeric(nobs)
  beta_1 <- numeric(nobs)
  y <- numeric(nobs)
  xi_1[1L] <- 0
  mu[1L] <- 5
  beta_1[1L] <- 0.5

  set.seed(1111)
  # Covariate
  x1 <- rexp(nobs, 1)

  # Random errors
  nu <- rnorm(n = nobs, sd = sd_y)
  omega_mu <- rnorm(n = nobs, sd = sd_mu)
  omega_beta_1 <- rnorm(n = nobs, sd = sd_beta_1)

  # First observation
  y[1L] <- mu[1L] + xi_1[1L] + beta_1[1L] * x1[1L] + nu[1L]

  for (t in seq_len(nobs)[-1]) {
    xi_1[t] <- true_phi_1 * xi_1[t - 1] + true_phi_2 * xi_2[t - 1] + nu[t]
    xi_2[t] <- xi_1[t - 1]
    beta_1[t] <- beta_1[t - 1] + omega_beta_1[t]
    mu[t] <- mu[t - 1] + omega_mu[t]
    y[t] <- mu[t] + xi_1[t] + beta_1[t] * x1[t]
  }

  mod <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    regressor = list(xreg = matrix(x1, ncol = 1), discount_factor = 0.998),
    autoregressive = list(order = 2, discount_factor = 0.998)
  )
  m0 <- matrix(c(y[1L], 0, 0, 0, 0, 0), ncol = 1)
  C0 <- diag(x = c(10, 4, 2, 2, 1, 1), nrow = 6L)
  out <- fit(model = mod, y = y, m0 = m0, C0 = C0)

  expect_equal(out$filtered$m[5L, nobs] - true_phi_1, 0.0799, tolerance = 1e-3)
  expect_equal(out$filtered$m[6L, nobs] - true_phi_2, -0.0487, tolerance = 1e-3)

  # plot(beta_1)
  # lines(out$filtered$m[2L, ], col = "blue")
  # lines(out$smoothed$ak[2L, ], col = "red")

  # plot(y)
  # lines(out$filtered$f[, 1L], col = "blue")
  # lines(out$smoothed$fk[, 1L], col = "red")
  #
  # plot(out$filtered$m[5L, ]); abline(h = true_phi_1)
  # lines(out$smoothed$ak[5L, ], col = "blue")
  # plot(out$filtered$m[6L, ]); abline(h = true_phi_2)
  # lines(out$smoothed$ak[6L, ], col = "blue")



})

test_that("TF(1) for simulated data", {


  n_obs <- 100
  sd_y <- 1
  sd_E <- 0.02
  sd_psi <- 0.002
  lambda <- 0.7


  E <- numeric(n_obs)
  psi <- numeric(n_obs)
  y <- numeric(n_obs)
  E[1L] <- 0
  psi[1L] <- 2.5

  # Covariate
  set.seed(1111)
  x <- rweibull(n_obs, 2, 0.4)

  # Random errors
  nu <- rnorm(n = n_obs, sd = sd_y)
  omega_psi <- rnorm(n = n_obs, sd = sd_psi)
  omega_E <- rnorm(n = n_obs, sd = sd_E)

  # First observation
  y[1L] <- E[1L] + nu[1L]

  for (t in seq_len(n_obs)[-1]) {
    psi[t] <- psi[t - 1] + omega_psi[t]
    E[t] <- lambda * E[t - 1] + psi[t] * x[t] + omega_E[t]
    y[t] <- E[t] + nu[t]
  }

  # Estimate lambda by using 1st Taylor approximation
  mod_taylor <- dlm(transfer_function = list(order = 1, xreg = x,
                                             discount_factor = 1))

  m0 <- c(0, 0, 0)
  C0 <- diag(x = c(1, 0.5, 10), nrow = 3)
  out <- fit(model = mod_taylor, y = y, m0 = m0, C0 = C0)
  plot(out$filtered$m[2L, ])

  # Estimate for fixed value of lambda
  mod_fixed <- dlm(transfer_function = list(order = 1, xreg = x,
                                            lambda = 0.7,
                                            discount_factor = 1))
  out_fixed <- fit(model = mod_fixed, y = y, m0 = matrix(c(0, 0), ncol = 1),
                   C0 = diag(c(1, 10)))
  plot(out_fixed$filtered$m[2L, ])
  lines(out$filtered$m[3L, ], col = "blue")

  # Estimate different models for fixed lambdas and evaluate the LPL
  lambda_seq <- seq(-1, 1, l = 10)
  list_fits <- lapply(lambda_seq, function(lbd) {
    mod_fixed <- dlm(transfer_function = list(order = 1, xreg = x,
                                              lambda = lbd,
                                              discount_factor = 1))
    fit(model = mod_fixed, y = y,
        m0 = matrix(c(0, 0), ncol = 1),
        C0 = diag(c(1, 10)))
  })
  names(list_fits) <- lambda_seq
  sort(sapply(list_fits, logLik))

})


test_that("TF(2) model with simulated data", {

  # Simulate the data
  npulse <- 5
  nobs <- 50
  sd_y <- 1
  sd_E1 <- 0.0004
  sd_psi <- 0.002
  true_lambda_1 <- 0.8
  true_lambda_2 <- -0.2

  E1 <- numeric(nobs)
  E2 <- numeric(nobs)
  psi <- numeric(nobs)
  y <- numeric(nobs)
  E1[1L] <- 0
  psi[1L] <- 2.5

  # Covariate
  x <- numeric(nobs) # rpois(nobs, lambda = 1)
  impulse <- seq(1, nobs, by = nobs / npulse)
  x[impulse] <- 1

  # Random errors
  set.seed(1111)
  nu <- rnorm(n = nobs, sd = sd_y)
  omega_psi <- rnorm(n = nobs, sd = sd_psi)
  omega_E1 <- rnorm(n = nobs, sd = sd_E1)

  # First observation
  y[1L] <- E1[1L] + nu[1L]

  for (t in seq_len(nobs)[-1]) {
    psi[t] <- psi[t - 1] + omega_psi[t]
    E1[t] <- true_lambda_1 * E1[t - 1] + true_lambda_2 * E2[t - 1] + psi[t] * x[t] + omega_E1[t]
    E2[t] <- E1[t - 1]
    y[t] <- E1[t] + nu[t]
  }

  mod <- dlm(transfer_function = list(order = 2, xreg = x,
                                      discount_factor = 0.998))
  m <- c(0, 0, 0, 0, 0)
  C <- diag(x = c(100, 0.01, 0.5, 0.5, 10), nrow = 5)
  out <- fit(model = mod, y = y, m0 = m, C0 = C)

  expect_equal(out$filtered$m[3L, nobs] - true_lambda_1, 0.0081, tolerance = 1e-4)
  expect_equal(out$filtered$m[4L, nobs] - true_lambda_2, -0.0998, tolerance = 1e-4)

  # plot(y)
  # lines(out$filtered$f[, 1L], col = "blue")

})

test_that("TF(2) + level model with simulated data", {

  # Simulate the data
  nobs <- 1000
  npulse <- 0.25 * nobs
  sd_y <- 1
  sd_mu <- 0.05
  sd_E1 <- 0.0004
  sd_psi <- 0.002
  true_lambda_1 <- 0.8
  true_lambda_2 <- -0.2

  mu <- numeric(nobs)
  E1 <- numeric(nobs)
  E2 <- numeric(nobs)
  psi <- numeric(nobs)
  y <- numeric(nobs)
  E1[1L] <- 0
  psi[1L] <- 2.5
  mu[1L] <- 1

  # Covariate
  x <- numeric(nobs) # rpois(nobs, lambda = 1)
  impulse <- seq(1, nobs, by = nobs / npulse)
  x[impulse] <- 1
  prop.table(table(x))

  # Random errors
  set.seed(1111)
  nu <- rnorm(n = nobs, sd = sd_y)
  omega_psi <- rnorm(n = nobs, sd = sd_psi)
  omega_E1 <- rnorm(n = nobs, sd = sd_E1)
  omega_mu <- rnorm(n = nobs, sd = sd_mu)

  # First observation
  y[1L] <- mu[1L] + E1[1L] + nu[1L]

  for (t in seq_len(nobs)[-1]) {
    psi[t] <- psi[t - 1] + omega_psi[t]
    E1[t] <- true_lambda_1 * E1[t - 1] + true_lambda_2 * E2[t - 1] + psi[t] * x[t] + omega_E1[t]
    E2[t] <- E1[t - 1]
    mu[t] <- mu[t - 1] + omega_mu[t]
    y[t] <- mu[t] + E1[t] + nu[t]
  }

  mod <- dlm(polynomial = list(order = 1, discount_factor = 1),
             transfer_function = list(order = 2, xreg = x,
                                      discount_factor = 1))
  p <- nrow(mod$FF)
  or_tf <- 2L
  prior_m0 <- matrix(data = c(mean(mu), rep(0, p - 1L)))
  prior1_C0 <- diag(x = c(rep(10, 1L),
                          9, 0.001, # E_i
                          rep(0.5, or_tf), # lambda_i
                          10))
  prior10_C0 <- diag(x = c(rep(10, 1L),
                           9, 0.001, # E_i
                           rep(10, or_tf), # lambda_i
                           10))
  fit1 <- fit(model = mod, y = y, m0 = prior_m0, C0 = prior1_C0)
  fit10 <- fit(model = mod, y = y, m0 = prior_m0, C0 = prior10_C0)

  # out
  # expect_s3_class(fit1, "dlm.fit")

  # x11()
  # plot(mu)
  # plot(fit1$filtered$m[1L, ], col = 4)
  # lines(fit10$filtered$m[1L, ], col = 2)

  # x11()
  # plot(fit1$filtered$m[4L, ], col = 4, type = "l", ylim = c(-1.5, 1.5))
  # lines(fit10$filtered$m[4L, ], col = 2)
  # abline(h = true_lambda_1)

  # x11()
  # plot(fit1$filtered$m[5L, ], col = 4, type = "l", ylim = c(-1.5, 1))
  # lines(fit10$filtered$m[5L, ], col = 2)
  # abline(h = true_lambda_2)

  d1 <- extract(fit1, component = "state")
  d10 <- extract(fit10, component = "state")

  # dplot <- d1 |>
  #   dplyr::filter(parameter %in% c("lambda_1", "lambda_2")) |>
  #   dplyr::mutate(prior = "N[0, 1]") |>
  #   dplyr::bind_rows(
  #     d10 |>
  #       dplyr::filter(parameter %in% c("lambda_1", "lambda_2")) |>
  #       dplyr::mutate(prior = "N[0, 10]")
  #   ) |>
  #   dplyr::as_tibble() |>
  #   dplyr::mutate(true = ifelse(parameter == "lambda_1",
  #                               true_lambda_1,
  #                               true_lambda_2))
  # library(ggplot2)
  # data_pulses <- data.frame(t = impulse)
  # p <- ggplot(data = dplyr::filter(dplot, t > 10),
  #             aes(x = t, y = mean)) +
  #   facet_wrap(~parameter, ncol = 1, scales = "free_y") +
  #   geom_line(aes(col = prior)) +
  #   geom_hline(aes(yintercept = true)) +
  #   geom_ribbon(aes(ymin = ci_lower__95, ymax = ci_upper__95, fill = prior),
  #               alpha = 0.4) +
  #   geom_rug(data = dplyr::filter(data_pulses, t > 5),
  #            mapping =  aes(x = t, y = NA_real_))
  # p


  # plot(y)
  # lines(out$filtered$f[, 1L], col = "blue")
  # lines(out$filtered$m[1L, ], col = "red")
  # plot(mu)
  # lines(out$filtered$m[1L, ], col = "red")
  #
  # plot(out$filtered$m[4, ], col = "red"); abline(h = true_lambda_1)
  # plot(out$filtered$m[5, ], col = "red"); abline(h = true_lambda_2)


})
