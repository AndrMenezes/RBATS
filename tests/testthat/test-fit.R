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
  str(fitted_object)


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

  model$n_parms
  model$ar_order
  model$D
  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], 0, 0), ncol = 1),
                       C0 = diag(c(100, 2, 1)))
  str(fitted_object)
  # plot(y)
  # lines(fitted_object$filtered$f[, 1L], col = "blue")
  # lines(fitted_object$filtered$m[1L, ], col = "red")
  # fitted_object$filtered$m[3L, ]


  model <- dlm(
    polynomial = list(order = 1, discount_factor = 0.90),
    autoregressive = list(order = 2, discount_factor = 0.998)
  )

  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], 0, 0, 0, 0), ncol = 1),
                       C0 = diag(c(100, 2, 2, 1, 1)))

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
  out <- fit(model = mod, y = y, m = m0, C = C0)
  expect_equal(out$filtered$m[3L, nobs] - true_phi_1, -0.0792, tolerance = 1e-4)
  expect_equal(out$filtered$m[4L, nobs] - true_phi_2, 0.0597, tolerance = 1e-4)

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
             m = matrix(c(y[1L], rep(0, 4)), ncol = 1),
             C = diag(c(10, 2, 2, 1, 1), nrow = 5))

  expect_equal(out$filtered$m[4L, nobs] - true_phi_1, 0.0447, tolerance = 1e-4)
  expect_equal(out$filtered$m[5L, nobs] - true_phi_2, 0.0463, tolerance = 1e-4)

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
  out <- fit(model = mod, y = y, m = m0, C = C0)

  expect_equal(out$filtered$m[6L, nobs] - true_phi_1, -0.0149, tolerance = 1e-4)
  expect_equal(out$filtered$m[7L, nobs] - true_phi_2, 0.0334, tolerance = 1e-4)

  plot(out$filtered$m[6L, ]); abline(h = true_phi_1)
  lines(out$smoothed$ak[6L, ], col = "blue")
  plot(out$filtered$m[7L, ]); abline(h = true_phi_2)
  lines(out$smoothed$ak[7L, ], col = "blue")

})

