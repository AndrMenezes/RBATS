test_that("filter and smooth functions for growth model for Nile data works", {

  y <- c(Nile)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = c(0.90, 0.95))
  )

  # Fitting the model
  fitted_object <- fit(object = model_object,  y = y, prior_length = 10)

  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])

})

test_that("filter and smooth functions with missing value for Nile data works", {
  y <- c(Nile)
  y[c(10, 30)] <- NA_real_
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = c(0.90, 0.95))
  )

  # Fitting the model
  fitted_object <- fit(object = model_object,  y = y, prior_length = 10)

  # x11(); autoplot(fitted_object, what = "predictive")
  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])
})


test_that("filter and smooth for seasonal growth model for AirPassengers data works", {
  y <- c(AirPassengers)

  # Seasonal fourier (sin and cos)
  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "fourier", period = 12, harmonics = c(1, 2, 3)),
    discount_factors = list(polynomial = c(0.90), seasonal = 0.98)
  )
  fitted_object <- fit(object = model_object, y = y, prior_length = 40)
  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])

  # Seasonal free
  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "free", period = 12),
    discount_factors = list(polynomial = 0.90, seasonal = 0.98)
  )
  fitted_object <- fit(object = model_object, y = y, prior_length = 40)
  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])

})

test_that("filter and smooth for seasonal growth model for us_retail_employment data works", {
  data(us_retail_employment)
  y <- us_retail_employment$employed

  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "fourier", period = 12, harmonics = c(1, 2, 3)),
    discount_factors = list(polynomial = 0.90, seasonal = 0.98)
  )
  fitted_object <- fit(object = model_object, y = y, prior_length = 40)
  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])

  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "free", period = 12),
    discount_factors = list(polynomial = 0.90, seasonal = 0.98)
  )
  fitted_object <- fit(object = model_object, y = y, prior_length = 40)
  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])

})

test_that("filter and smooth functions for regression model for Nile data works", {
  y <- c(Nile)
  # Fake covariate
  X <- matrix(rnorm(2 * length(y)), ncol = 2)
  # Define the model
  model_object <- dlm(
    polynomial_order = 1, xreg = X,
    discount_factors = list(polynomial = c(0.90, 0.95), regressors = 1)
  )

  # Fitting the model
  fitted_object <- fit(object = model_object,  y = y, prior_length = 10)
  expect_true(fitted_object[["filtered_parameters"]])
  expect_true(fitted_object[["smooth_parameters"]])
})

test_that("fit method for dgegm object for simulated data works", {
  set.seed(6669)
  n <- 60
  theta_1 <- 50
  theta_2 <- 45
  theta_3 <- 0.95
  t <- seq.int(1, n)
  mu <- theta_1 - theta_2 * theta_3^t
  y <- mu + rnorm(n, sd = 0.05 * mu)

  lambdas <- c(-1, -1/2, 0, 1/2, 1)
  out <- lapply(lambdas, function(l) {
    # Define the model
    model_object <- dgegm(lambda = -1, discount_factors = c(0.90, 0.90, 0.95),
                          variance_law = list(type = "power", p = 2))
    # Fitting the model
    fit(object = model_object, y = y)
  })
  expect_length(out, length(lambdas))
})

test_that("fit method for dgegm object for AIDS data works", {
  y <- aids_brasil$cases
  lambdas <- c(-1, -1/2, 0, 1/2, 1)
  out <- lapply(lambdas, function(l) {
    # Define the model
    model_object <- dgegm(lambda = l, discount_factors = c(0.90, 0.90, 0.95),
                          variance_law = list(type = "poisson"))
    # Fitting the model
    fit(object = model_object, y = y)
  })
  expect_length(out, length(lambdas))
})
