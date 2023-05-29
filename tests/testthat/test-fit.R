test_that("filter and smooth functions for growth model for Nile data works", {
  y <- c(Nile)

  # Define the model
  model <- dlm(
    polynomial = list(order = 1, discount_factor = 0.90)
  )

  # Fitting the model
  fitted_object <- fit_cpp.dlm(model = model, y = y,
                               a = matrix(c(y[1]), ncol = 1),
                               R = diag(100, 1))

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
  fitted_object <- fit_cpp.dlm(model = model, y = y,
                               a = matrix(c(y[1]), ncol = 1),
                               R = diag(100, 1))
  str(fitted_object)
  fitted_object$model$posterior$C


  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
})



test_that("filter and smooth for seasonal growth model for AirPassengers data works", {
  y <- c(AirPassengers)

  # Seasonal fourier (sin and cos)
  model_object <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 12, harmonics = c(1, 2, 3),
                    discount_factor = 0.98)
  )
  fitted_object <- fit_cpp.dlm(model = model_object, y = y,
                               a = matrix(c(y[1], rep(0, 7)), ncol = 1),
                               R = diag(100, 8))
  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")

  # Seasonal free
  model_object <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "free", period = 12,
                    discount_factor = 0.99)
  )
  fitted_object <- fit_cpp.dlm(model = model_object, y = y,
                               a = matrix(c(y[1], rep(0, 12)), ncol = 1),
                               R = diag(100, 13))
  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")
})

test_that("filter and smooth for seasonal growth model for us_retail_employment data works", {
  data(us_retail_employment)
  y <- us_retail_employment$employed

  model_object <- dlm(
    polynomial = list(order = 2, discount_factor = 0.90),
    seasonal = list(type = "fourier", period = 12, harmonics = c(1, 2, 3),
                    discount_factor = 0.98)
  )
  fitted_object <- fit_cpp.dlm(model = model_object, y = y,
                               a = matrix(c(y[1], rep(0, 7)), ncol = 1),
                               R = diag(100, 8))
  expect_s3_class(fitted_object, "dlm.fit")
  expect_s3_class(fitted_object$model, "dlm")


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
    model_object <- dgegm(lambda = l, discount_factors = c(0.95, 0.95, 0.998))
    # Fitting the model
    fit(object = model_object, y = y)
  })
  plot(y)
  lines(out[[1]]$data_predictive$mean)

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
