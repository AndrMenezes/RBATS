test_that("forecast method works for level and growth models", {

  y <- c(Nile)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = 0.9),
    variance_law = list(type = "power", power = 1/2)
  )
  # Fitting the model
  fitted_model <- fit(object = model_object, y = y, prior_length = 10)
  model_object <- fitted_model$model
  # Forecasting
  fcast_pred <- forecast(model_object, horizon = 10)
  fcast_pred_parms <- forecast(model_object, horizon = 10, state_parameters = TRUE)

  expect_s3_class(fcast_pred$predictive, "data.frame")
  expect_null(fcast_pred$state_parameters)
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
})

test_that("forecast method works for seasonal models", {
  y <- c(AirPassengers)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "fourier", period = 12, harmonics = 1:3),
    discount_factors = list(polynomial = c(0.90, 0.95), seasonal = 0.95)
  )
  # Fitting the model
  fitted_model <- fit(object = model_object, y = y, prior_length = 20)
  model_object <- fitted_model$model
  # Forecasting
  fcast_pred <- forecast(model_object, horizon = 10)
  fcast_pred_parms <- forecast(model_object, horizon = 10, state_parameters = TRUE)

  expect_s3_class(fcast_pred$predictive, "data.frame")
  expect_null(fcast_pred$state_parameters)
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
})

test_that("forecast method works for regression models", {
  y <- c(Nile)
  # Fake covariate
  X <- matrix(rnorm(2 * length(y)), ncol = 2)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2, xreg = X,
    seasonal = list(type = "fourier", period = 3, harmonics = 2),
    discount_factors = list(polynomial = c(0.90, 0.95), regressors = 0.95,
                            seasonal = 0.98)
  )
  # Fitting the model
  fitted_model <- fit(object = model_object, y = y, prior_length = 20)
  model_object <- fitted_model$model
  # Forecasting
  h <- 10
  X_future <- matrix(rnorm(2 * h), ncol = 2)
  fcast_pred <- forecast(model_object, horizon = h, xreg = X_future)
  fcast_pred_parms <- forecast(model_object, horizon = h, xreg = X_future,
                               state_parameters = TRUE)

  expect_s3_class(fcast_pred$predictive, "data.frame")
  expect_null(fcast_pred$state_parameters)
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
})

test_that("forecast method works for dgegm class", {
  y <- aids_brasil$cases

  # Define the model
  model_object <- dgegm(lambda = 0, discount_factors = c(0.90, 0.90, 0.998),
                        variance_law = list(type = "poisson"))
  # Fitting the model
  fitted_model <- fit(object = model_object, y = y)
  model_object <- fitted_model$model

  # Forecasting
  h <- 5
  fcast_pred <- forecast(model_object, horizon = h)
  fcast_pred_parms <- forecast(model_object, horizon = h, state_parameters = TRUE)

  expect_s3_class(fcast_pred$predictive, "data.frame")
  expect_null(fcast_pred$state_parameters)
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
  expect_s3_class(fcast_pred_parms$state_parameters, "data.frame")
})
