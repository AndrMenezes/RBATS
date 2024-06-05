rm(list = ls())
# devtools::load_all()
test_that("test forecast level and growth", {
  y <- c(Nile)

  # Define the model
  model <- dlm(polynomial = list(order = 1, discount_factor = 1))
  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1]), ncol = 1),
                       C0 = diag(10, 1))
  # Forecast
  f <- forecast(x = fitted_object, horizon = 10)
  expect_length(unique(f$mean), 1L)

  # Define the model
  model <- dlm(polynomial = list(order = 2, discount_factor = 0.95))
  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], 0), ncol = 1),
                       C0 = diag(10, 2))
  # Forecast
  f <- forecast(x = fitted_object, horizon = 20)
  # plot(y, xlim = c(0, max(f$t)))
  # lines(f$t, f$mean, col = "red")
  # lines(f$t, f$ci_lower__95, col = "red")
  # lines(f$t, f$ci_upper__95, col = "red")

})

test_that("test forecast seasonal", {
  y <- c(AirPassengers)

  # Define the model
  model <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
               seasonal = list(type = "fourier", period = 12, harmonics = 1:2,
                               discount_factor = 0.98))
  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], rep(0, 5L)), ncol = 1),
                       C0 = diag(10, 6))
  # Forecast
  f <- forecast(x = fitted_object, horizon = 24)
  # t <- seq_along(y)
  # rmv <- -(1:12)
  # plot(y, xlim = c(0, max(f$t)), ylim = range(y, f$ci_upper__80))
  # lines(t[rmv], fitted_object$filtered$f[rmv], col = "blue")
  # lines(f$t, f$mean, col = "red")
  # lines(f$t, f$ci_lower__95, col = "red", lty = 2)
  # lines(f$t, f$ci_upper__95, col = "red", lty = 2)

})


test_that("test forecast autoregressive", {
  y <- c(Nile)

  # Define the model
  model <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
               autoregressive = list(order = 1, discount_factor = 0.998))
  model$FF
  # Fitting the model
  fitted_object <- fit(model = model, y = y,
                       m0 = matrix(c(y[1], 0, 0,  0), ncol = 1),
                       C0 = diag(10, 4))
  class(fitted_object)
  f1 <- forecast(x = fitted_object, horizon = 4)
  f2 <- forecast(x = fitted_object, horizon = 4, discount = TRUE)
  expect_type(f1, "list")
  expect_type(f2, "list")
  expect_equal(f1$mean, f2$mean)
  expect_type(
    forecast(x = fitted_object, horizon = 4, state_parameters = TRUE),
    "list"
  )
  expect_type(
    forecast(x = fitted_object, horizon = 4, discount = TRUE,
             state_parameters = TRUE),
    "list"
  )


  # rmv <- -(1:8)
  # t <- seq_along(y)
  # t_h <- f1$t
  # x11()
  # plot(t[rmv], y[rmv])
  # lines(t[rmv], fitted_object$filtered$f[rmv], col = "blue")
  # lines(f1$t, f1$mean, col = "red")
  # lines(f1$t, f1$ci_lower__95, col = "red")
  # lines(f1$t, f1$ci_upper__95, col = "red")
  # lines(f1$t, f2$ci_lower__95, col = "red", lty = 2)
  # lines(f1$t, f2$ci_upper__95, col = "red", lty = 2)

})

