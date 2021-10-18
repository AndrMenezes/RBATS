test_that("autoplot for growth model works", {

  y <- c(Nile)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = 0.9),
    variance_law = list(type = "power", power = 1/2)
  )
  # Fitting the model
  fitted_model <- fit(object = model_object, y = y, prior_length = 10)

  # Predictive
  autoplot(fitted_model, what = "predictive")
  autoplot(fitted_model, what = "predictive", type = "filter")
  autoplot(fitted_model, what = "predictive", type = "filter", interval = FALSE)
  autoplot(fitted_model, what = "predictive", type = "filter", ts_geom = "line")
  autoplot(fitted_model, what = "predictive", type = "filter", ts_geom = "line",
           interval_geom = "line")

  autoplot(fitted_model, what = "predictive", type = "smooth")
  autoplot(fitted_model, what = "predictive", type = "smooth", interval = FALSE)
  autoplot(fitted_model, what = "predictive", type = "smooth", ts_geom = "line")

  autoplot(fitted_model, what = "predictive", type = "both")
  autoplot(fitted_model, what = "predictive", type = "both", ts_geom = "point",
           ts_size = 3)

  # Posterior
  autoplot(fitted_model, what = "posterior")
  autoplot(fitted_model, what = "posterior", type = "filter", interval_geom = "line")
  autoplot(fitted_model, what = "posterior", type = "filter", interval = FALSE)

  autoplot(fitted_model, what = "posterior", type = "smooth")

  autoplot(fitted_model, what = "posterior", type = "both")
})

test_that("autoplot for seasonal fourier model works", {
  y <- c(AirPassengers)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "free", period = 12, harmonics = 1:3),
    discount_factors = list(polynomial = c(0.90, 0.95), seasonal = 0.95)
  )

  # Fitting the model
  fitted_model <- fit(object = model_object, y = y, prior_length = 20)
  head(fitted_model$data_posterior_smooth, 10)
  tail(fitted_model$data_posterior_smooth)

  # Predictive
  autoplot(fitted_model, what = "predictive")
  autoplot(fitted_model, what = "predictive", type = "filter")
  autoplot(fitted_model, what = "predictive", type = "filter", interval = FALSE)
  autoplot(fitted_model, what = "predictive", type = "filter", ts_geom = "point")
  autoplot(fitted_model, what = "predictive", type = "filter", interval = FALSE,
           ts_geom = "point",  ts_size = 3)

  autoplot(fitted_model, what = "predictive", type = "smooth")
  autoplot(fitted_model, what = "predictive", type = "both")

  # Posterior
  autoplot(fitted_model, what = "posterior")
  autoplot(fitted_model, what = "posterior", type = "filter", interval_geom = "line")
  autoplot(fitted_model, what = "posterior", type = "filter", interval = FALSE)

  autoplot(fitted_model, what = "posterior", type = "smooth")
  autoplot(fitted_model, what = "posterior", type = "smooth", interval_geom = "line")
  autoplot(fitted_model, what = "posterior", type = "both")
})
