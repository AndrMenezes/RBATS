test_that("plot for growth model works", {

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
  plot(fitted_model, what = "predictive")
  plot(fitted_model, what = "predictive", type = "smooth")
  plot(fitted_model, what = "predictive", type = "filter", interval = FALSE)
  plot(fitted_model, what = "predictive", type = "both")

  # Posterior
  plot(fitted_model, what = "posterior")
  plot(fitted_model, what = "posterior", exclude_prior = TRUE)
  plot(fitted_model, what = "posterior", type = "smooth")
  plot(fitted_model, what = "posterior", type = "filter", interval = FALSE)
  plot(fitted_model, what = "posterior", type = "smooth", interval = FALSE)

  plot(fitted_model, what = "posterior", type = "smooth", parm = "level")
  plot(fitted_model, what = "posterior", type = "smooth", parm = "growth")
  plot(fitted_model, what = "posterior", type = "both", parm = "growth")
  plot(fitted_model, what = "posterior", type = "both", parm = "level")
  plot(fitted_model, what = "posterior", type = "both")
  plot(fitted_model, what = "posterior", type = "both", exclude_prior = TRUE)
})

test_that("plot for seasonal fourier model works", {
  y <- c(AirPassengers)
  # Define the model
  model_object <- dlm(
    polynomial_order = 2,
    seasonal = list(type = "free", period = 12, harmonics = 1:3),
    discount_factors = list(polynomial = c(0.90, 0.95), seasonal = 0.95)
  )

  # Fitting the model
  fitted_model <- fit(object = model_object, y = y, prior_length = 20)
  # Predictive
  plot(fitted_model, what = "predictive")
  plot(fitted_model, what = "predictive", type = "smooth")
  plot(fitted_model, what = "predictive", type = "filter", interval = FALSE)
  plot(fitted_model, what = "predictive", type = "both")

  # Posterior
  plot(fitted_model, what = "posterior")
  plot(fitted_model, what = "posterior", type = "smooth")

  plot(fitted_model, what = "posterior", type = "smooth", parm = "growth")
  plot(fitted_model, what = "posterior", type = "smooth", parm = "seas_11")
  plot(fitted_model, what = "posterior", type = "smooth",
       parm = c("level", "growth", "seas_1"))
  plot(fitted_model, what = "posterior", type = "both", parm = "growth")
  plot(fitted_model, what = "posterior", type = "both", parm = "level")
  plot(fitted_model, what = "posterior", type = "both",
       parm = c("level", "growth", "seas_1"))
  plot(fitted_model, what = "posterior", type = "both",
       parm = c("level", "growth", "seas_1"), exclude_prior = TRUE)
})
