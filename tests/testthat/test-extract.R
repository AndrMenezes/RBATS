test_that("extract method for dlm objects", {

  y <- c(AirPassengers)
  m <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
           seasonal = list(type = "free", period = 12,
                           harmonics = 1:3, discount_factor = 0.98))
  x <- fit(model = m, y = y,
           m0 = matrix(c(110, rep(0, 12)), ncol = 1),
           C0 = diag(1000, 13))

  e_s <- extract(x = x, prob_interval = c(0.05, 0.20), distribution = "smooth",
                 component = "state")
  head(e_s)
  sliced_e_s <- e_s[e_s$parameter == "sum_seasonality", ]
  head(sliced_e_s)
  # plot(sliced_e_s$mean, type = "l", col = "blue")
  # lines(sliced_e_s$ci_lower__95, col = "blue", lty = 2)
  # lines(sliced_e_s$ci_upper__95, col = "blue", lty = 2)

})
