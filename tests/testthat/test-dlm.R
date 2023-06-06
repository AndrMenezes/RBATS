test_that("polynomial model works", {
  (m1 <- dlm(polynomial = list(order = 2, discount_factor = c(0.95, 0.98))))
  expect_equal(unname(diag(m1$D)), 1/c(0.95, 0.98))

  (m1_monitor <- dlm(
    polynomial = list(order = 2, discount_factor = c(0.95, 0.98)),
    monitor = list(
      execute = TRUE,
      bf_threshold = 0.135,
      location_shift = 4,
      scale_shift = 1,
      discount_factors = list(polynomial = c(0.20, 0.50)))
    ))
  expect_equal(unname(diag(m1_monitor$monitor$D)), 1/c(0.20, 0.50))

  (m2 <- dlm(polynomial = list(order = 3,
                               discount_factor = c(0.95, 0.98, 0.99))))
  expect_equal(unname(diag(m2$D)), 1/c(0.95, 0.98, 0.99))

  (m2_monitor <- dlm(polynomial = list(order = 3,
                                       discount_factor = c(0.95, 0.98, 0.99)),
                     monitor = list(
                       execute = TRUE,
                       bf_threshold = 0.135,
                       location_shift = 4,
                       scale_shift = 1,
                       discount_factors = list(polynomial = 0.20))))
  expect_equal(unname(diag(m2_monitor$monitor$D)), 1/c(0.20, 0.98, 0.99))


})

test_that("seasonal free model works",{


  p <- 3

  (m1 <- dlm(polynomial = list(order = 1, discount_factor = 0.95),
             seasonal = list(type = "free", period = p,
                             discount_factor = 0.98)))
  expect_equal(unname(diag(m1$D)), 1/c(0.95, rep(0.98, p - 1)))

  (m1_monitor <- dlm(polynomial = list(order = 1, discount_factor = 0.95),
             seasonal = list(type = "free", period = p,
                             discount_factor = 0.98),
             monitor = list(
               execute = TRUE,
               bf_threshold = 0.135,
               location_shift = 4,
               scale_shift = 1,
               discount_factors = list(polynomial = 0.10))))
  expect_equal(unname(diag(m1_monitor$monitor$D)), 1/c(0.10, rep(0.98, p - 1)))

  (m2 <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
             seasonal = list(type = "free", period = p,
                             discount_factor = 0.98)))
  expect_equal(unname(diag(m2$D)), 1/c(0.95, 0.95, rep(0.98, p - 1)))

  (m2_monitor <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
                     seasonal = list(type = "free", period = p,
                                     discount_factor = 0.98),
                     monitor = list(
                       execute = TRUE,
                       bf_threshold = 0.135,
                       location_shift = 4,
                       scale_shift = 1,
                       discount_factors = list(polynomial = c(0.10, 0.5)))))
  expect_equal(unname(diag(m2_monitor$monitor$D)), 1/c(0.10, 0.5, rep(0.98, p - 1)))

  (m2_monitor <- dlm(polynomial = list(order = 2, discount_factor = 0.95),
                     seasonal = list(type = "free", period = p,
                                     discount_factor = 0.98),
                     monitor = list(
                       execute = TRUE,
                       bf_threshold = 0.135,
                       location_shift = 4,
                       scale_shift = 1,
                       discount_factors = list(polynomial = c(0.10, 0.5),
                                               seasonal = 0.80))))
  expect_equal(unname(diag(m2_monitor$monitor$D)), 1/c(0.10, 0.5, rep(0.80, p - 1)))

})

test_that("seasonal fourier model works",{

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 7, harmonics = 1:3,
                    discount_factor = 0.98)
  ))

  (m2 <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 7, harmonics = 1:3,
                    discount_factor = 0.98)
  ))

  (m3 <- dlm(
    polynomial = list(order = 3, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98)
  ))

})

test_that("dynamic linear regression model works",{

  (m1 <- dlm(
    polynomial_order = 1,
    xreg = matrix(rnorm(100), ncol = 4)
  ))
  expect_equal(ncol(m1$GG), 5L)

  (m2 <- dlm(
    polynomial_order = 2,
    xreg = matrix(rnorm(100), ncol = 4)
  ))
  expect_equal(ncol(m2$GG), 6L)

  (m3 <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = 0.95, seasonal = 0.98, regressors = 0.96),
    seasonal = list(type = "fourier", period = 3, harmonics = 1),
    xreg = matrix(rnorm(100), ncol = 4)
  ))
  expect_equal(ncol(m3$GG), 8L)

})
