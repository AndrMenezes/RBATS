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

test_that("seasonal free model works", {

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

test_that("seasonal fourier model works", {

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 7, harmonics = 1:3,
                    discount_factor = 0.98)
  ))

  expect_equal(nrow(m1$FF), 7L)

  (m2 <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 7, harmonics = 1:3,
                    discount_factor = 0.98)
  ))
  expect_equal(nrow(m2$FF), 8L)

  (m3 <- dlm(
    polynomial = list(order = 3, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98)
  ))
  expect_equal(nrow(m3$FF), 5L)

})

test_that("dynamic linear regression model works",{

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    regressor = list(xreg = matrix(runif(10), ncol = 2), discount_factor = 0.99)
  ))
  expect_equal(nrow(m1$GG), 3L)
  expect_equal(m1$i_polynomial, 1L)
  expect_equal(m1$i_regressor, c(2L, 3L))

  (m2 <- dlm(
    polynomial = list(order = 2, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98),
    regressor = list(xreg = matrix(runif(10), ncol = 2), discount_factor = 0.99)
  ))
  expect_equal(nrow(m2$GG), 6L)
  expect_equal(m2$i_polynomial, c(1L, 2L))
  expect_equal(m2$i_seasonal, c(3L, 4L))
  expect_equal(m2$i_regressor, c(5L, 6L))

})

test_that("dlm with trend + seasonal + regression + autoregressive",{

  (m1 <- dlm(
    polynomial = list(order = 1, discount_factor = 0.95),
    seasonal = list(type = "fourier", period = 3, harmonics = 1,
                    discount_factor = 0.98),
    regressor = list(xreg = matrix(runif(10), ncol = 2), discount_factor = 0.99),
    autoregressive = list(order = 2, discount_factor = 0.998)
  ))
  expect_equal(rownames(m1$FF), c("level", "period_3__fourier_cos_1",
                                  "period_3__fourier_sin_1", "regressor_1",
                                  "regressor_2", "xi_1", "xi_2", "phi_1",
                                  "phi_2"))
  expect_equal(nrow(m1$GG), 9L)
  expect_equal(m1$i_polynomial, 1L)
  expect_equal(m1$i_seasonal, c(2L, 3L))
  expect_equal(m1$i_regressor, c(4L, 5L))
  expect_equal(m1$i_autoregressive, 6L:9L)

})
