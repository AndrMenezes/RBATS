test_that("polynomial model works", {
  (m1 <- dlm(polynomial_order = 2, discount_factors = list(polynomial = c(0.95, 0.98))))

  (m2 <- dlm(polynomial_order = 3,
             discount_factors = list(polynomial = c(0.95, 0.98, 0.99))))

  (m3 <- dlm(polynomial_order = 4,
             discount_factors = list(polynomial = c(0.95, 0.98, 0.99, 0.998))))

})

test_that("seasonal free model works",{

  (m1 <- dlm(polynomial_order = 1,
             discount_factors = list(polynomial = 0.95, seasonal = 0.98),
             seasonal = list(type = "free", period = 12)))

  (m2 <- dlm(polynomial_order = 2,
             discount_factors = list(polynomial = c(0.95, 0.98), seasonal = 0.98),
             seasonal = list(type = "free", period = 12)))

  (m3 <- dlm(polynomial_order = 1,
             discount_factors = list(polynomial = 0.95, seasonal = 0.98),
             seasonal = list(type = "free", period = 6)))
  (m4 <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = c(0.95, 0.98), seasonal = 0.98),
    seasonal = list(type = "free", period = 6)
  ))


})

test_that("seasonal fourier model works",{

  (m1 <- dlm(
    polynomial_order = 1,
    discount_factors = list(polynomial = 0.95, seasonal = 0.98),
    seasonal = list(type = "fourier", period = 12, harmonics = 1:2)
  ))

  (m2 <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = 0.95, seasonal = 0.98),
    seasonal = list(type = "fourier", period = 12, harmonics = 1:3)
  ))

  (m3 <- dlm(
    polynomial_order = 2,
    discount_factors = list(polynomial = 0.95, seasonal = 0.98),
    seasonal = list(type = "fourier", period = 3, harmonics = 1)
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
